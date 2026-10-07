"""Build the vaccine panel from frozen, corrected osteosarc inputs, offline."""

from __future__ import annotations

import json
from collections import defaultdict
from typing import Any

from .config import RESULTS_DIR, TIMEPOINT_ORDER
from .osteosarc_inputs import (
    INPUTS_DIR,
    digest,
    frozen_inputs,
    input_id,
    manifest,
    require_version,
)
from .vaccine_peptides import write_vaccine_peptide_subset

# Shared by the reference downloader; catalogue construction does no HTTP I/O.
USER_AGENT = "DoesExactoWorkYet/1.0 (+https://github.com/pirl-unc/DoesExactoWorkYet)"


def load_vaccine_epitopes() -> dict[tuple[str, int], list[dict]]:
    path = INPUTS_DIR / "epitopes.json"
    if digest(path) != manifest()["epitopes_sha256"]:
        raise ValueError("Frozen osteosarc epitope checksum changed")
    return {
        (row["chrom"], row["pos"]): row["epitopes"]
        for row in json.loads(path.read_text())
    }


def _number(value, cast):
    return cast(value) if value not in (None, "", "NA") else None


def _assay_support(rows: list[dict]) -> list[dict]:
    support = []
    for row in rows:
        support.append(
            {
                **{
                    key: row.get(key) or None
                    for key in (
                        "sample_label",
                        "data_source",
                        "pipeline",
                        "assay_type",
                        "timepoint",
                        "tissue",
                        "sample_date",
                        "bam_file",
                        "corrections",
                    )
                },
                **{
                    key: _number(row.get(key), int)
                    for key in ("ref_reads", "alt_reads", "total_reads")
                },
                "vaf": _number(row.get("vaf"), float),
            }
        )
    return sorted(
        support, key=lambda row: (row["assay_type"] or "", row["timepoint"] or "")
    )


def _ont_expectation(support: list[dict]) -> dict:
    rows = {
        row["timepoint"]: {
            key: row[key] for key in ("ref_reads", "alt_reads", "total_reads", "vaf")
        }
        for row in support
        if row["assay_type"] == "scRNA_ONT" and row["tissue"] == "tumor"
    }
    return {timepoint: rows.get(timepoint) for timepoint in TIMEPOINT_ORDER}


def _locus(chrom, pos) -> tuple[str, int]:
    return str(chrom).removeprefix("chr").replace("MT", "M"), int(pos)


def _vaccine_membership(variants: list[dict], vaccine_rows: list[dict]):
    """Keep unresolved assertions without claiming a resolved allele join."""
    from osteosarc.corpus import vaccine_membership
    from osteosarc.errors import IntegrityError

    ready_loci = {
        _locus(*variant["alleles"][0][:2])
        for variant in variants
        if variant["status"] == "ready"
    }
    unresolved_loci = defaultdict(list)
    for variant in variants:
        if variant["status"] == "ready":
            continue
        annotation = variant["annotations"]
        location = annotation.get("index", {}).get("location", "")
        chrom, separator, pos = location.partition(":")
        if not separator:
            source = annotation.get("source_record", {})
            chrom, pos = source.get("chr"), source.get("pos")
        if chrom and pos:
            unresolved_loci[_locus(chrom, pos)].append(variant["id"])

    ready_numbers, unresolved_rows = [], []
    for number, row in enumerate(vaccine_rows):
        locus = _locus(row["chrom"], row["pos"])
        candidates = unresolved_loci.get(locus, [])
        if locus not in ready_loci and candidates:
            if len(candidates) != 1:
                raise IntegrityError(
                    f"Vaccine row {number} ({row['gene']}) joins "
                    f"{len(candidates)} unresolved catalogue entries"
                )
            unresolved_rows.append((number, candidates[0]))
        else:
            # Unknown loci and ambiguous ready alleles still fail the strict join.
            ready_numbers.append(number)

    memberships, joined = vaccine_membership(
        variants, [vaccine_rows[number] for number in ready_numbers]
    )
    overlap = {
        ready_numbers[row["row"]]: dict(row, row=ready_numbers[row["row"]])
        for row in joined
    }
    for entry in memberships.values():
        entry["overlap_rows"] = [
            ready_numbers[number] for number in entry["overlap_rows"]
        ]
    for number, variant_id in unresolved_rows:
        entry = memberships[variant_id]
        entry["overlap_rows"].append(number)
        entry["overlap_join_status"] = "unresolved_allele"
        overlap[number] = dict(vaccine_rows[number], variant_id=variant_id, row=number)
    for variant_id in {variant_id for _, variant_id in unresolved_rows}:
        entry = memberships[variant_id]
        entry["overlap_rows"].sort()
        names = {
            name
            for number in entry["overlap_rows"]
            for name, used in vaccine_rows[number]["vaccines"].items()
            if used
        }
        entry["overlap_union"] = sorted(names)
        entry["vaccines_union"] = sorted(
            names | set(entry["source_variants"]) | set(entry["parsed_overlap"])
        )
        entry["included"] = bool(
            entry["vaccines_union"] or (entry["index_count"] or 0) > 0
        )
        entry["membership_disagreement"] = set(entry["source_variants"]) != names
    return memberships, [overlap[number] for number in range(len(vaccine_rows))]


def build_variant_records(inputs: dict | None = None) -> dict[str, Any]:
    # Use osteosarc's allele-aware membership join, including assertions that
    # exist only in the catalogue/index and preserving duplicate TECPR1 rows.
    require_version()
    inputs = frozen_inputs() if inputs is None else inputs
    memberships, overlap = _vaccine_membership(
        inputs["variants"], inputs["vaccine_rows"]
    )
    epitopes = load_vaccine_epitopes()
    counts = defaultdict(list)
    for row in inputs["counts"]:
        counts[row["variant_id"]].append(row)
    variants = []
    for variant in inputs["variants"]:
        membership = memberships[variant["id"]]
        if not membership["included"]:
            continue
        annotation = variant["annotations"]
        source = annotation.get("source_record", {})
        rows = [overlap[index] for index in membership["overlap_rows"]]
        allele = variant["alleles"][0] if variant["status"] == "ready" else None
        chrom, pos, ref, alt = allele or (None, None, None, None)
        support = _assay_support(counts[variant["id"]])
        tested = [row for row in rows if row.get("elispot_tested")]
        elispot = (
            tested[0] if tested else (rows[0] if rows else source.get("validation", {}))
        )
        variants.append(
            {
                "variant_id": variant["id"],
                "gene": variant["gene"],
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "assembly": variant["assembly"],
                "allele_status": variant["status"],
                "consequence": annotation.get("consequence"),
                "protein_change": annotation.get("protein_change")
                or (rows[0]["mutation"] if rows else None),
                "variant_type": "snv"
                if allele and len(ref) == len(alt) == 1
                else "indel",
                "vaccines": membership["vaccines_union"],
                "vaccine_membership": membership,
                "vaccine_label": rows[0]["mutation"]
                if rows
                else annotation.get("protein_change"),
                "impact": rows[0].get("impact") if rows else None,
                "peptide_classes": [row["note"] for row in rows if row.get("note")],
                "elispot": {
                    "tested": bool(elispot.get("elispot_tested")),
                    "status": elispot.get("elispot_status", "not_tested"),
                    "response": elispot.get("elispot_response"),
                },
                "vaf_trend": rows[0].get("vafs", []) if rows else [],
                "vaccine_epitopes": epitopes.get((chrom, pos), []),
                "published_vaccine_peptides": source.get("vaccine_peptides", []),
                "corrections": list(annotation.get("corrections", [])),
                "count_corrections": list(annotation.get("count_corrections", [])),
                "assay_support": support,
                "ont_expectation": _ont_expectation(support),
            }
        )
    variants.sort(
        key=lambda row: (row["chrom"] or "", row["pos"] or 0, row["variant_id"])
    )
    vaccine_names = sorted(
        {name for variant in variants for name in variant["vaccines"]}
    )
    return {
        "source": {
            "repository": "https://github.com/iskandr/osteosarc",
            "input_id": input_id(),
            "snapshot": inputs["snapshot"],
            "osteosarc_version": manifest()["osteosarc_version"],
            "membership_policy": manifest()["panel_policy"],
            "provenance": inputs["provenance"],
        },
        "vaccine_names": vaccine_names,
        "vaccine_set_sizes": {
            name: sum(name in v["vaccines"] for v in variants) for name in vaccine_names
        },
        "n_peptide_entries": len(overlap),
        "n_variants": len(variants),
        "n_ready": sum(v["allele_status"] == "ready" for v in variants),
        "variants": variants,
    }


def main() -> None:
    payload = build_variant_records()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "vaccine_variants.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")
    subset = write_vaccine_peptide_subset(payload, RESULTS_DIR)
    print(
        f"{payload['n_variants']} osteosarc vaccine variants ({payload['n_ready']} ready), "
        f"snapshot {payload['source']['snapshot']['name']} -> {out}"
    )
    print(
        f"{subset['n_variants']} targets with {subset['n_peptides']} recorded vaccine "
        f"peptides ({subset['n_unique_sequences']} distinct sequences) "
        f"-> {RESULTS_DIR / 'vaccine_peptide_subset.json'} and .fasta"
    )


if __name__ == "__main__":
    main()
