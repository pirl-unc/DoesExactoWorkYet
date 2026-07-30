"""Pull the vaccine-neoantigen and somatic-variant tables from osteosarc.com.

Produces ``results/vaccine_variants.json``: one record per somatic variant that
made it into at least one of Sid's five personalised vaccines, joined against
the portal's per-assay read counts so we know, before running anything, which
variants the ONT data even covers.
"""

from __future__ import annotations

import csv
import io
import json
import urllib.request
from collections import defaultdict
from typing import Any

from .config import (
    PVACTOOLS_EPITOPE_URLS,
    RESULTS_DIR,
    TIMEPOINTS,
    VACCINE_OVERLAP_URL,
    VARIANT_VAFS_COLUMNS_URL,
    VARIANT_VAFS_URL,
)

# The portal is behind Cloudflare and 403s the default urllib agent.
USER_AGENT = "DoesExactoWorkYet/1.0 (+https://github.com/pirl-unc/DoesExactoWorkYet)"

# Read counts in variant_vafs_long.tsv are keyed by this assay label for the
# ONT long-read single-cell data — the same BAMs this repo runs Exacto on.
ONT_ASSAY = "scRNA_ONT"


def _get(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read()


def _get_text(url: str) -> str:
    return _get(url).decode("utf-8")


def _int_or_none(value: str) -> int | None:
    return int(value) if value not in ("", "NA") else None


def _float_or_none(value: str) -> float | None:
    return float(value) if value not in ("", "NA") else None


def load_vaf_table() -> list[dict[str, str]]:
    """Long-format per-assay read counts for every curated somatic variant."""
    text = _get_text(VARIANT_VAFS_URL)
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def load_vaccine_overlap() -> dict[str, Any]:
    return json.loads(_get_text(VACCINE_OVERLAP_URL))


def load_vaccine_epitopes() -> dict[tuple[str, int], list[dict[str, Any]]]:
    """The mutant epitope sequences pVACtools predicted, keyed by locus.

    These are what the vaccine designs were picked from, so a proteoform that
    contains one is a literal match to the peptide that was made — a stronger
    claim than "carries the right amino-acid change".

    pVACtools reports Start/Stop in the VCF's coordinates for SNVs but shifts by
    one for indels depending on the caller, so each row is matched against a
    small window of candidate positions rather than a single one.
    """
    epitopes: dict[tuple[str, int], dict[str, dict[str, Any]]] = {}
    for mhc_class, url in PVACTOOLS_EPITOPE_URLS.items():
        reader = csv.DictReader(io.StringIO(_get_text(url)), delimiter="\t")
        for row in reader:
            chrom = row["Chromosome"]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            candidates = {int(row["Start"]), int(row["Start"]) + 1, int(row["Stop"])}
            sequence = row["MT Epitope Seq"]
            if not sequence:
                continue
            for position in candidates:
                key = (chrom, position)
                bucket = epitopes.setdefault(key, {})
                entry = bucket.setdefault(
                    sequence,
                    {
                        "sequence": sequence,
                        "wild_type": row.get("WT Epitope Seq") or None,
                        "mhc_class": mhc_class,
                        "gene": row.get("Gene Name"),
                        "hgvsp": row.get("HGVSp") or None,
                        "alleles": set(),
                    },
                )
                entry["alleles"].add(row.get("HLA Allele", ""))

    resolved: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for key, bucket in epitopes.items():
        rows = []
        for entry in bucket.values():
            rows.append({**entry, "alleles": sorted(a for a in entry["alleles"] if a)})
        rows.sort(key=lambda item: (len(item["sequence"]), item["sequence"]))
        resolved[key] = rows
    return resolved


def _allele_info(rows: list[dict[str, str]]) -> dict[str, Any]:
    """Locus-level fields that are constant across a variant's assay rows."""
    first = rows[0]
    return {
        "variant_id": first["variant_id"],
        "gene": first["gene"],
        "chrom": first["chrom"],
        "pos": int(first["pos"]),
        "ref": first["ref"],
        "alt": first["alt"],
        "variant_type": first["variant_type"],
        "protein_change": first["protein_change"] or None,
        "consequence": first["consequence"] or None,
    }


def _assay_support(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Every assay/timepoint the portal genotyped this locus in."""
    support = []
    for row in rows:
        support.append(
            {
                "sample_label": row["sample_label"],
                "data_source": row["data_source"],
                "pipeline": row["pipeline"] or None,
                "assay_type": row["assay_type"],
                "timepoint": row["timepoint"] or None,
                "tissue": row["tissue"],
                "sample_date": row["sample_date"] or None,
                "ref_reads": _int_or_none(row["ref_reads"]),
                "alt_reads": _int_or_none(row["alt_reads"]),
                "total_reads": _int_or_none(row["total_reads"]),
                "vaf": _float_or_none(row["vaf"]),
            }
        )
    support.sort(key=lambda item: (item["assay_type"], item["timepoint"] or ""))
    return support


def _ont_expectation(support: list[dict[str, Any]]) -> dict[str, Any]:
    """What the portal's own genotyping of the ONT BAMs says.

    This is the yardstick for the Exacto run: a variant with zero ONT alt reads
    at a timepoint cannot be recovered from that timepoint by anyone, so
    counting it as an Exacto failure would be unfair.
    """
    by_timepoint = {}
    for entry in support:
        if entry["assay_type"] != ONT_ASSAY or entry["timepoint"] is None:
            continue
        by_timepoint[entry["timepoint"]] = {
            "ref_reads": entry["ref_reads"],
            "alt_reads": entry["alt_reads"],
            "total_reads": entry["total_reads"],
            "vaf": entry["vaf"],
        }
    return {
        timepoint.name: by_timepoint.get(timepoint.name) for timepoint in TIMEPOINTS
    }


def build_variant_records() -> dict[str, Any]:
    """Merge the vaccine list with the per-assay read counts."""
    overlap = load_vaccine_overlap()
    vaf_rows = load_vaf_table()
    epitopes = load_vaccine_epitopes()

    by_locus: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in vaf_rows:
        by_locus[(row["chrom"], row["pos"])].append(row)

    # A single genomic variant can appear more than once in the vaccine table —
    # TECPR1's frameshift is carried as both an MHC I and an MHC II peptide —
    # so merge on locus and take the union of the vaccines.
    merged: dict[tuple[str, str], dict[str, Any]] = {}
    order: list[tuple[str, str]] = []

    for mutation in overlap["mutations"]:
        key = (mutation["chrom"], str(mutation["pos"]))
        rows = by_locus.get(key)
        if not rows:
            raise SystemExit(
                f"{mutation['gene']} {mutation['chrom']}:{mutation['pos']} is in the "
                "vaccine table but absent from variant_vafs_long.tsv — the portal's "
                "two tables have drifted apart."
            )

        vaccines = sorted(name for name, used in mutation["vaccines"].items() if used)

        if key not in merged:
            order.append(key)
            support = _assay_support(rows)
            merged[key] = {
                **_allele_info(rows),
                "vaccine_label": mutation["mutation"],
                "impact": mutation["impact"],
                "vaccines": vaccines,
                "peptide_classes": [],
                "elispot": {
                    "tested": mutation["elispot_tested"],
                    "status": mutation["elispot_status"],
                    "response": mutation["elispot_response"],
                },
                "vaf_trend": mutation["vafs"],
                "vaccine_epitopes": epitopes.get(
                    (mutation["chrom"], int(mutation["pos"])), []
                ),
                "assay_support": support,
                "ont_expectation": _ont_expectation(support),
            }
        else:
            record = merged[key]
            record["vaccines"] = sorted(set(record["vaccines"]) | set(vaccines))

        if mutation.get("note"):
            merged[key]["peptide_classes"].append(mutation["note"])

    variants = [merged[key] for key in order]
    variants.sort(key=lambda v: (v["chrom"], v["pos"]))

    return {
        "source": {
            "vaccine_overlap_url": VACCINE_OVERLAP_URL,
            "variant_vafs_url": VARIANT_VAFS_URL,
            "variant_vafs_columns_url": VARIANT_VAFS_COLUMNS_URL,
            "pvactools_epitope_urls": PVACTOOLS_EPITOPE_URLS,
        },
        "vaccine_names": overlap["vaccine_names"],
        "vaccine_set_sizes": overlap["set_sizes"],
        "n_peptide_entries": len(overlap["mutations"]),
        "n_variants": len(variants),
        "variants": variants,
    }


def main() -> None:
    payload = build_variant_records()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "vaccine_variants.json"
    out_path.write_text(json.dumps(payload, indent=2) + "\n")

    covered = sum(
        1
        for variant in payload["variants"]
        if any(
            (entry or {}).get("alt_reads", 0)
            for entry in variant["ont_expectation"].values()
        )
    )
    print(
        f"{payload['n_variants']} vaccine variants "
        f"({payload['n_peptide_entries']} peptide entries) -> {out_path}"
    )
    print(f"{covered} have >=1 ONT alt read at some timepoint per osteosarc.com")
    with_epitopes = sum(
        1 for variant in payload["variants"] if variant["vaccine_epitopes"]
    )
    print(f"{with_epitopes} have published pVACtools epitope sequences to match against")


if __name__ == "__main__":
    main()
