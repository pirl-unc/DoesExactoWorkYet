"""Re-score archived Exacto 0.5 proteins against recorded vaccine peptides.

Uses complete protein FASTAs and exact-allele RNA calls, not the three candidate
examples retained in the compact benchmark JSON. Does not run Exacto again.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from .config import RESULTS_DIR, SAMPLES
from .evaluate import read_tsv
from .run_exacto import as_graph_operation
from .vaccine_peptides import (
    PEPTIDE_SCORING_POLICY,
    peptide_scoring_region,
    results_fingerprint,
)

ORF_HEADER = re.compile(r"^(.+)\|orf_(\d+)-(\d+)$")


def read_proteins(path: Path):
    name, sequence = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield protein_record(name, "".join(sequence))
                name, sequence = line[1:].strip().split()[0], []
            else:
                sequence.append(line.strip())
    if name is not None:
        yield protein_record(name, "".join(sequence))


def protein_record(name: str, sequence: str) -> dict:
    match = ORF_HEADER.fullmatch(name)
    if not match:
        raise ValueError(f"Unsupported Exacto protein header: {name}")
    transcript, start, end = match.groups()
    start, end = int(start), int(end)
    if end - start + 1 != len(sequence) * 3:
        raise ValueError(f"Protein length disagrees with ORF coordinates: {name}")
    return {
        "protein_id": name, "transcript": transcript,
        "orf_start": start, "orf_end": end, "sequence": sequence,
    }


def exact_calls(path: Path, variants: list[dict]) -> dict:
    wanted = {}
    for variant in variants:
        start, end, _, kind, sequence = as_graph_operation(variant)
        wanted[(variant["chrom"], start, end, kind, sequence.upper())] = variant
    calls = defaultdict(list)
    for row in read_tsv(path):
        key = (
            row["chromosome_1"], int(row["position_1"]), int(row["position_2"]),
            row["variant_type"], (row.get("sequence") or "").strip('"').upper(),
        )
        variant = wanted.get(key)
        if variant is not None:
            calls[row["assembled_transcript_name"]].append({
                "variant_id": variant["variant_id"],
                "rna_call_id": row["variant_id"],
                "read_start": int(row["read_start"]),
                "read_end": int(row["read_end"]),
                "variant_type": row["variant_type"],
                "frameshift": (len(variant["alt"]) - len(variant["ref"])) % 3 != 0,
            })
    return calls


def peptide_positions(protein: dict, call: dict, peptide: str) -> list[int]:
    """Exact matches overlapping the variant codon or its frameshifted tail.

    Exacto's deletion span brackets the junction; the event is assigned to the
    first base after it (read_end). ORF and RNA-call coordinates are zero-based,
    inclusive. Return zero-based amino-acid positions, including repeated hits.
    """
    start = call["read_end"] if call["variant_type"] == "DEL" else call["read_start"]
    end = call["read_end"]
    orf_start, orf_end = protein["orf_start"], protein["orf_end"]
    if not peptide or end < orf_start or start > orf_end:
        return []
    first = (max(start, orf_start) - orf_start) // 3
    last = (min(end, orf_end) - orf_start) // 3
    positions, offset = [], 0
    while (index := protein["sequence"].find(peptide, offset)) >= 0:
        if index + len(peptide) > first and (call["frameshift"] or index <= last):
            positions.append(index)
        offset = index + 1
    return positions


def summarize(variants: list[dict], run_keys: set[tuple[str, str]] | None = None) -> dict:
    def matched(peptide):
        return any(
            run_keys is None or (hit["sample"], hit["arm"]) in run_keys
            for hit in peptide["matches"]
        )

    peptides = [p for v in variants for p in v["published_vaccine_peptides"]]
    return {
        "n_variants": len(variants),
        "n_variants_any_peptide_matched": sum(
            any(matched(p) for p in v["published_vaccine_peptides"]) for v in variants
        ),
        "n_variants_non_minimal_peptide_matched": sum(
            any(matched(p) and not p["is_mrna_minimal_epitope"]
                for p in v["published_vaccine_peptides"]) for v in variants
        ),
        "n_variants_minimal_epitope_matched": sum(
            any(matched(p) and p["is_mrna_minimal_epitope"]
                for p in v["published_vaccine_peptides"]) for v in variants
        ),
        "n_variants_all_peptides_matched": sum(
            all(matched(p) for p in v["published_vaccine_peptides"]) for v in variants
        ),
        "n_peptide_entries": len(peptides),
        "n_peptide_entries_matched": sum(matched(p) for p in peptides),
        "n_unique_sequences": len({p["sequence"] for p in peptides}),
        "n_unique_sequences_matched": len({p["sequence"] for p in peptides if matched(p)}),
        "n_non_minimal_peptide_entries": sum(not p["is_mrna_minimal_epitope"] for p in peptides),
        "n_non_minimal_peptide_entries_matched": sum(
            matched(p) and not p["is_mrna_minimal_epitope"] for p in peptides
        ),
        "n_minimal_epitope_entries": sum(p["is_mrna_minimal_epitope"] for p in peptides),
        "n_minimal_epitope_entries_matched": sum(
            matched(p) and p["is_mrna_minimal_epitope"] for p in peptides
        ),
    }


def analyze(subset: dict, previous: dict, outputs_dir: Path, *, run_url: str, version: str) -> dict:
    if subset["source"]["input_id"] != previous["input_id"]:
        raise ValueError("Subset and benchmark use different frozen inputs")
    variants = json.loads(json.dumps(subset["variants"]))
    by_id = {v["variant_id"]: v for v in variants}
    previous_variants = {v["variant_id"]: v for v in previous["variants"]}
    for variant in variants:
        old = previous_variants[variant["variant_id"]]
        variant["candidate_outcome"] = old["outcome"]
        variant["residue_confirmed"] = old["residue_confirmed"]
        variant["reconstructed_candidates"] = []
        variant["exacto_runs"] = []
        for peptide in variant["published_vaccine_peptides"]:
            peptide["matches"] = []
            peptide["scoring_region"] = peptide_scoring_region(peptide["sequence"])

    completed = {(run["sample"], run["arm"]): run for run in previous["runs"] if run["status"] == "ok"}
    seen, artifacts, reconstructions = set(), [], {}
    for run_path in sorted(outputs_dir.rglob("run.json")):
        run = json.loads(run_path.read_text())
        key = run["sample"], run["arm"]
        if key not in completed:
            continue
        if key in seen:
            raise ValueError(f"Duplicate artifact for {key}")
        if (run["input_id"] != previous["input_id"]
                or run.get("reads_input_id") != completed[key].get("reads_input_id")
                or run["status"] != "ok" or run.get("exacto_table_schema") != "0.5"):
            raise ValueError(f"Artifact identity/status mismatch for {key}")
        directory = run_path.parent
        rna_path = directory / "rna_vars" / Path(run["outputs"]["rna_variant_calls"]).name
        fasta_path = directory / Path(run["outputs"]["primary_structures_fasta"]).name
        # Missing files must fail the analysis, not become negative results.
        calls = exact_calls(rna_path, variants)
        for variant in variants:
            variant_calls = [
                (transcript, call)
                for transcript, entries in calls.items() for call in entries
                if call["variant_id"] == variant["variant_id"]
            ]
            variant["exacto_runs"].append({
                "sample": key[0], "arm": key[1],
                "n_variant_rnas": len({transcript for transcript, _ in variant_calls}),
            })
        for protein in read_proteins(fasta_path):
            for call in calls.get(protein["transcript"], []):
                start = call["read_end"] if call["variant_type"] == "DEL" else call["read_start"]
                if call["read_end"] < protein["orf_start"] or start > protein["orf_end"]:
                    continue
                reconstruction_id = f"{key[0]}/{key[1]}/{protein['protein_id']}"
                reconstructions[reconstruction_id] = protein["sequence"]
                candidate = {
                    "sample": key[0], "arm": key[1],
                    "reconstruction_id": reconstruction_id,
                    "protein_id": protein["protein_id"], "rna_call_id": call["rna_call_id"],
                    "variant_amino_acid_start": (max(start, protein["orf_start"]) - protein["orf_start"]) // 3 + 1,
                }
                by_id[call["variant_id"]]["reconstructed_candidates"].append(candidate)
                for peptide in by_id[call["variant_id"]]["published_vaccine_peptides"]:
                    positions = peptide_positions(protein, call, peptide["scoring_region"]["sequence"])
                    if positions:
                        peptide["matches"].append({
                            "sample": key[0], "arm": key[1],
                            "reconstruction_id": reconstruction_id,
                            "protein_id": protein["protein_id"],
                            "rna_call_id": call["rna_call_id"],
                            "amino_acid_starts": [position + 1 for position in positions],
                        })
        artifacts.append({
            "sample": key[0], "arm": key[1], "reads_input_id": run["reads_input_id"],
            "files": {
                str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (run_path, rna_path, fasta_path)
            },
        })
        seen.add(key)
    if missing := completed.keys() - seen:
        raise ValueError(f"Missing completed-run artifacts: {sorted(missing)}")
    if not seen:
        raise ValueError("No completed Exacto runs available")

    summary = summarize(variants)
    summary.update({
        "candidate_outcome_counts": dict(Counter(v["candidate_outcome"] for v in variants)),
        "n_translated_candidate_variants": sum(v["candidate_outcome"] in ("peptide", "proteoform") for v in variants),
        "n_residue_checkable": sum(v["residue_confirmed"] is not None for v in variants),
        "n_residue_confirmed": sum(v["residue_confirmed"] is True for v in variants),
    })
    return {
        "source": subset["source"],
        "analysis": {
            "run_url": run_url, "exacto_version": version,
            "mode": "rescore_saved_outputs",
            "results_sha256": results_fingerprint(previous),
            "metric_label": "Vaccine sequence contained in reconstruction",
            "peptide_scoring_policy": PEPTIDE_SCORING_POLICY,
            "n_completed_methods": len(seen),
            "samples_available": sorted({sample for sample, _ in seen}),
            "samples_unavailable": sorted({s.name for s in SAMPLES} - {sample for sample, _ in seen}),
            "matching_policy": (
                "Complete terminal runs of 1–4 lysines at either end of every "
                "recorded vaccine peptide are treated as suspected solubility "
                "additions outside the ORF and excluded from scoring by benchmark "
                "convention, including minimal epitopes. Internal lysines and "
                "longer terminal runs are retained. The nonempty remaining sequence "
                "must be a contiguous substring "
                "of an exported Exacto protein sequence; the reconstruction may "
                "extend on either side of the vaccine sequence. It must be "
                "linked to an RNA call with the target's exact "
                "genomic locus and allele. The matched peptide must overlap the "
                "call's translated codon/junction or frameshifted tail. Counts "
                "are any-candidate unions over the available methods and samples; "
                "all-peptides matched can combine different candidates. "
                "Predicted pVACtools epitopes are excluded. Full recorded sequences "
                "are preserved; scoring_region uses zero-based, end-exclusive "
                "recorded-peptide coordinates. amino_acid_starts locates the "
                "scored region in the protein (one-based), not any excluded N-terminal tag."
            ),
            "artifacts": artifacts,
        },
        "summary": summary,
        "reconstructions": reconstructions,
        "by_method": {
            arm: summarize(variants, {key for key in seen if key[1] == arm})
            for arm in sorted({arm for _, arm in seen})
        },
        "by_sample": {
            sample: summarize(variants, {key for key in seen if key[0] == sample})
            for sample in sorted({sample for sample, _ in seen})
        },
        "variants": variants,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subset", type=Path, default=RESULTS_DIR / "vaccine_peptide_subset.json")
    parser.add_argument("--results", type=Path, required=True, help="Merged results JSON from the same run")
    parser.add_argument("--outputs-dir", type=Path, required=True, help="Extracted exacto-outputs-* artifacts")
    parser.add_argument("--run-url", required=True)
    parser.add_argument("--exacto-version", required=True)
    parser.add_argument("--output", type=Path, default=RESULTS_DIR / "vaccine_peptide_analysis.json")
    args = parser.parse_args()
    result = analyze(
        json.loads(args.subset.read_text()), json.loads(args.results.read_text()),
        args.outputs_dir, run_url=args.run_url, version=args.exacto_version,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"], indent=2))
    print(f"-> {args.output}")


if __name__ == "__main__":
    main()
