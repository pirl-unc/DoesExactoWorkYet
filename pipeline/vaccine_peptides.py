"""Export recorded vaccine peptides, retaining their variant and vaccine labels."""

from __future__ import annotations

import hashlib
import json
import subprocess
from functools import lru_cache
from pathlib import Path

PEPTIDE_SCORING_POLICY = "reference_supported_terminal_tags_v2"
CONTEXT_PATH = Path(__file__).resolve().parent.parent / "data/vaccine_peptide_context.json"


def context_fingerprint() -> str:
    return hashlib.sha256(CONTEXT_PATH.read_bytes()).hexdigest()


@lru_cache
def terminal_contexts() -> dict:
    return json.loads(CONTEXT_PATH.read_text())


def score_primary_windows(report: dict) -> dict:
    """Share window grouping/ranking with the browser; refresh it on every build."""
    result = subprocess.run(
        ["node", str(Path(__file__).resolve().parent.parent / "scripts/score_primary_windows.js")],
        input=json.dumps(report), text=True, capture_output=True, check=True,
    )
    return json.loads(result.stdout)


def peptide_context(variant_id: str, sequence: str, input_id: str) -> dict | None:
    contexts = terminal_contexts()
    if contexts["input_id"] != input_id:
        return None
    return next((entry for entry in contexts["entries"]
                 if entry["variant_id"] == variant_id and entry["sequence"] == sequence), None)


def peptide_scoring_region(sequence: str, context: dict | None = None) -> dict:
    """Exclude terminal lysines only when independent sequence context supports it.

    Coordinates are zero-based, end-exclusive within the recorded peptide.
    A terminal K alone is not evidence of a tag. Encoded mRNA sequences, native
    parent/reference lysines, conflicting contexts and unknowns remain scored.
    """
    leading = len(sequence) - len(sequence.lstrip("K"))
    trailing = len(sequence) - len(sequence.rstrip("K"))
    start, end = 0, len(sequence)
    context = context or {}
    matches = context.get("reference_matches", [])
    if "mRNA" not in context.get("in_vaccines", []) and matches:
        for terminus, count, key in (("N", leading, "n_flank"), ("C", trailing, "c_flank")):
            if not 1 <= count <= 4 or not all(len(match[key]) == count for match in matches):
                continue
            # A longer parent that includes this end protects it, even when a
            # different transcript/reference would imply an addition.
            supported = any(
                "mRNA" in parent["in_vaccines"] or (
                    parent["offset"] > len(parent["sequence"]) - len(parent["sequence"].lstrip("K"))
                    if terminus == "N" else
                    parent["offset"] + len(sequence) < len(parent["sequence"].rstrip("K"))
                ) for parent in context.get("parent_matches", [])
            )
            if supported:
                continue
            indices = range(count) if terminus == "N" else range(count - 1, -1, -1)
            for index in indices:
                if any(match[key][index] == "K" for match in matches):
                    break
                if terminus == "N":
                    start += 1
                else:
                    end -= 1
    tags = []
    if start:
        tags.append({"terminus": "N", "start": 0, "end": start, "sequence": sequence[:start]})
    if end < len(sequence):
        tags.append({"terminus": "C", "start": end, "end": len(sequence), "sequence": sequence[end:]})
    return {"sequence": sequence[start:end], "start": start, "end": end, "terminal_tags": tags}


def results_fingerprint(results: dict) -> str:
    """Bind a sequence reanalysis to the exact scored run, not just its panel."""
    encoded = json.dumps(results, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_vaccine_peptide_subset(catalogue: dict) -> dict:
    variants = []
    for variant in catalogue["variants"]:
        peptides = [
            {
                "peptide_id": f"{variant['variant_id']}.peptide-{index}",
                **peptide,
                "scoring_region": peptide_scoring_region(peptide["sequence"], peptide_context(
                    variant["variant_id"], peptide["sequence"], catalogue["source"]["input_id"])),
            }
            for index, peptide in enumerate(
                variant.get("published_vaccine_peptides", []), start=1
            )
            if peptide.get("sequence") and peptide.get("in_vaccines")
        ]
        if not peptides:
            continue
        variants.append(
            {
                **{
                    key: variant.get(key)
                    for key in (
                        "variant_id", "gene", "chrom", "pos", "ref", "alt",
                        "assembly", "allele_status", "protein_change", "consequence",
                    )
                },
                "vaccines": sorted(
                    {name for peptide in peptides for name in peptide["in_vaccines"]}
                ),
                "published_vaccine_peptides": peptides,
            }
        )
    peptides = [
        peptide
        for variant in variants
        for peptide in variant["published_vaccine_peptides"]
    ]
    return {
        "source": catalogue["source"],
        "selection_policy": (
            "Keep published_vaccine_peptides entries with a nonempty sequence and "
            "nonempty in_vaccines. Include a variant if it has at least one such "
            "entry. Preserve full peptides and mRNA minimal epitopes as recorded; "
            "do not include pVACtools predictions or peptides without vaccine "
            "membership. Membership does not establish expression or RNA recovery."
        ),
        "peptide_scoring_policy": PEPTIDE_SCORING_POLICY,
        "peptide_context_sha256": context_fingerprint(),
        "n_catalogue_variants": catalogue["n_variants"],
        "n_variants": len(variants),
        "n_peptides": len(peptides),
        "n_unique_sequences": len({peptide["sequence"] for peptide in peptides}),
        "variants": variants,
    }


def write_vaccine_peptide_subset(catalogue: dict, output_dir: Path) -> dict:
    subset = build_vaccine_peptide_subset(catalogue)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "vaccine_peptide_subset.json").write_text(
        json.dumps(subset, indent=2) + "\n"
    )
    with (output_dir / "vaccine_peptide_subset.fasta").open("w") as out:
        for variant in subset["variants"]:
            for peptide in variant["published_vaccine_peptides"]:
                out.write(
                    f">{peptide['peptide_id']} gene={variant['gene']} "
                    f"protein_change={variant['protein_change']} "
                    f"vaccines={json.dumps(peptide['in_vaccines'])} "
                    "is_mrna_minimal_epitope="
                    f"{str(peptide['is_mrna_minimal_epitope']).lower()}\n"
                )
                sequence = peptide["sequence"]
                for offset in range(0, len(sequence), 80):
                    out.write(sequence[offset:offset + 80] + "\n")
    return subset
