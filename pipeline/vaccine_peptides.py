"""Export recorded vaccine peptides, retaining their variant and vaccine labels."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


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
