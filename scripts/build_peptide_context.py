"""Freeze independent GENCODE/parent-peptide evidence for terminal lysines.

Run with --proteins pointing to the benchmark's vaccine_genes.proteins.fa.
Never use Exacto reconstructions to decide which reference residues are scored.
"""

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from pipeline.config import GENCODE_PROTEINS_URL, RESULTS_DIR
from pipeline.evaluate import expected_change


def build_context(catalogue, proteins_path):
    proteins = defaultdict(list)
    header, sequence = None, []
    for line in proteins_path.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if header:
                fields = header.split("|")
                proteins[fields[6]].append((fields[0], "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line)
    entries = []
    for variant in catalogue["variants"]:
        expected = expected_change(variant)
        peptides = variant.get("published_vaccine_peptides", [])
        for peptide in peptides:
            sequence = peptide["sequence"]
            leading = len(sequence) - len(sequence.lstrip("K"))
            trailing = len(sequence) - len(sequence.rstrip("K"))
            if not (leading or trailing) or not peptide.get("in_vaccines"):
                continue
            start = leading if leading <= 4 else 0
            end = len(sequence) - trailing if trailing <= 4 else len(sequence)
            core = sequence[start:end]
            matches = []
            if len(core) >= 6:
                for protein_id, reference in proteins[variant["gene"]]:
                    for offset in range(len(reference) - len(core) + 1):
                        chunk = reference[offset:offset + len(core)]
                        differences = [(a, b) for a, b in zip(chunk, core) if a != b]
                        if differences and not (
                            expected["kind"] == "missense"
                            and differences == [(expected["ref_aa"], expected["alt_aa"])]
                        ):
                            continue
                        matches.append({
                            "protein_id": protein_id, "core_start": offset,
                            "reference_core": chunk,
                            "n_flank": reference[max(0, offset - start):offset],
                            "c_flank": reference[offset + len(core):offset + len(core) + len(sequence) - end],
                        })
            parents = []
            for parent in peptides:
                context = parent["sequence"]
                if len(context) <= len(sequence) or not parent.get("in_vaccines"):
                    continue
                offset = context.find(sequence)
                if offset >= 0:
                    parents.append({"sequence": context, "offset": offset,
                                    "in_vaccines": parent["in_vaccines"]})
            entries.append({
                "variant_id": variant["variant_id"], "sequence": sequence,
                "in_vaccines": peptide["in_vaccines"],
                "reference_matches": matches, "parent_matches": parents,
            })
    return {
        "input_id": catalogue["source"]["input_id"],
        "source_url": GENCODE_PROTEINS_URL,
        "protein_subset_sha256": hashlib.sha256(proteins_path.read_bytes()).hexdigest(),
        "policy": "Same-gene reference core matches are exact except for the annotated missense substitution. Coordinates are zero-based. Parent matches retain the complete recorded peptide.",
        "entries": entries,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("data/vaccine_peptide_context.json"))
    args = parser.parse_args()
    catalogue = json.loads((RESULTS_DIR / "vaccine_variants.json").read_text())
    args.output.write_text(json.dumps(build_context(catalogue, args.proteins), indent=2) + "\n")
