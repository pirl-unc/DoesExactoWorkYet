"""Small, offline runtime check for the installed toolchain (opt-in locally).

Run with DEWY_EXACTO_INTEGRATION=1 in environment.yml after installing Exacto.
The reference and reads are synthetic. The spades case starts with an existing
contig, exercising the assembly-to-Exacto boundary without testing SPAdes here.
Clustering, correction, isoform assembly, alignment, every Exacto command and
scoring execute real binaries; no patient data or downloads are needed.
"""

import gzip
import json
import os
import random

import pytest

from pipeline import run_exacto
from pipeline.config import SAMPLES_BY_NAME
from pipeline.evaluate import evaluate_arm

pytestmark = pytest.mark.skipif(
    os.environ.get("DEWY_EXACTO_INTEGRATION") != "1",
    reason="requires the installed Exacto/toolchain; exercised in compatibility CI",
)


def synthetic_inputs(folder):
    import pysam

    rng = random.Random(150)
    codons = dict(zip(
        ["GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA",
         "TTA", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"],
        "ACDEFGHIKLNPQRSTVWY",
    ))
    coding = ["ATG"] + [rng.choice(list(codons)) for _ in range(299)]
    coding[150] = "GCT"
    protein = "M" + "".join(codons[codon] for codon in coding[1:])
    transcript = "TAA" * 16 + "TA" + "".join(coding) + "TAA" + "TAA" * 16 + "TA"
    mutant = transcript[:501] + "T" + transcript[502:]
    genome = list("".join(rng.choices("ACGT", k=4000)))
    exons = [(1001, 1350), (1801, 2100), (2551, 2903)]
    offset = 0
    for start, end in exons:
        length = end - start + 1
        genome[start - 1:end] = transcript[offset:offset + length]
        genome[start - 3:start - 1] = "AG"
        genome[end:end + 2] = "GT"
        offset += length
    assert offset == len(transcript)
    reference = folder / "reference.fa"
    reference.write_text(">chr1\n" + "".join(genome) + "\n")
    pysam.faidx(str(reference))
    annotation = folder / "annotation.gtf"
    common = 'gene_id "ENSG00000000001.1"; gene_name "TEST"; gene_type "protein_coding"; level 2;'
    tx = common + ' transcript_id "ENST00000000001.1"; transcript_name "TEST-201"; transcript_type "protein_coding"; tag "basic";'
    rows = []

    def row(kind, start, end, attrs, phase="."):
        rows.append(f"chr1\tHAVANA\t{kind}\t{start}\t{end}\t.\t+\t{phase}\t{attrs}\n")

    row("gene", 1001, 2903, common)
    row("transcript", 1001, 2903, tx)
    for i, (start, end) in enumerate(exons, 1):
        attrs = tx + f' exon_number {i}; exon_id "ENSE0000000000{i}.1";'
        row("exon", start, end, attrs)
        row("CDS", max(start, 1051), min(end, 2850), attrs, "0")
    row("start_codon", 1051, 1053, tx, "0")
    row("stop_codon", 2851, 2853, tx, "0")
    annotation.write_text("".join(rows))
    proteins = folder / "proteins.fa"
    proteins.write_text(">TEST\n" + protein + "\n")
    reads = folder / "reads.fastq.gz"
    with gzip.open(reads, "wt") as sink:
        for i in range(12):
            sink.write(f"@read{i}\n{mutant}\n+\n{'I' * len(mutant)}\n")
    contig = folder / "contig.fa"
    contig.write_text(">contig\n" + mutant + "\n")
    variant = {"variant_id": "synthetic-A151V", "gene": "TEST", "chrom": "chr1",
                   "pos": 1952, "ref": "C", "alt": "T", "protein_change": "p.Ala151Val"}
    return reference, annotation, proteins, reads, contig, variant


@pytest.mark.parametrize("arm", ["reads", "spades", "corrected", "isonform"])
def test_installed_pipeline_recovers_synthetic_mutation(tmp_path, monkeypatch, arm):
    reference, annotation, proteins, reads, contig, variant = synthetic_inputs(tmp_path)
    monkeypatch.setattr(run_exacto, "EXACTO_DIR", tmp_path / "exacto")
    monkeypatch.setattr(run_exacto, "MASKED_FASTA", reference)
    monkeypatch.setattr(run_exacto, "GENE_PROTEINS", proteins)
    monkeypatch.setattr(run_exacto, "ANNOTATION_ARGS", [
        "--reference-gene-annotation-file", str(annotation),
        *run_exacto.ANNOTATION_ARGS[2:],
    ])
    monkeypatch.setattr(run_exacto, "reads_arm_fastq", lambda sample: reads)
    monkeypatch.setattr(run_exacto, "extraction_outputs", lambda sample: [reads])
    monkeypatch.setattr(run_exacto, "run_rnaspades", lambda *args: contig)
    monkeypatch.setattr(run_exacto, "reads_input_id", lambda stats: "synthetic")
    stats = tmp_path / "stats.json"
    stats.write_text(json.dumps({
        "input_id": run_exacto.input_id(),
        "files": {reads.name: run_exacto.digest(reads)},
        "n_spanning_reads": 12,
    }))
    monkeypatch.setattr(run_exacto, "stats_path", lambda sample: stats)
    sample = SAMPLES_BY_NAME["T2-ILMN" if arm == "spades" else "T1-ONT"]
    result = run_exacto.run_arm(sample, arm, [variant], threads=2)
    assert result["status"] == "ok", json.dumps(result, indent=2)
    assert "peptide_error" not in result, json.dumps(result, indent=2)
    score = evaluate_arm(result, [variant], {variant["variant_id"]: 12})[variant["variant_id"]]
    assert score["outcome"] == "peptide", score
    assert score["residue_confirmed"] is True, score
    if arm == "spades":
        assert result["counts"]["after_unspliced_filter"] == 1
