"""0.5 table contracts, including allele and peptide provenance boundaries."""

import csv
import sys

import pytest

from pipeline.evaluate import (
    integrated_pairs,
    peptides_by_rna_call,
    proteoforms_by_rna_call,
    rna_calls_by_variant,
    split_ids,
)
from pipeline.run_exacto import (
    Runner,
    StepFailed,
    write_somatic_tsv,
    write_transcript_support,
)


def table(path, rows):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def rna_row(id="1", **overrides):
    return dict(variant_id=id, assembled_transcript_name="read1",
                reference_transcript_id="ENST1", chromosome_1="chr1", position_1="999",
                position_2="1001", variant_type="SNV", sequence="T", **overrides)


def test_modern_dna_and_rna_match_exact_allele(tmp_path):
    variant = {"variant_id": "v1", "gene": "TEST", "chrom": "chr1", "pos": 1000, "ref": "C", "alt": "T"}
    dna = tmp_path / "dna.tsv"
    write_somatic_tsv([variant], dna, transcript_cli=True)
    row = next(csv.DictReader(dna.open(), delimiter="\t"))
    assert row["variant_id"] == "1" and row["origin"] == "somatic"
    rna = table(tmp_path / "rna.tsv", [rna_row(), {**rna_row("2"), "sequence": "A"}])
    support = table(tmp_path / "support.tsv", [{"assembled_transcript_name": "read1", "read_names": "a;b"}])
    hits = rna_calls_by_variant(rna, [variant], support)["v1"]
    assert [hit["rna_variant_call_id"] for hit in hits] == ["1"]
    assert hits[0]["n_supporting_reads"] == 2
    pairs = table(tmp_path / "integrated.tsv", [{"dna_variant_id": "1", "rna_variant_id": "2"}])
    assert integrated_pairs(pairs) == {"1": {"2"}}
    assert split_ids("1;2,3") == {"1", "2", "3"}


@pytest.mark.parametrize("sequence,position_2,shift", [("", "1003", False), ("", "1002", True), ("AGGG", "1032", True)])
def test_event_only_indel_is_translated_and_frame_uses_net_change(tmp_path, sequence, position_2, shift):
    rna = table(tmp_path / "rna.tsv", [{**rna_row(), "variant_type": "INS" if sequence else "DEL",
                                       "sequence": sequence, "position_2": position_2}])
    rows = []
    for i, aa in enumerate("MACDEFGHIKLMNPQRSTVWY"):
        for codon in range(3):
            rows.append({"proteoform_id": "1", "assembled_transcript_name": "read1",
                             "amino_acid_index": str(i), "amino_acid": aa, "codon_index": str(codon),
                             "assembled_transcript_variant_id": "", "dna_variant_ids": "9",
                             "preceding_event_assembled_transcript_variant_id": "1" if i == 4 and codon == 2 else "",
                             "is_amino_acid_variant": "true" if i == 4 else "false"})
    path = table(tmp_path / "nucleotides.tsv", rows)
    found = proteoforms_by_rna_call(path, {"1", "9"}, rna=rna)
    assert set(found) == {"1"}, "integrated DNA IDs cannot supply RNA evidence"
    assert found["1"][0]["mutant_residue_indices"] == [4]
    assert found["1"][0]["mutant_residues"] == "E"
    assert found["1"][0]["frameshift"] is shift


def test_peptides_need_position_specific_rna_evidence(tmp_path):
    peptides = table(tmp_path / "peptides.tsv", [{
        "proteoform_id": "1", "assembled_transcript_variant_ids": "1;2",
        "mutant_peptide_sequence": "MACDEFGH", "k": "8", "amino_acid_index_start": "0", "amino_acid_index_end": "7",
    }])
    forms = {
        "1": [{"peptide_id": 1, "mutant_residue_indices": [4], "frameshift": False}],
        "2": [{"peptide_id": 1, "mutant_residue_indices": [15], "frameshift": False}],
    }
    assert set(peptides_by_rna_call(peptides, {"1", "2"}, forms)) == {"1"}
    # Deletions may only appear as preceding events and be absent from the
    # proteoform-wide ID list. Nucleotide provenance still anchors the peptide.
    forms["3"] = [{"peptide_id": 1, "mutant_residue_indices": [7], "frameshift": False}]
    assert set(peptides_by_rna_call(peptides, {"1", "2", "3"}, forms)) == {"1", "3"}


def test_transcript_support_preserves_sequence_without_inventing_assembly_reads(tmp_path):
    pytest.importorskip("pysam")
    fasta = tmp_path / "sequences.fa"
    fasta.write_text(">read1 description\nACGT\nTT\n")
    dest = tmp_path / "support.tsv"
    write_transcript_support(fasta, dest, reads=True)
    row = next(csv.DictReader(dest.open(), delimiter="\t"))
    assert row == {"assembled_transcript_name": "read1", "sequence": "ACGTTT", "read_names": "read1"}
    write_transcript_support(fasta, dest, reads=False)
    assert next(csv.DictReader(dest.open(), delimiter="\t"))["read_names"] == ""
    nexus = table(tmp_path / "nexus.tsv", [
        {"transcript_id": "read1", "read_name": name} for name in ["r2", "r1", "r2"]
    ])
    write_transcript_support(fasta, dest, reads=False, nexus_reads=nexus)
    assert next(csv.DictReader(dest.open(), delimiter="\t"))["read_names"] == "r1;r2"


def test_worker_error_is_in_recorded_failure(tmp_path):
    workers = tmp_path / "workers" / "1"
    workers.mkdir(parents=True)
    (workers / "stderr.txt").write_text("FileNotFoundError: spoa\n")
    runner = Runner(tmp_path / "logs", 1)
    with pytest.raises(StepFailed):
        runner.run("worker", [sys.executable, "-c", "raise SystemExit(1)"], worker_log_dir=workers.parent)
    assert "FileNotFoundError: spoa" in runner.steps[0]["log_tail"]
