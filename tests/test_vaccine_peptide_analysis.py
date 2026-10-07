"""Guard exact-sequence recovery against unrelated loci and incomplete archives."""

import csv
import json

import pytest

from pipeline import build_site
from pipeline.score_vaccine_peptides import analyze, peptide_positions, protein_record
from pipeline.vaccine_peptides import results_fingerprint


def call(**overrides):
    return {"read_start": 6, "read_end": 6, "variant_type": "SNV", "frameshift": False, **overrides}


def test_repeated_peptide_must_overlap_the_target_codon():
    protein = protein_record("read|orf_0-23", "ACDQQACD")
    assert peptide_positions(protein, call(read_start=18, read_end=18), "ACD") == [5]
    assert peptide_positions(protein, call(), "ACDQ") == [0]
    assert peptide_positions(protein, call(), "ACNQ") == []


def test_deletion_uses_junction_and_frameshift_allows_downstream_peptides():
    protein = protein_record("read|orf_0-17", "MAKQTT")
    deletion = call(read_start=5, variant_type="DEL")
    assert peptide_positions(protein, deletion, "AK") == [1]
    assert peptide_positions(protein, deletion, "MA") == []
    assert peptide_positions(protein, deletion, "QTT") == []
    assert peptide_positions(protein, {**deletion, "frameshift": True}, "QTT") == [3]


def test_orf_that_does_not_cover_the_mutation_cannot_match():
    protein = protein_record("read|orf_9-26", "MAKQTT")
    assert peptide_positions(protein, call(), "AKQ") == []


@pytest.mark.parametrize("header", ["read", "read|orf_0-20"])
def test_incompatible_fasta_fails_instead_of_silently_losing_matches(header):
    with pytest.raises(ValueError):
        protein_record(header, "MAKQTT")


def archive(tmp_path, *, allele="C", position=103, transcript="target"):
    variant = {
        "variant_id": "TEST", "chrom": "chr1", "pos": 103, "ref": "A", "alt": "C",
        "published_vaccine_peptides": [{
            "peptide_id": "TEST.peptide-1", "sequence": "AKQ",
            "in_vaccines": ["mRNA"], "is_mrna_minimal_epitope": False,
        }],
    }
    subset = {"source": {"input_id": "frozen"}, "variants": [variant]}
    previous = {
        "input_id": "frozen",
        "variants": [{"variant_id": "TEST", "outcome": "peptide", "residue_confirmed": True}],
        "runs": [{"sample": "T1-ONT", "arm": "reads", "status": "ok", "reads_input_id": "reads"}],
    }
    run = {
        **previous["runs"][0], "input_id": "frozen", "exacto_table_schema": "0.5",
        "outputs": {"rna_variant_calls": "calls.tsv", "primary_structures_fasta": "proteins.fasta"},
    }
    (tmp_path / "run.json").write_text(json.dumps(run))
    (tmp_path / "proteins.fasta").write_text(f">{transcript}|orf_0-17\nMAKQTT\n")
    (tmp_path / "rna_vars").mkdir()
    row = {
        "variant_id": "42", "chromosome_1": "chr1", "position_1": position - 1,
        "position_2": position + 1, "variant_type": "SNV", "sequence": allele,
        "assembled_transcript_name": "target", "read_start": 6, "read_end": 6,
    }
    with (tmp_path / "rna_vars/calls.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    return subset, previous


@pytest.mark.parametrize("changes,expected", [
    ({}, 1), ({"allele": "G"}, 0), ({"position": 104}, 0), ({"transcript": "decoy"}, 0),
])
def test_sequence_match_requires_same_transcript_and_exact_genomic_allele(tmp_path, changes, expected):
    subset, previous = archive(tmp_path, **changes)
    result = analyze(subset, previous, tmp_path, run_url="run", version="0.5.0a1")
    assert result["summary"]["n_variants_any_peptide_matched"] == expected
    if expected:
        match = result["variants"][0]["published_vaccine_peptides"][0]["matches"][0]
        assert match["amino_acid_starts"] == [2]
        assert match["rna_call_id"] == "42"


@pytest.mark.parametrize("problem", ["missing_file", "missing_run", "different_inputs", "different_reads"])
def test_incomplete_or_mismatched_archive_is_not_scored_as_a_miss(tmp_path, problem):
    subset, previous = archive(tmp_path)
    if problem == "missing_file":
        (tmp_path / "proteins.fasta").unlink()
    elif problem == "missing_run":
        (tmp_path / "run.json").unlink()
    elif problem == "different_inputs":
        subset["source"]["input_id"] = "different"
    else:
        previous["runs"][0]["reads_input_id"] = "different"
    with pytest.raises((ValueError, FileNotFoundError)):
        analyze(subset, previous, tmp_path, run_url="run", version="0.5.0a1")


def test_site_hides_sequence_report_from_a_different_run_or_catalogue(tmp_path, monkeypatch):
    results = {"runs": [{"sample": "T1-ONT", "status": "ok"}]}
    report = {"source": {"input_id": "panel"}, "analysis": {"results_sha256": results_fingerprint(results)}}
    (tmp_path / "vaccine_peptide_analysis.json").write_text(json.dumps(report))
    monkeypatch.setattr(build_site, "RESULTS_DIR", tmp_path)
    catalogue = {"source": {"input_id": "panel"}}
    assert build_site.current_peptide_report(catalogue, results) == report
    assert build_site.current_peptide_report(catalogue, {"runs": []}) is None
    assert build_site.current_peptide_report({"source": {"input_id": "other"}}, results) is None
    assert build_site.current_peptide_report(catalogue, None) is None
