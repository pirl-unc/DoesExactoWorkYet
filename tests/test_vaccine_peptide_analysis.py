"""Guard exact-sequence recovery against unrelated loci and incomplete archives."""

import csv
import json

import pytest

from pipeline import build_site
from pipeline.score_vaccine_peptides import analyze, peptide_positions, protein_record
from pipeline.vaccine_peptides import (
    PEPTIDE_SCORING_POLICY,
    peptide_scoring_region,
    results_fingerprint,
)


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


@pytest.mark.parametrize("n_tag", ["", "K", "KK", "KKK", "KKKK"])
@pytest.mark.parametrize("c_tag", ["", "K", "KK", "KKK", "KKKK"])
def test_terminal_tags_are_excluded_but_core_and_allele_link_are_required(tmp_path, n_tag, c_tag):
    subset, previous = archive(tmp_path)
    peptide = subset["variants"][0]["published_vaccine_peptides"][0]
    peptide["sequence"] = n_tag + "AKQ" + c_tag
    result = analyze(subset, previous, tmp_path, run_url="run", version="0.5.0a1")
    scored = result["variants"][0]["published_vaccine_peptides"][0]
    assert scored["sequence"] == peptide["sequence"]
    assert scored["scoring_region"]["sequence"] == "AKQ"  # internal K remains required
    assert scored["scoring_region"]["start"] == len(n_tag)
    assert scored["scoring_region"]["end"] == len(n_tag) + 3
    assert scored["matches"][0]["amino_acid_starts"] == [2]  # core position, not tag position
    assert scored["matches"][0]["rna_call_id"] == "42"
    assert result["summary"]["n_peptide_entries_matched"] == 1
    assert result["analysis"]["peptide_scoring_policy"] == PEPTIDE_SCORING_POLICY


@pytest.mark.parametrize("sequence", ["KANQKK", "KAKQTTGKK", "KQ", "KKKK"])
def test_tags_do_not_rescue_core_mismatches_missing_core_or_unlinked_core(tmp_path, sequence):
    subset, previous = archive(tmp_path)
    subset["variants"][0]["published_vaccine_peptides"][0]["sequence"] = sequence
    result = analyze(subset, previous, tmp_path, run_url="run", version="0.5.0a1")
    assert result["summary"]["n_peptide_entries_matched"] == 0


def test_internal_and_longer_lysine_runs_are_not_partially_stripped():
    assert peptide_scoring_region("AKKKQA")["sequence"] == "AKKKQA"
    assert peptide_scoring_region("KKKKKAKQKKKKK")["sequence"] == "KKKKKAKQKKKKK"
    assert peptide_scoring_region("KKKK")["sequence"] == ""


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
    report = {"source": {"input_id": "panel"}, "analysis": {
        "results_sha256": results_fingerprint(results), "peptide_scoring_policy": PEPTIDE_SCORING_POLICY,
    }}
    (tmp_path / "vaccine_peptide_analysis.json").write_text(json.dumps(report))
    monkeypatch.setattr(build_site, "RESULTS_DIR", tmp_path)
    catalogue = {"source": {"input_id": "panel"}}
    assert build_site.current_peptide_report(catalogue, results) == report
    assert build_site.current_peptide_report(catalogue, {"runs": []}) is None
    assert build_site.current_peptide_report({"source": {"input_id": "other"}}, results) is None
    assert build_site.current_peptide_report(catalogue, None) is None
    del report["analysis"]["peptide_scoring_policy"]
    (tmp_path / "vaccine_peptide_analysis.json").write_text(json.dumps(report))
    assert build_site.current_peptide_report(catalogue, results) is None


def test_unmatched_proteins_and_rna_without_translation_remain_reported(tmp_path):
    subset, previous = archive(tmp_path)
    subset["variants"][0]["published_vaccine_peptides"][0]["sequence"] = "DIFFERENT"
    result = analyze(subset, previous, tmp_path, run_url="run", version="0.5.0a1")
    variant = result["variants"][0]
    assert not variant["published_vaccine_peptides"][0]["matches"]
    candidate = variant["reconstructed_candidates"][0]
    assert result["reconstructions"][candidate["reconstruction_id"]] == "MAKQTT"
    assert candidate["variant_amino_acid_start"] == 3
    assert variant["exacto_runs"][0]["n_variant_rnas"] == 1
    (tmp_path / "proteins.fasta").write_text("")
    result = analyze(subset, previous, tmp_path, run_url="run", version="0.5.0a1")
    assert result["variants"][0]["reconstructed_candidates"] == []
    assert result["variants"][0]["exacto_runs"][0]["n_variant_rnas"] == 1


def test_source_rna_counts_match_the_actual_bam_and_do_not_pool_libraries():
    report = {"variants": [{"variant_id": "TEST"}]}
    rows = [
        {"tissue": "tumor", "assay_type": "scRNA_ONT", "bam_file": "IPISRC044_T2_sclrs_ONT_dedup.bam", "alt_reads": 4, "total_reads": 20},
        {"tissue": "tumor", "assay_type": "RNA", "bam_file": "other-T2.bam", "alt_reads": 10, "total_reads": 40},
        {"tissue": "blood", "assay_type": "RNA", "bam_file": "normal.bam", "alt_reads": 1, "total_reads": 10},
        {"tissue": "tumor", "assay_type": "WGS", "bam_file": "dna.bam", "alt_reads": 30, "total_reads": 60},
    ]
    catalogue = {"variants": [{"variant_id": "TEST", "assay_support": rows}]}
    enriched = build_site.peptide_report_with_rna_evidence(report, catalogue)
    evidence = enriched["variants"][0]["source_rna_support"]
    assert len(evidence) == 2
    assert evidence[0]["benchmark_sample"] == "T2-ONT"
    assert evidence[0]["alt_reads"] == 4
    assert evidence[1]["benchmark_sample"] is None
    assert "source_rna_support" not in report["variants"][0]


def test_rna_support_tiers_use_maximum_without_pooling_or_imputing_missing_samples():
    rows = [
        {"benchmark_sample": "T1-ONT", "alt_reads": 1},
        {"benchmark_sample": "T2-ONT", "alt_reads": 1},
        {"benchmark_sample": "T2-ONT", "alt_reads": 1},
        {"benchmark_sample": None, "alt_reads": 100},
    ]
    tier = build_site.rna_support_tier(rows)
    assert tier["category"] == "single"
    assert tier["max_alt_reads"] == 1
    assert tier["by_sample"]["T2-ONT"]["category"] == "single"
    assert "T1-PacBio" in tier["missing_samples"]
    rows.append({"benchmark_sample": "T3-ONT", "alt_reads": 2})
    assert build_site.rna_support_tier(rows)["category"] == "multiple"
    assert build_site.rna_support_tier([])["category"] == "unknown"
    assert build_site.rna_support_tier([{"benchmark_sample": "T1-ONT", "alt_reads": None}])["category"] == "unknown"
    assert build_site.rna_support_tier([{"benchmark_sample": "T1-ONT", "alt_reads": 0}])["category"] == "zero"


def test_raw_rna_evidence_promotes_a_tier_without_double_counting_calls_or_assemblies():
    evidence = [{"benchmark_sample": "T2-ONT", "alt_reads": 1}]
    calls = [{"transcript_model_id": "read-a"}, {"transcript_model_id": "read-a"}]
    recovery = {"samples": {"T1-PacBio": {"arms": {
        "reads": {"rna_variant_calls": calls},
        "corrected": {"rna_variant_calls": [{"transcript_model_id": str(i)} for i in range(20)]},
    }}}}
    assert build_site.rna_support_tier(evidence, recovery)["category"] == "single"
    calls.append({"transcript_model_id": "read-b"})
    tier = build_site.rna_support_tier(evidence, recovery)
    assert tier["category"] == "multiple"
    assert tier["max_alt_reads"] == 2
    assert tier["by_sample"]["T1-PacBio"]["sid_alt_reads"] is None
    assert tier["by_sample"]["T1-PacBio"]["exacto_raw_rnas"] == 2
    assert tier["by_sample"]["T2-ONT"]["category"] == "single"
