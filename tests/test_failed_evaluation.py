"""An absent mutation score must not turn into a coverage claim or a verdict."""

from __future__ import annotations

import json
import sys

import pytest

from pipeline import build_site, evaluate, run_exacto
from pipeline.config import SAMPLES_BY_NAME

VARIANTS = [
    {"variant_id": "A", "gene": "GENEA", "protein_change": "p.Ala1Thr"},
    {"variant_id": "B", "gene": "GENEB", "protein_change": "p.Ala2Thr"},
]


@pytest.fixture
def scores(tmp_path, monkeypatch):
    monkeypatch.setattr(evaluate, "SCORED_DIR", tmp_path / "scored")
    monkeypatch.setattr(evaluate, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(evaluate, "load_variants", lambda: VARIANTS)
    evaluate.SCORED_DIR.mkdir()

    def write(arm, variants, status="ok"):
        run = {"sample": "T1-ONT", "arm": arm, "status": status, "variants": variants}
        (evaluate.SCORED_DIR / f"T1-ONT.{arm}.json").write_text(json.dumps({
            "sample": "T1-ONT", "runs": [run],
        }))
    return write


def test_failed_only_merge_has_unscored_mutations_not_no_reads(scores):
    scores("corrected", {}, "failed")
    payload = evaluate.merge()
    assert payload["evaluation_status"] == "failed"
    assert payload["n_evaluated"] == payload["n_testable"] == payload["n_recovered"] == 0
    assert payload["outcome_counts"]["no_reads"] == 0
    assert payload["outcome_counts"]["not_run"] == 2
    assert {v["outcome"] for v in payload["variants"]} == {"not_run"}


def test_completed_but_uncovered_is_distinct_from_failure(scores):
    scores("reads", {v["variant_id"]: {"outcome": "no_reads"} for v in VARIANTS})
    payload = evaluate.merge()
    assert payload["evaluation_status"] == "no_coverage"
    assert payload["n_evaluated"] == 2
    assert payload["n_testable"] == 0
    assert payload["outcome_counts"]["no_reads"] == 2


def test_failed_arm_does_not_erase_successful_measurement(scores):
    scores("reads", {"A": {"outcome": "proteoform"}, "B": {"outcome": "no_call"}})
    scores("assembly", {}, "failed")
    payload = evaluate.merge()
    assert payload["evaluation_status"] == "ok"
    assert payload["n_evaluated"] == payload["n_testable"] == 2
    assert payload["n_recovered"] == 1


def test_missing_variant_is_excluded_from_the_denominator(scores):
    scores("reads", {"A": {"outcome": "no_call"}})
    payload = evaluate.merge()
    assert payload["n_testable"] == payload["n_evaluated"] == 1
    assert payload["outcome_counts"]["not_run"] == 1
    assert payload["outcome_counts"]["no_reads"] == 0


def test_failed_run_cannot_contribute_partial_or_stale_scores(scores):
    scores("reads", {"A": {"outcome": "proteoform"}}, "failed")
    payload = evaluate.merge()
    assert payload["evaluation_status"] == "failed"
    assert payload["n_recovered"] == 0


def test_bad_allele_is_recorded_before_any_external_tool_runs(tmp_path, monkeypatch):
    monkeypatch.setattr(run_exacto, "EXACTO_DIR", tmp_path / "exacto")
    stats = tmp_path / "extraction.json"
    stats.write_text("{}")
    monkeypatch.setattr(run_exacto, "stats_path", lambda sample: stats)
    monkeypatch.setattr(evaluate, "stats_path", lambda sample: stats)
    monkeypatch.setattr(evaluate, "SCORED_DIR", tmp_path / "scored")
    bad_allele = {**VARIANTS[0], "chrom": "chr1", "pos": 1, "ref": "A", "alt": "<DEL>"}
    monkeypatch.setattr(evaluate, "load_variants", lambda: [bad_allele])

    def unexpected_tool(*args, **kwargs):
        pytest.fail("Unsupported input must fail before running external tools")
    monkeypatch.setattr(run_exacto.Runner, "run", unexpected_tool)
    run = run_exacto.run_arm(SAMPLES_BY_NAME["T1-ONT"], "reads", [bad_allele], 1)
    assert run["status"] == "failed"
    assert "cannot encode GENEA" in run["error"]
    assert run["steps"] == []
    saved = tmp_path / "exacto" / "T1-ONT" / "reads" / "run.json"
    assert json.loads(saved.read_text()) == run

    graded = evaluate.score_samples(["T1-ONT"])
    assert graded[0]["variants"] == {}
    assert graded[0]["error"] == run["error"]
    assert (evaluate.SCORED_DIR / "T1-ONT.reads.json").exists()


def test_failed_pipeline_exits_nonzero(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["run_exacto", "--samples", "T1-ONT", "--arms", "reads"])
    monkeypatch.setattr(run_exacto, "load_variants", lambda: VARIANTS)
    monkeypatch.setattr(run_exacto, "run_arm", lambda *args: {
        "sample": "T1-ONT", "arm": "reads", "status": "failed", "error": "bad input",
    })
    with pytest.raises(SystemExit) as error:
        run_exacto.main()
    assert error.value.code == 1


def test_site_recognises_the_old_failed_zero_over_zero_payload(tmp_path, monkeypatch, scores):
    scores("corrected", {}, "failed")
    payload = evaluate.merge()
    payload.pop("n_evaluated")
    payload.pop("evaluation_status")
    payload["outcome_counts"] = {"no_reads": 2, "no_call": 0, "rna_only": 0, "proteoform": 0, "peptide": 0}
    for v in payload["variants"]:
        v["outcome"] = "no_reads"
    (tmp_path / "exacto_results.json").write_text(json.dumps(payload))
    (tmp_path / "vaccine_variants.json").write_text(json.dumps({
        "variants": VARIANTS, "n_variants": 2, "n_peptide_entries": 2,
        "source": {}, "vaccine_names": [], "vaccine_set_sizes": {},
    }))
    history_path = tmp_path / "history.json"
    history = [{"n_recovered": 1, "n_testable": 2, "exacto_version": "test"}]
    history_path.write_text(json.dumps(history))
    monkeypatch.setattr(build_site, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(build_site, "HISTORY_PATH", history_path)
    monkeypatch.setattr(build_site, "git", lambda *args: None)
    site = build_site.build_payload()
    assert site["summary"]["evaluation_status"] == "failed"
    assert site["summary"]["outcome_counts"]["no_reads"] == 0
    assert site["history"] == history
    assert json.loads(history_path.read_text()) == history


@pytest.mark.parametrize("status", ["failed", "no_coverage"])
def test_unmeasured_run_does_not_replace_last_history_measurement(tmp_path, monkeypatch, status):
    path = tmp_path / "history.json"
    previous = [{"n_testable": 37, "n_recovered": 22, "exacto_version": "0.4.6a1"}]
    path.write_text(json.dumps(previous))
    monkeypatch.setattr(build_site, "HISTORY_PATH", path)
    assert build_site.update_history({
        "n_testable": 0, "n_recovered": 0, "evaluation_status": status,
    }) == previous
    assert json.loads(path.read_text()) == previous
