"""Old timepoint-keyed results must still render after the switch to samples.

A run that was already in flight when the run unit changed writes ``T1`` where
the pipeline now writes ``T1-ONT``. Those results are good; only their labels
are stale, and the site has to keep showing them rather than quietly reporting
that nothing ran.
"""

from __future__ import annotations

import json

from pipeline import evaluate
from pipeline.build_site import migrate_name, migrate_payload


def test_legacy_timepoint_maps_to_the_ont_sample():
    # ONT was the only platform a timepoint-keyed run could have used, so this
    # is exact rather than a guess.
    assert migrate_name("T1") == "T1-ONT"
    assert migrate_name("T3") == "T3-ONT"


def test_current_sample_names_pass_through_untouched():
    assert migrate_name("T1-ONT") == "T1-ONT"
    assert migrate_name("T1-PacBio") == "T1-PacBio"


def test_migrates_a_whole_legacy_payload():
    payload = migrate_payload(
        {
            "runs": [{"timepoint": "T2", "arm": "assembly", "status": "ok"}],
            "extraction": {"T2": {"n_reads": 10}},
            "variants": [
                {
                    "variant_id": "GENE-chr1-1",
                    "timepoints": {"T2": {"outcome": "proteoform", "arms": {}}},
                }
            ],
        }
    )

    assert payload["runs"][0]["sample"] == "T2-ONT"
    assert payload["runs"][0]["platform"] == "ONT"
    assert payload["extraction"] == {"T2-ONT": {"n_reads": 10}}
    variant = payload["variants"][0]
    assert variant["samples"]["T2-ONT"]["outcome"] == "proteoform"
    # The stale key is removed, not left alongside to be read by accident.
    assert "timepoints" not in variant


def test_leaves_a_current_payload_alone():
    payload = migrate_payload(
        {
            "runs": [{"sample": "T1-PacBio", "platform": "PacBio", "arm": "reads"}],
            "variants": [{"samples": {"T1-PacBio": {"outcome": "no_call"}}}],
        }
    )

    assert payload["runs"][0]["sample"] == "T1-PacBio"
    assert payload["runs"][0]["platform"] == "PacBio"
    assert payload["variants"][0]["samples"]["T1-PacBio"]["outcome"] == "no_call"


# --------------------------------------------------------------------------
# Scored files are per sample *and* arm, because CI runs them as separate jobs
# --------------------------------------------------------------------------


def _scored(sample, arm, outcome):
    return {
        "sample": sample,
        "arm": arm,
        "extraction": {"n_reads": 10},
        "runs": [
            {
                "sample": sample,
                "arm": arm,
                "status": "ok",
                "variants": {"GENE-chr1-1": {"outcome": outcome}},
            }
        ],
    }


def test_merge_keeps_both_arms_of_one_sample(tmp_path, monkeypatch):
    """Two legs of the same sample must not overwrite each other.

    Splitting the CI matrix by arm means T1-ONT/reads and T1-ONT/assembly are
    separate jobs writing separate artifacts. Before they were named per arm
    both wrote results/scored/T1-ONT.json, and download-artifact's
    merge-multiple would silently keep whichever landed second.
    """
    monkeypatch.setattr(evaluate, "SCORED_DIR", tmp_path)
    monkeypatch.setattr(evaluate, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(
        evaluate,
        "load_variants",
        lambda: [
            {
                "variant_id": "GENE-chr1-1",
                "gene": "GENE",
                "protein_change": "p.Ala1Thr",
                "consequence": "missense_variant",
            }
        ],
    )

    (tmp_path / "T1-ONT.reads.json").write_text(
        json.dumps(_scored("T1-ONT", "reads", "rna_only"))
    )
    (tmp_path / "T1-ONT.assembly.json").write_text(
        json.dumps(_scored("T1-ONT", "assembly", "proteoform"))
    )

    payload = evaluate.merge()

    arms = payload["variants"][0]["samples"]["T1-ONT"]["arms"]
    assert sorted(arms) == ["assembly", "reads"]
    # The verdict is the better of the two, not whichever file sorted last.
    assert payload["variants"][0]["samples"]["T1-ONT"]["outcome"] == "proteoform"
    assert payload["n_recovered"] == 1
