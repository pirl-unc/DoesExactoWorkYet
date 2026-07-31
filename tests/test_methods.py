"""Method definitions: what reaches a tool's command line, and what must not.

Written after every assembly-unspliced leg died in its first step.
`unspliced-filter` steers this harness, not the Nexus filter, but it lived in
`params` and run_exacto forwards every params entry as `--key value`. The filter
exits 2 on an unrecognised argument, so the arm failed before assembling
anything — and the CI job still went green, because run_exacto records an arm
failure rather than raising.
"""

from __future__ import annotations

import inspect

from pipeline import run_exacto
from pipeline.methods import METHODS, NEXUS_FILTER_DEFAULTS, methods_for

# Every flag nexus_filter_rnabloom2_transcripts actually accepts, from its
# argument parser.
NEXUS_FILTER_FLAGS = {
    "min-mapping-quality", "min-read-support", "min-fraction-match",
    "base-quality",
}


def test_assembly_params_are_all_flags_the_filter_accepts():
    for method in METHODS:
        if method.family != "assembly":
            continue
        unknown = set(method.params) - NEXUS_FILTER_FLAGS
        assert not unknown, (
            f"{method.name} would pass {sorted(unknown)} to "
            "nexus_filter_rnabloom2_transcripts, which exits 2 on unrecognised "
            "arguments — harness settings belong in controls, not params"
        )


def test_harness_settings_live_in_controls_not_params():
    unspliced = next(m for m in METHODS if m.name == "assembly-unspliced")
    isonform = next(m for m in METHODS if m.name == "isonform")

    assert unspliced.controls["unspliced-filter"] == "off"
    assert "unspliced-filter" not in unspliced.params
    assert isonform.controls["assemble"] == "isonform"
    assert "assemble" not in isonform.params


def test_the_canonical_method_uses_nexus_defaults_unchanged():
    canonical = next(m for m in METHODS if m.name == "assembly")
    assert canonical.params == NEXUS_FILTER_DEFAULTS
    assert not canonical.controls


def test_only_params_are_forwarded_to_the_command_line():
    source = inspect.getsource(run_exacto.run_arm)
    assert "method.params.items()" in source
    assert "method.controls.items()" not in source, (
        "controls steer the harness and must never reach a tool's argv"
    )


def test_sweeps_differ_from_the_canonical_method_in_exactly_one_setting():
    canonical = next(m for m in METHODS if m.name == "assembly")
    for name, key in (("assembly-ms1", "min-read-support"),
                      ("assembly-unspliced", "unspliced-filter")):
        method = next(m for m in METHODS if m.name == name)
        combined = {**method.params, **method.controls}
        base = {**canonical.params, **canonical.controls}
        differing = {k for k in combined | base.keys()
                     if combined.get(k) != base.get(k)}
        assert differing == {key}, f"{name} should vary only {key}, varies {differing}"


def test_short_reads_get_only_the_short_read_method():
    assert [m.name for m in methods_for("short")] == ["spades"]
    assert "spades" not in [m.name for m in methods_for("long")]
