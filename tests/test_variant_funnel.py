"""The per-sample funnel: where mutations stop being recoverable.

Written after a wrong number reached the site. Changing the counting loop from
``.items()`` to ``.values()`` left the body referencing ``variant_id``, and it
did not raise — Python leaks the loop variable, so an earlier loop's last value
was still bound and the allele stage silently counted 0 instead of 24. A wrong
number rather than a traceback, caught only because the output looked
implausible. These tests make that class of failure loud.
"""

from __future__ import annotations

from pipeline.paths import variant_funnel


def _run(sample, method, platform, timepoint, variants, status="ok"):
    return {
        "sample": sample, "arm": method, "platform": platform,
        "timepoint": timepoint, "status": status,
        "method": {"label": method, "family": method},
        "variants": variants,
    }


def _entry(outcome="no_call", spanning=10, rna_calls=0, residue=None, **extra):
    return {
        "outcome": outcome,
        "spanning_reads": spanning,
        "rna_variant_calls": [{"rna_variant_call_id": str(i)}
                              for i in range(rna_calls)],
        "residue_confirmed": residue,
        **extra,
    }


PORTAL = [{
    "variant_id": "A", "ont_expectation": {"T1": {"alt_reads": 12, "vaf": 0.2}},
}, {
    "variant_id": "B", "ont_expectation": {"T1": {"alt_reads": 0, "vaf": 0.0}},
}]


def _stage(row, key):
    return next(s for s in row["stages"] if s["key"] == key)


def test_allele_stage_counts_variants_not_zero():
    """The regression. Two variants, one with portal alt reads."""
    payload = {"runs": [_run("T1-ONT", "reads", "ONT", "T1", {
        "A": _entry(outcome="peptide", rna_calls=1, residue=True),
        "B": _entry(outcome="no_call"),
    })]}

    row = variant_funnel(payload, PORTAL)[0]

    assert _stage(row, "covered")["n"] == 2
    assert _stage(row, "allele")["n"] == 1     # A only — B has no alt reads
    assert _stage(row, "called")["n"] == 1
    assert _stage(row, "translated")["n"] == 1
    assert _stage(row, "residue")["n"] == 1


def test_allele_presence_is_a_property_of_the_sample_not_the_method():
    """PacBio has no portal genotyping, so allele support must pool methods.

    Otherwise a sample shows a different number of alleles depending on which
    method looked, which is incoherent: the RNA either carries it or not.
    """
    payload = {"runs": [
        # This method found the allele.
        _run("T1-PB", "reads", "PacBio", "T1", {"A": _entry(rna_calls=1),
                                                 "B": _entry()}),
        # This one did not, and must still report the sample's allele count.
        _run("T1-PB", "assembly", "PacBio", "T1", {"A": _entry(), "B": _entry()}),
    ]}

    rows = variant_funnel(payload, PORTAL)

    assert {_stage(row, "allele")["n"] for row in rows} == {1}


def test_absent_fields_are_pending_not_zero():
    """A stage never computed must not read as a stage that found nothing."""
    payload = {"runs": [_run("T1-ONT", "reads", "ONT", "T1", {"A": _entry()})]}

    row = variant_funnel(payload, PORTAL)[0]

    assert _stage(row, "residue_pick")["pending"] is True
    assert _stage(row, "residue_pick")["n"] is None
    assert _stage(row, "covered")["pending"] is False


def test_single_pick_counted_separately_from_the_ceiling():
    """"Any candidate carries it" and "the chosen one carries it" differ."""
    payload = {"runs": [_run("T1-ONT", "reads", "ONT", "T1", {
        "A": _entry(outcome="peptide", rna_calls=1, residue=True,
                    consensus_residue_confirmed=False),
    })]}

    row = variant_funnel(payload, PORTAL)[0]

    assert _stage(row, "residue")["n"] == 1        # ceiling
    assert _stage(row, "residue_pick")["n"] == 0   # what a caller would get


def test_failed_runs_are_excluded():
    payload = {"runs": [_run("T1-ONT", "reads", "ONT", "T1",
                             {"A": _entry()}, status="failed")]}

    assert variant_funnel(payload, PORTAL) == []
