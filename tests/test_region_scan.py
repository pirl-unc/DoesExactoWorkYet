"""Checks on the region scan: sampling caps, dedup, and retry behaviour.

The scan is the one place that talks to a 37 GB remote BAM, so it has to be
both reproducible — same reads every run, or the whole test stops being a
regression check — and tolerant of a dropped connection partway through a
region.
"""

from __future__ import annotations

import pytest

from pipeline import extract_reads
from pipeline.config import SAMPLES_BY_NAME


class FakeRead:
    def __init__(self, name, start, end, flags=()):
        self.query_name = name
        self.reference_start = start - 1  # pysam is 0-based here
        self.reference_end = end
        self.is_secondary = "secondary" in flags
        self.is_supplementary = "supplementary" in flags
        self.is_unmapped = "unmapped" in flags
        self._sequence = "ACGT" * 25

    def get_forward_sequence(self):
        return self._sequence

    def get_forward_qualities(self):
        return [40] * len(self._sequence)


class FakeBam:
    """Hands back a fixed read list, optionally failing the first N fetches."""

    def __init__(self, reads, fail_times=0, fail_after=0):
        self.reads = reads
        self.fail_times = fail_times
        self.fail_after = fail_after
        self.fetch_calls = 0

    def fetch(self, chrom, start, end):
        self.fetch_calls += 1
        should_fail = self.fetch_calls <= self.fail_times
        for index, read in enumerate(self.reads):
            if should_fail and index >= self.fail_after:
                raise OSError("Error in the HTTP2 framing layer")
            yield read


SAMPLE = SAMPLES_BY_NAME["T1-ONT"]
REGION = {"chrom": "chr1", "start": 1, "end": 10_000, "genes": ["GENE"]}
VARIANT = {"variant_id": "GENE-chr1-500", "chrom": "chr1", "pos": 500, "ref": "C"}
SPANS = [(VARIANT, 500, 500)]


@pytest.fixture(autouse=True)
def small_caps(monkeypatch):
    monkeypatch.setattr(extract_reads, "SPANNING_READS_PER_VARIANT", 5)
    monkeypatch.setattr(extract_reads, "CONTEXT_READS_PER_REGION", 4)
    monkeypatch.setattr(extract_reads, "REGION_FETCH_BACKOFF_SECONDS", 0)


def spanning_reads(count, offset=0):
    return [FakeRead(f"span{i + offset}", 400, 600) for i in range(count)]


def context_reads(count, offset=0):
    return [FakeRead(f"ctx{i + offset}", 2000, 2200) for i in range(count)]


def scan(bam, already_seen=None):
    return extract_reads.scan_region(
        bam, SAMPLE, REGION, SPANS, already_seen or set()
    )


def test_keeps_everything_under_the_caps():
    result = scan(FakeBam(spanning_reads(3) + context_reads(2)))

    assert len(result["spanning"][VARIANT["variant_id"]]) == 3
    assert result["spanning_seen"][VARIANT["variant_id"]] == 3
    assert len(result["context"]) == 2
    assert result["context_seen"] == 2


def test_caps_are_enforced_but_the_true_count_is_still_reported():
    result = scan(FakeBam(spanning_reads(20) + context_reads(30)))

    assert len(result["spanning"][VARIANT["variant_id"]]) == 5
    assert result["spanning_seen"][VARIANT["variant_id"]] == 20
    assert len(result["context"]) == 4
    assert result["context_seen"] == 30


def test_sampling_is_reproducible():
    reads = spanning_reads(20) + context_reads(30)
    first = scan(FakeBam(list(reads)))
    second = scan(FakeBam(list(reads)))

    assert first["spanning"] == second["spanning"]
    assert first["context"] == second["context"]


def test_skips_secondary_supplementary_and_unmapped():
    reads = [
        FakeRead("primary", 400, 600),
        FakeRead("secondary", 400, 600, flags=("secondary",)),
        FakeRead("supp", 400, 600, flags=("supplementary",)),
        FakeRead("unmapped", 400, 600, flags=("unmapped",)),
    ]

    result = scan(FakeBam(reads))

    assert result["names"] == {"primary"}


def test_skips_reads_already_taken_by_a_neighbouring_region():
    result = scan(FakeBam(spanning_reads(3)), already_seen={"span0"})

    assert result["names"] == {"span1", "span2"}
    assert result["spanning_seen"][VARIANT["variant_id"]] == 2


def test_a_dropped_connection_is_retried_from_a_clean_slate():
    bam = FakeBam(spanning_reads(4) + context_reads(3), fail_times=1, fail_after=3)

    result = scan(bam)

    assert bam.fetch_calls == 2
    # The partial first attempt must not have leaked into the counts.
    assert result["spanning_seen"][VARIANT["variant_id"]] == 4
    assert result["context_seen"] == 3
    assert len(result["names"]) == 7


def test_a_retried_region_samples_identically_to_a_clean_one():
    reads = spanning_reads(20) + context_reads(30)
    clean = scan(FakeBam(list(reads)))
    flaky = scan(FakeBam(list(reads), fail_times=2, fail_after=5))

    assert flaky["spanning"] == clean["spanning"]
    assert flaky["context"] == clean["context"]


def test_gives_up_after_the_attempt_limit():
    bam = FakeBam(spanning_reads(4), fail_times=99, fail_after=1)

    with pytest.raises(OSError):
        scan(bam)

    assert bam.fetch_calls == extract_reads.REGION_FETCH_ATTEMPTS
