"""The routes a read can take to become a transcript Exacto will call on.

An "arm" used to mean one of two hard-coded pipelines. It now means one entry
here: a family, the tool that implements it, and the parameters it runs with.
That distinction matters because the most consequential number this test has
produced is not about Exacto at all — 0.6% of the reads entering the assembly
route survive to variant calling, and roughly half of that loss is a threshold
we inherited rather than chose. ``min-read-support 3`` is Nexus's default, and
Nexus is right to have one; whether it is right *here*, where the mutations of
interest are subclonal by nature, is an empirical question nobody has asked.

So methods are data. Adding a threshold sweep or a new assembler is an entry in
``METHODS``, a CI matrix leg, and a row in the comparison — not a code path.

The families:

``reads``
    No assembly. Every variant-spanning read is its own transcript. Highest
    recall measured so far, and the least specific: a median of 5 candidate
    proteins per recovered mutation, a quarter of them frameshifted.

``corrected``
    isONclust groups reads by shared transcript structure, isONcorrect polishes
    each read against its own group. One corrected read per input read, so a
    minority allele is never outvoted — unlike assembly, which averages it away.

``assembly``
    Reads are collapsed into contigs first. Cleanest proteoforms measured (no
    frameshifts at all) and by far the worst recall (3 of 37).
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Nexus's own defaults for filter_rnabloom2_transcripts, verified against the
# script rather than transcribed from documentation.
NEXUS_FILTER_DEFAULTS = {
    "min-mapping-quality": "30",
    "min-read-support": "3",
    "min-fraction-match": "0.5",
    "base-quality": "30",
}


@dataclass(frozen=True)
class Method:
    """One way of turning reads into transcripts, with its parameters."""

    name: str
    family: str
    label: str
    read_types: tuple[str, ...]
    tool: str | None = None
    # Passed through to whichever tool the family uses. Kept as strings because
    # that is how they reach a command line.
    params: dict[str, str] = field(default_factory=dict)
    note: str = ""

    def applies_to(self, read_type: str) -> bool:
        return read_type in self.read_types


LONG = ("long",)
SHORT = ("short",)

METHODS = (
    Method(
        name="reads",
        family="reads",
        label="Reads, uncorrected",
        read_types=LONG,
        note="Each variant-spanning read handed over as its own transcript. "
             "The control that turned out to have the best recall.",
    ),
    Method(
        name="corrected",
        family="corrected",
        label="isONcorrect",
        read_types=LONG,
        tool="isONclust 0.0.6.1 + isONcorrect 0.1.3.5",
        note="Cluster by shared transcript structure, polish each read against "
             "its own cluster. Keeps one read per molecule, so a subclonal "
             "allele survives; reference-free, so a novel junction defines its "
             "own cluster instead of being measured against an annotation.",
    ),
    Method(
        name="isonform",
        family="corrected",
        label="isONform",
        read_types=LONG,
        tool="isONclust + isONcorrect + isONform 0.3.9",
        params={"assemble": "isonform"},
        note="The same clustering, then assembled per cluster rather than left "
             "as reads. Assembly that happens inside an allele-separated "
             "cluster should not average the allele away, which is the "
             "mechanism that costs the global assembler its recall.",
    ),
    Method(
        name="assembly",
        family="assembly",
        label="RNA-Bloom2, Nexus defaults",
        read_types=LONG,
        tool="RNA-Bloom2 2.0.1",
        params=dict(NEXUS_FILTER_DEFAULTS),
        note="The canonical route: Andy Lee's PEPTIDE_PREDICTION_EXACTO "
             "subworkflow, at the thresholds Nexus itself defaults to.",
    ),
    Method(
        name="assembly-ms1",
        family="assembly",
        label="RNA-Bloom2, min read support 1",
        read_types=LONG,
        tool="RNA-Bloom2 2.0.1",
        params={**NEXUS_FILTER_DEFAULTS, "min-read-support": "1"},
        note="Identical but for one threshold. Nexus drops any contig with "
             "fewer than three supporting reads, which costs about half of "
             "them — and falls hardest on exactly the low-VAF mutant contigs "
             "this test is looking for. This asks what that filter is worth.",
    ),
    Method(
        name="spades",
        family="assembly",
        label="rnaSPAdes",
        read_types=SHORT,
        tool="rnaSPAdes 4.0.0",
        note="The only route by which short reads become transcripts. Exacto "
             "never asks where a sequence came from, so this tests whether its "
             "long-read requirement is about the reads or about the assembly "
             "they enable.",
    ),
)

METHODS_BY_NAME = {method.name: method for method in METHODS}
METHOD_NAMES = tuple(method.name for method in METHODS)


def methods_for(read_type: str, requested: list[str] | None = None) -> list[Method]:
    """Methods that can actually run on this kind of read.

    Filtering here rather than in the CI matrix means an impossible pairing is
    skipped loudly instead of producing an empty result, which would read as a
    negative finding rather than as a question never asked.
    """
    names = requested or METHOD_NAMES
    unknown = [name for name in names if name not in METHODS_BY_NAME]
    if unknown:
        raise SystemExit(
            f"unknown method(s): {', '.join(unknown)} — "
            f"known: {', '.join(METHOD_NAMES)}"
        )
    return [
        METHODS_BY_NAME[name]
        for name in names
        if METHODS_BY_NAME[name].applies_to(read_type)
    ]
