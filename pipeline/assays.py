"""Collapse the portal's per-sample read counts into a VAF matrix.

``variant_vafs_long.tsv`` carries one row per variant per sequencing sample —
200 rows for each of the 37 vaccine variants. The site wants a compact grid:
one cell per assay type and timepoint, so DNA and RNA support can be read across
the disease course at a glance.

Several samples can share an assay and timepoint (three WGS providers at T0, two
WES at T1, three bulk RNA runs at T1). Those are pooled by summing reads, which
is the honest summary when the same locus is sequenced repeatedly; the
per-sample rows stay available in ``results/vaccine_variants.json``.
"""

from __future__ import annotations

# The assay under test leads, then the rest of the RNA, then the DNA that
# establishes each variant is really there.
ASSAY_ORDER = ["scRNA_ONT", "scRNA_PacBio", "RNA", "scRNA", "WGS", "WES"]

# Every column states its platform and whether it is single-cell, rather than
# leaving either to be inferred from a terse name. "Bulk RNA" in particular
# reads as though it might be long-read and is not: every bulk sample is
# STAR- or oncoanalyser-aligned Illumina from BostonGene, Tempus, Personalis or
# UCLA. Of the five groups exactly one is long-read.
# The platform is part of the name, never left implicit and never collapsed into
# a bare "Bulk RNA" that could be read as any technology. Two groups labelled
# only "scRNA" would be worse still — one is Illumina, one is ONT, and the whole
# point of the test is the difference between them.
ASSAY_META = {
    "scRNA_ONT": {"label": "scRNA · ONT", "platform": "ONT",
                  "read": "long read", "single_cell": True, "kind": "rna"},
    # No portal rows today — kept so that if osteosarc.com ever genotypes these
    # mutations against the Iso-Seq BAM, the column names itself correctly
    # instead of falling through to the "?" label.
    "scRNA_PacBio": {"label": "scRNA · PacBio", "platform": "PacBio",
                     "read": "long read", "single_cell": True, "kind": "rna"},
    "RNA":       {"label": "Bulk RNA · ILMN", "platform": "Illumina",
                  "read": "short read", "single_cell": False, "kind": "rna"},
    "scRNA":     {"label": "scRNA · ILMN 10x", "platform": "Illumina 10x",
                  "read": "short read", "single_cell": True, "kind": "rna"},
    "WGS":       {"label": "WGS · ILMN", "platform": "Illumina",
                  "read": "short read", "single_cell": False, "kind": "dna"},
    "WES":       {"label": "WES · ILMN", "platform": "Illumina",
                  "read": "short read", "single_cell": False, "kind": "dna"},
}
FALLBACK_META = {
    "label": "?", "platform": "unknown", "read": "unknown",
    "single_cell": False, "kind": "rna",
}

# The long-read assays this test actually runs Exacto on.
TESTED_ASSAYS = {"scRNA_ONT", "scRNA_PacBio"}

TIMEPOINT_ORDER = ["T0", "T1", "T1-organoid", "T2", "T3"]


def _sort_key(assay: str, timepoint: str | None) -> tuple[int, int, str]:
    assay_rank = ASSAY_ORDER.index(assay) if assay in ASSAY_ORDER else len(ASSAY_ORDER)
    label = timepoint or ""
    tp_rank = (
        TIMEPOINT_ORDER.index(label) if label in TIMEPOINT_ORDER else len(TIMEPOINT_ORDER)
    )
    return (assay_rank, tp_rank, label)


def columns(variants: list[dict]) -> list[dict]:
    """The assay/timepoint pairs the portal actually has tumour data for.

    Derived from the data rather than hard-coded, so a new provider or timepoint
    on osteosarc.com shows up as a new column instead of being dropped.
    """
    seen: set[tuple[str, str | None]] = set()
    for variant in variants:
        for entry in variant.get("assay_support", []):
            if entry["tissue"] != "tumor":
                continue
            seen.add((entry["assay_type"], entry["timepoint"]))

    columns = []
    for assay, timepoint in sorted(seen, key=lambda item: _sort_key(*item)):
        meta = ASSAY_META.get(assay, {**FALLBACK_META, "label": assay})
        columns.append(
            {
                "key": f"{assay}|{timepoint or ''}",
                "assay": assay,
                "assay_label": meta["label"],
                "platform": meta["platform"],
                "read_length": meta["read"],
                "single_cell": meta["single_cell"],
                "long_read": meta["read"] == "long read",
                "kind": meta["kind"],
                "timepoint": timepoint or "—",
                "tested": assay in TESTED_ASSAYS,
            }
        )
    return columns


# Platforms the portal holds sequencing for but never genotyped these variants
# against, so they have no VAF to show. Stated explicitly: a column that simply
# is not there reads as an oversight.
ABSENT_PLATFORMS = [
    {
        "platform": "PacBio",
        "reason": (
            "The portal's variant table contains no PacBio rows, so it publishes "
            "no VAF for these 37 mutations on that platform and there is nothing "
            "to put in this grid. That is a gap in the published genotyping, not "
            "in the data: the portal does hold a PacBio Iso-Seq BAM for T1 "
            "(IPISRC044_T1_sclrs_live_pbmm2_mapped.bam), and this test now runs "
            "Exacto on it as the T1-PacBio sample. Its numbers are derived here "
            "rather than read off the portal — see the Exacto columns."
        ),
    },
]


def absent_platforms(variants: list[dict]) -> list[dict]:
    """Only report a platform as absent if it really is absent from the data."""
    present = {
        entry["assay_type"].lower()
        for variant in variants
        for entry in variant.get("assay_support", [])
    }
    return [
        item
        for item in ABSENT_PLATFORMS
        if not any(item["platform"].lower() in name for name in present)
    ]


def matrix(variant: dict) -> dict[str, dict]:
    """Pooled tumour VAF and read counts for one variant, keyed by column."""
    pooled: dict[str, dict] = {}
    for entry in variant.get("assay_support", []):
        if entry["tissue"] != "tumor":
            continue
        key = f"{entry['assay_type']}|{entry['timepoint'] or ''}"
        cell = pooled.setdefault(
            key, {"alt": 0, "total": 0, "samples": 0, "covered": 0}
        )
        cell["samples"] += 1
        if entry["total_reads"]:
            cell["alt"] += entry["alt_reads"] or 0
            cell["total"] += entry["total_reads"]
            cell["covered"] += 1

    for cell in pooled.values():
        cell["vaf"] = round(cell["alt"] / cell["total"], 4) if cell["total"] else None
    return pooled


def germline(variant: dict) -> dict[str, dict]:
    """The same, for the matched blood normals.

    Worth showing in the detail: a somatic call with alt reads in blood is a
    germline variant that slipped through, and the table should not hide that.
    """
    pooled: dict[str, dict] = {}
    for entry in variant.get("assay_support", []):
        if entry["tissue"] == "tumor":
            continue
        key = f"{entry['assay_type']}|{entry['timepoint'] or ''}"
        cell = pooled.setdefault(key, {"alt": 0, "total": 0, "samples": 0})
        cell["samples"] += 1
        if entry["total_reads"]:
            cell["alt"] += entry["alt_reads"] or 0
            cell["total"] += entry["total_reads"]

    for cell in pooled.values():
        cell["vaf"] = round(cell["alt"] / cell["total"], 4) if cell["total"] else None
    return pooled
