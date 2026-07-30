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

# Ordered so the table reads DNA → bulk RNA → single-cell RNA → long read,
# which is roughly increasing relevance to what this test actually measures.
ASSAY_ORDER = ["WGS", "WES", "RNA", "scRNA", "scRNA_ONT"]

ASSAY_LABELS = {
    "WGS": "WGS",
    "WES": "WES",
    "RNA": "Bulk RNA",
    "scRNA": "scRNA 10x",
    "scRNA_ONT": "scRNA ONT",
}

# Spelled out because "Bulk RNA" reads as though it might be long-read, and it is
# not: every bulk sample is STAR- or oncoanalyser-aligned Illumina, from
# BostonGene, Tempus, Personalis or UCLA. Of the five groups only scRNA ONT is
# long-read, which is why it is the one Exacto is pointed at.
ASSAY_PLATFORM = {
    "WGS": "Illumina, short read",
    "WES": "Illumina, short read",
    "RNA": "Illumina, short read",
    "scRNA": "Illumina 10x, short read",
    "scRNA_ONT": "Oxford Nanopore, long read",
}

ASSAY_KIND = {
    "WGS": "dna",
    "WES": "dna",
    "RNA": "rna",
    "scRNA": "rna",
    "scRNA_ONT": "rna",
}

# The long-read assay this test actually runs Exacto on.
TESTED_ASSAY = "scRNA_ONT"

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

    return [
        {
            "key": f"{assay}|{timepoint or ''}",
            "assay": assay,
            "assay_label": ASSAY_LABELS.get(assay, assay),
            "platform": ASSAY_PLATFORM.get(assay, "unknown"),
            "kind": ASSAY_KIND.get(assay, "rna"),
            "timepoint": timepoint or "—",
            "tested": assay == TESTED_ASSAY,
        }
        for assay, timepoint in sorted(seen, key=lambda item: _sort_key(*item))
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
