"""Everything the run touched, described for the report.

The site's data-sources page is generated from here rather than written by hand,
so a moved bucket path or a bumped GENCODE release cannot leave the published
provenance quietly wrong.
"""

from __future__ import annotations

from . import config
from .run_exacto import ASSEMBLY_BASE_QUALITY


def _bytes(value: int | None) -> str | None:
    if not value:
        return None
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.1f} {unit}".replace(".0 ", " ")
        value /= 1024
    return None


def data_sources(extraction: dict) -> list[dict]:
    """Every input, grouped, with the exact URL it came from."""
    timepoint_rows = []
    for timepoint in config.TIMEPOINTS:
        stats = extraction.get(timepoint.name, {})
        detail = (
            f"Single-cell long-read RNA-seq of the {timepoint.name} tumour biopsy "
            f"({timepoint.biopsy_date}), sample {timepoint.sample}. UMI-deduplicated "
            "and genome-aligned by ONT's wf-single-cell — the same BAM the portal's "
            "own variant table is genotyped from."
        )
        taken = None
        if stats:
            taken = (
                f"{stats.get('n_spanning_reads', 0):,} variant-spanning + "
                f"{stats.get('n_context_reads', 0):,} context reads read out over "
                "HTTP byte ranges; the file is never downloaded."
            )
        timepoint_rows.append(
            {
                "label": f"{timepoint.name} ONT BAM",
                "url": timepoint.bam_url,
                "size": _bytes(stats.get("bam_bytes")),
                "detail": detail,
                "taken": taken,
            }
        )
        timepoint_rows.append(
            {
                "label": f"{timepoint.name} BAM index",
                "url": timepoint.bai_url,
                "detail": "Downloaded in full so htslib knows which blocks to request.",
            }
        )

    return [
        {
            "group": "Long-read RNA-seq — the data under test",
            "origin": "osteosarc.com (Backblaze B2)",
            "entries": timepoint_rows,
        },
        {
            "group": "Vaccine and variant tables",
            "origin": "osteosarc.com",
            "entries": [
                {
                    "label": "Vaccine neoantigen overlap",
                    "url": config.VACCINE_OVERLAP_URL,
                    "detail": "Which mutations went into which of the five personalised "
                    "vaccines, with ELISPOT status and a VAF trend. Defines the 37 "
                    "variants under test.",
                },
                {
                    "label": "Somatic variants, per assay",
                    "url": config.VARIANT_VAFS_URL,
                    "detail": "Ref/alt alleles, consequence and per-assay read counts "
                    "for every curated somatic variant. Supplies the alleles and the "
                    "portal's own genotyping of these same ONT BAMs.",
                },
                {
                    "label": "…its column documentation",
                    "url": config.VARIANT_VAFS_COLUMNS_URL,
                    "detail": "Read to confirm what each column means.",
                },
                {
                    "label": "pVACtools epitopes, MHC class I",
                    "url": config.PVACTOOLS_EPITOPE_URLS["I"],
                    "detail": "The curated neoantigen prediction the vaccine designs "
                    "were picked from. Its MT Epitope Seq column is the peptide that "
                    "was manufactured — searched for verbatim inside Exacto's "
                    "proteoforms. Covers 10 of the 37 variants.",
                },
                {
                    "label": "pVACtools epitopes, MHC class II",
                    "url": config.PVACTOOLS_EPITOPE_URLS["II"],
                    "detail": "Same, for class II peptides.",
                },
            ],
        },
        {
            "group": "Reference sequence and annotation",
            "origin": "osteosarc.com and GENCODE",
            "entries": [
                {
                    "label": "hg38 (GATK Homo_sapiens_assembly38)",
                    "url": config.HG38_FASTA_URL,
                    "detail": "Only the vaccine-gene windows are pulled, over HTTP byte "
                    "ranges — about 8 Mb of the 3.1 Gb genome. Written back at true "
                    "hg38 offsets inside otherwise all-N chromosomes.",
                },
                {
                    "label": "…its .fai index",
                    "url": config.HG38_FASTA_URL + ".fai",
                    "detail": "Tells samtools where each region starts.",
                },
                {
                    "label": f"GENCODE v{config.GENCODE_RELEASE} annotation",
                    "url": config.GENCODE_GTF_URL,
                    "detail": "Matches the 10x refdata-gex-GRCh38-2024-A the source "
                    "BAMs were aligned against, so gene models stay consistent. "
                    "Subset to the analysed windows.",
                },
                {
                    "label": f"GENCODE v{config.GENCODE_RELEASE} translations",
                    "url": config.GENCODE_PROTEINS_URL,
                    "detail": "Wild-type background for call-peptide-vars, subset to "
                    "the tested genes.",
                },
            ],
        },
        {
            "group": "Available on the portal, deliberately not used",
            "origin": "osteosarc.com",
            "entries": [
                {
                    "label": "PacBio single-cell long-read RNA, T1 only",
                    "url": f"{config.B2}/pacbio/"
                    "IPISRC044_T1_sclrs_live_pbmm2_mapped.bam",
                    "detail": "The portal's only PacBio RNA data — T1 alone, also "
                    "single-cell, pbmm2-mapped, with Iso-Seq intermediates under "
                    "ucsf/T1/pacbio_bams/. A useful second platform at one timepoint, "
                    "but it cannot answer the across-timepoints question, so this test "
                    "sticks to ONT.",
                },
                {
                    "label": "Bulk RNA-seq (BostonGene, Tempus, UCLA)",
                    "url": f"{config.B2}/rna-seq/fastq/",
                    "detail": "All Illumina short-read. Exacto is a long-read tool, so "
                    "there is no bulk long-read option here — every long-read RNA "
                    "dataset on the portal is single-cell.",
                },
                {
                    "label": "Raw ONT FASTQs and POD5",
                    "url": f"{config.B2}/ucsf/T1/IPISRC044_T1_sclrs_ONT/",
                    "detail": "Pre-basecalling and pre-alignment. The dedup BAMs are "
                    "used instead so read support is directly comparable to the "
                    "portal's own genotyping.",
                },
            ],
        },
        {
            "group": "The tool under test",
            "origin": "GitHub",
            "entries": [
                {
                    "label": "Exacto",
                    "url": "https://github.com/pirl-unc/exacto",
                    "detail": "Installed from the newest published release unless "
                    "EXACTO_VERSION names a git ref. The resolved version is recorded "
                    "in results/environment.json and shown throughout.",
                },
                {
                    "label": "This harness",
                    "url": "https://github.com/pirl-unc/DoesExactoWorkYet",
                    "detail": "Every step below is a module in pipeline/.",
                },
            ],
        },
    ]


def parameters() -> list[dict]:
    """The knobs, with the reasoning, so a number on the page can be traced."""
    return [
        {
            "step": "build_reference",
            "settings": [
                ("Gene window padding", f"±{config.GENE_FLANK_BP:,} bp"),
                (
                    "Window expansion",
                    "grown until every overlapping GENCODE transcript sits entirely "
                    "on real sequence, because Exacto drops reads whose candidate "
                    "transcript touches an N",
                ),
                ("Reference format", "uncompressed FASTA — bgzip costs ~50% runtime"),
            ],
        },
        {
            "step": "extract_reads",
            "settings": [
                (
                    "Variant-spanning reads",
                    f"kept up to {config.SPANNING_READS_PER_VARIANT:,} per variant, "
                    "seeded reservoir sample",
                ),
                (
                    "Context reads",
                    f"kept up to {config.CONTEXT_READS_PER_REGION:,} per region",
                ),
                ("Excluded", "secondary, supplementary and unmapped records"),
                ("Retry", "4 attempts per region, region restarted from scratch"),
            ],
        },
        {
            "step": "alignment (minimap2)",
            "settings": [
                ("Assembly arm preset", f"-ax {config.MINIMAP2_PRESET['assembly']}"),
                ("Reads arm preset", f"-ax {config.MINIMAP2_PRESET['reads']}"),
                ("Flags", " ".join(config.MINIMAP2_COMMON_FLAGS)),
                (
                    "Post-processing",
                    "samtools view -F 4 to drop unmapped records, which crash Exacto",
                ),
                (
                    "Assembled contig quality",
                    f"flat Phred '{ASSEMBLY_BASE_QUALITY}' (Q40) — Exacto panics on a "
                    "BAM with no QUAL, and RNA-Bloom2 emits FASTA",
                ),
            ],
        },
        {
            "step": "Exacto",
            "settings": [
                ("Gene types", ", ".join(config.GENE_TYPES)),
                (
                    "Gene / transcript levels",
                    ", ".join(config.GENE_LEVELS)
                    + " — level 3 included so the mitochondrial genes, and MT-ND5, "
                    "are not silently dropped",
                ),
                (
                    "integrate-vars tolerances",
                    "exon offset 2; transcript-boundary and intergenic offsets set to "
                    "0, against defaults of 10 kb and 100 kb that pair nearly "
                    "everything with everything",
                ),
                ("translate-structs strategy", "longest_orf"),
                ("call-peptide-vars k", "8–11"),
            ],
        },
    ]
