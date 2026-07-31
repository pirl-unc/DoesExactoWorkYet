"""Everything the run touched, described for the report.

The site's data-sources page is generated from here rather than written by hand,
so a moved bucket path or a bumped GENCODE release cannot leave the published
provenance quietly wrong.
"""

from __future__ import annotations

from . import config, extract_reads
from .config import NEXUS_VERSION, RNABLOOM_FILTER


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
    sample_rows = []
    for sample in config.SAMPLES:
        stats = extraction.get(sample.name, {})
        genotyping = (
            "The same BAM the portal's own variant table is genotyped from, so "
            "its read counts are directly comparable."
            if sample.portal_genotyped
            else "The portal publishes no genotyping against this BAM — every "
            "number this test reports for it is derived here by Exacto."
        )
        detail = (
            f"{sample.library.capitalize()} RNA-seq of the {sample.timepoint} "
            f"tumour biopsy ({sample.biopsy_date}), sample {sample.biosample}, on "
            f"{sample.platform}. {sample.provenance} {genotyping}"
        )
        taken = None
        if stats:
            taken = (
                f"{stats.get('n_spanning_reads', 0):,} variant-spanning + "
                f"{stats.get('n_context_reads', 0):,} context reads read out over "
                "HTTP byte ranges; the file is never downloaded."
            )
        sample_rows.append(
            {
                "label": f"{sample.name} BAM",
                "url": sample.bam_url,
                "size": _bytes(stats.get("bam_bytes")),
                "detail": detail,
                "taken": taken,
            }
        )
        sample_rows.append(
            {
                "label": f"{sample.name} BAM index",
                "url": sample.bai_url,
                "detail": "Downloaded in full so htslib knows which blocks to request.",
            }
        )

    return [
        {
            "group": "Long-read RNA-seq — the data under test",
            "origin": "osteosarc.com (Backblaze B2)",
            "entries": sample_rows,
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
                {
                    "label": "PacBio Iso-Seq intermediates",
                    "url": f"{config.B2}/ucsf/T1/pacbio_bams/IPISRC044_T1_sclrs/",
                    "detail": "The lima/refine/cluster stages behind the T1-PacBio "
                    "sample. The deduplicated pbmm2-mapped BAM is used instead, for "
                    "the same reason as ONT: it is the finished alignment.",
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


def reproduction() -> list[dict]:
    """The exact commands, in order, to reproduce a run from nothing."""
    return [
        {
            "stage": "Set up",
            "commands": [
                "git clone https://github.com/pirl-unc/DoesExactoWorkYet.git",
                "cd DoesExactoWorkYet",
                "micromamba env create -f environment.yml",
                "micromamba activate does-exacto-work-yet",
                "# EXACTO_VERSION=<git ref> to test something unreleased",
                "bash scripts/install_exacto.sh",
            ],
        },
        {
            "stage": "Run",
            "commands": [
                "export DEWY_WORK_DIR=$PWD/work",
                "python -m pipeline.fetch_osteosarc",
                "python -m pipeline.build_reference",
                "python -m pipeline.extract_reads --samples T1-ONT T2-ONT T3-ONT T1-PacBio",
                (
                    'python -m pipeline.run_exacto --threads "$(nproc)"'
                    " --samples T1-ONT T2-ONT T3-ONT T1-PacBio"
                    " --arms assembly reads"
                ),
                "python -m pipeline.evaluate",
                "python -m pipeline.build_site",
                "python -m http.server -d site 8000",
            ],
        },
        {
            "stage": "Or run it on GitHub",
            "commands": [
                "gh workflow run 'Exacto test' -R pirl-unc/DoesExactoWorkYet",
                (
                    "gh workflow run 'Exacto test' -R pirl-unc/DoesExactoWorkYet"
                    " -f exacto_version=v0.4.6a1"
                ),
                "# otherwise it runs itself every Monday at 04:00 UTC",
            ],
        },
    ]


def configuration() -> list[dict]:
    """Every option this pipeline sets, and how it compares to canonical.

    "Canonical" means Andy Lee's PEPTIDE_PREDICTION_EXACTO subworkflow in Nexus,
    or Exacto's own defaults where Nexus passes nothing. Each row is marked so a
    reader can see at a glance what is stock, what was changed on purpose, and
    what was checked and deliberately left alone.
    """
    return [
        {
            "step": "What a run is",
            "settings": [
                {
                    "name": "Run unit",
                    "value": "one sequencing sample: "
                    + ", ".join(sample.name for sample in config.SAMPLES),
                    "status": "addition",
                    "canonical": "Nexus runs one sample per invocation too",
                    "why": (
                        "Not the biopsy. T1 was sequenced on ONT and again on "
                        "PacBio, so a run keyed on \"T1\" would have to pick one "
                        "platform and drop the other. Each sample is its own CI "
                        "job, its own Exacto run and its own column on the site; "
                        "the timepoint is kept alongside for grouping."
                    ),
                },
                {
                    "name": "Platforms run",
                    "value": ", ".join(
                        sorted({sample.platform for sample in config.SAMPLES})
                    ),
                    "status": "addition",
                    "canonical": "Exacto documents ONT and PacBio long reads",
                    "why": (
                        "Every long-read RNA dataset the portal holds, both of them "
                        "single-cell. The portal never genotyped these mutations "
                        "against the PacBio BAM, so its numbers here are Exacto's "
                        "own rather than a comparison against a published VAF."
                    ),
                },
                {
                    "name": "Reads-arm depth cap",
                    "value": (
                        f"{extract_reads.READS_ARM_READS_PER_VARIANT} reads per "
                        "variant, against "
                        f"{config.SPANNING_READS_PER_VARIANT} for the assembly arm"
                    ),
                    "status": "addition",
                    "canonical": "no such arm — Nexus always assembles first",
                    "why": (
                        "The assembly arm collapses every read into a few hundred "
                        "contigs before Exacto sees them; the reads arm has no such "
                        "step, so call-rna-vars meets the raw depth. Measured on "
                        "T2-ONT it accumulates ~0.55 MB of resident memory per read "
                        "and then allocates a further ~4 GB at the end: at 18,818 "
                        "reads that exceeded the 16 GB runner and killed it mid-step "
                        "after 1h45m. At 600 the arm keeps 42% of the reads, peaks "
                        "near 8 GB, and only the seven deepest variants lose any "
                        "depth — 600 reads still resolves a 5% VAF variant with 30 "
                        "supporting reads. Without this the arm cannot run at all."
                    ),
                },
                {
                    "name": "Base qualities, PacBio",
                    "value": (
                        "QUAL is \"*\" on every Iso-Seq record; a flat Q"
                        f"{extract_reads.SYNTHETIC_BASE_QUALITY} is written in its place"
                    ),
                    "status": "addition",
                    "canonical": (
                        "Nexus writes a flat quality for assembled contigs, which "
                        "have none either — same problem, same answer"
                    ),
                    "why": (
                        "isoseq groupdedup emits deduplicated consensus transcripts "
                        "without per-base quality, and Exacto's call-rna-vars panics "
                        "on a BAM that has none (alignment.rs:229). Dropping those "
                        "records would lose the entire platform — measured across "
                        "three vaccine loci, 206 of 207 PacBio records carry no "
                        "QUAL. These scores are manufactured: nothing downstream "
                        "should read them as measured base quality."
                    ),
                },
                {
                    "name": "Illumina short read",
                    "value": "not run",
                    "status": "checked",
                    "canonical": "Exacto is a long-read tool",
                    "why": (
                        "The bulk and 10x RNA on the portal is all Illumina. It "
                        "supplies the VAF columns in the variant table for context, "
                        "but Exacto's transcript-model construction assumes reads "
                        "that span a transcript, so running it on 100 bp reads "
                        "would test the wrong thing."
                    ),
                },
            ],
        },
        {
            "step": "Reference construction",
            "settings": [
                {
                    "name": "Reference genome",
                    "value": "hg38, masked to the vaccine genes, uncompressed",
                    "status": "addition",
                    "canonical": "a whole reference genome, bgzipped",
                    "why": (
                        "Only the vaccine loci are under test, so only they are "
                        "carved out — at their true hg38 coordinates, so nothing "
                        "downstream needs offset arithmetic. Uncompressed because "
                        "Exacto queries the reference one base at a time and a "
                        "bgzipped file pays a block decompression per query: "
                        "110 s versus 74 s on the same 1,179 reads."
                    ),
                },
                {
                    "name": "Window expansion",
                    "value": f"gene body ±{config.GENE_FLANK_BP:,} bp, then grown "
                    "until every overlapping transcript fits",
                    "status": "addition",
                    "canonical": "not applicable — no masking upstream",
                    "why": (
                        "Exacto drops any read whose candidate transcript touches a "
                        "non-ACGT base, silently. Growing the windows until no "
                        "transcript escapes removes that failure mode; the build "
                        "aborts if one still does."
                    ),
                },
                {
                    "name": "Gene annotation",
                    "value": f"GENCODE v{config.GENCODE_RELEASE}",
                    "status": "deviation",
                    "canonical": "v45 in Nexus's params.yaml",
                    "why": (
                        "v44 is what the source BAMs were aligned against "
                        "(10x refdata-gex-GRCh38-2024-A), so gene models stay "
                        "consistent with the data. Nexus's v45 is a placeholder "
                        "default, not a considered choice for this dataset."
                    ),
                },
            ],
        },
        {
            "step": "Read extraction",
            "settings": [
                {
                    "name": "Variant-spanning reads",
                    "value": f"≤ {config.SPANNING_READS_PER_VARIANT:,} per variant, "
                    "seeded reservoir",
                    "status": "addition",
                    "canonical": "whole-sample FASTQ",
                    "why": (
                        "Coverage spans five orders of magnitude across these genes "
                        "and Exacto's cost scales with read count. The cap sits far "
                        "above what any caller needs to decide a locus; the uncapped "
                        "count is recorded next to it."
                    ),
                },
                {
                    "name": "Context reads",
                    "value": f"≤ {config.CONTEXT_READS_PER_REGION:,} per region",
                    "status": "addition",
                    "canonical": "whole-sample FASTQ",
                    "why": "Filler for the assembler; interchangeable, so capped.",
                },
            ],
        },
        {
            "step": "RNA-Bloom2",
            "settings": [
                {
                    "name": "Extra args",
                    "value": "-chimera",
                    "status": "deviation",
                    "canonical": "-chimera -lrpb",
                    "why": (
                        "-lrpb means the long reads are PacBio. This data is ONT, "
                        "so the flag is dropped. Nexus's reference pipeline is "
                        "built around PacBio HiFi throughout."
                    ),
                },
            ],
        },
        {
            "step": f"nexus_filter_rnabloom2_transcripts ({NEXUS_VERSION})",
            "settings": [
                {
                    "name": "Thresholds",
                    "value": (
                        f"MAPQ ≥ {RNABLOOM_FILTER['min-mapping-quality']}, "
                        f"≥ {RNABLOOM_FILTER['min-read-support']} reads, "
                        f"≥ {RNABLOOM_FILTER['min-fraction-match']} fraction match"
                    ),
                    "status": "canonical",
                    "canonical": "Nexus passes no extra args — these are the defaults",
                    "why": None,
                },
                {
                    "name": "Output",
                    "value": f"FASTA and FASTQ, flat Phred "
                    f"{RNABLOOM_FILTER['base-quality']}",
                    "status": "canonical",
                    "canonical": (
                        "Nexus aligns the FASTQ so the BAM carries QUAL, because "
                        "call-rna-vars panics without it"
                    ),
                    "why": None,
                },
            ],
        },
        {
            "step": "minimap2",
            "settings": [
                {
                    "name": "Assembly arm",
                    "value": f"-ax {config.minimap2_preset('ONT', 'assembly')} "
                    + " ".join(config.MINIMAP2_COMMON_FLAGS),
                    "status": "canonical",
                    "canonical": "identical to Nexus minimap2_rna_args",
                    "why": None,
                },
                {
                    "name": "Reads arm, ONT",
                    "value": f"-ax {config.minimap2_preset('ONT', 'reads')} "
                    + " ".join(config.MINIMAP2_COMMON_FLAGS),
                    "status": "addition",
                    "canonical": "no such arm — Nexus always assembles first",
                    "why": (
                        "splice rather than splice:hq because raw ONT reads are not "
                        "as accurate as polished contigs. This arm exists to "
                        "separate an Exacto miss from an assembler miss."
                    ),
                },
                {
                    "name": "Reads arm, PacBio",
                    "value": f"-ax {config.minimap2_preset('PacBio', 'reads')} "
                    + " ".join(config.MINIMAP2_COMMON_FLAGS),
                    "status": "addition",
                    "canonical": "no such arm — Nexus always assembles first",
                    "why": (
                        "splice:hq here, unlike ONT: the Iso-Seq records are "
                        "deduplicated HiFi consensus transcripts, not raw subreads, "
                        "so the tolerant preset would cost accuracy for nothing."
                    ),
                },
                {
                    "name": "samtools calmd",
                    "value": "not run",
                    "status": "checked",
                    "canonical": "Nexus pipes through samtools calmd to add MD",
                    "why": (
                        "Verified unnecessary: Exacto reads the cs tag "
                        "(alignment.rs:417) and never reads MD anywhere in its "
                        "source."
                    ),
                },
                {
                    "name": "Drop unmapped records",
                    "value": "samtools view -F 4",
                    "status": "addition",
                    "canonical": "Nexus does not filter them",
                    "why": (
                        "remove-unspliced-rnas panics on an unmapped record. The "
                        "canonical pipeline is reachable by the same crash."
                    ),
                },
            ],
        },
        {
            "step": "Exacto",
            "settings": [
                {
                    "name": "Gene and transcript types",
                    "value": ", ".join(config.GENE_TYPES),
                    "status": "canonical",
                    "canonical": "Exacto's default",
                    "why": None,
                },
                {
                    "name": "Gene and transcript levels",
                    "value": ", ".join(config.GENE_LEVELS),
                    "status": "deviation",
                    "canonical": "1, 2 — Exacto's default; Nexus passes nothing",
                    "why": (
                        "GENCODE annotates the mitochondrial genes at level 3, so "
                        "the default would silently drop MT-ND5 — one of the 37 "
                        "vaccine targets."
                    ),
                },
                {
                    "name": "integrate-vars tolerances",
                    "value": "exon 2, transcript boundary 10 kb, intergenic 100 kb",
                    "status": "canonical",
                    "canonical": "Exacto's defaults; Nexus passes nothing",
                    "why": (
                        "Permissive enough that only 19 of 3,359 integrations were "
                        "exact in one measured arm. Left at stock anyway — the point "
                        "is to test Exacto as shipped — and the verdict is keyed on "
                        "exact RNA calls rather than on this table."
                    ),
                },
                {
                    "name": "translate-structs strategy",
                    "value": "longest_orf",
                    "status": "canonical",
                    "canonical": "Nexus params.yaml strategy: longest_orf",
                    "why": None,
                },
                {
                    "name": "call-peptide-vars k",
                    "value": "8 to 11",
                    "status": "canonical",
                    "canonical": "equals Exacto's own defaults",
                    "why": None,
                },
                {
                    "name": "--preset ont",
                    "value": "not set",
                    "status": "checked",
                    "canonical": "Nexus sets --preset pb on call-somatic-dna-vars",
                    "why": (
                        "The preset exists only on the DNA callers, which this "
                        "pipeline does not run. There is no equivalent on "
                        "call-rna-vars."
                    ),
                },
                {
                    "name": "Reference proteome",
                    "value": "GENCODE translations for the tested genes only",
                    "status": "deviation",
                    "canonical": "the whole proteome",
                    "why": (
                        "Keeps call-peptide-vars' k-mer set tractable. It means "
                        "'novel' is relative to the gene's own isoforms rather than "
                        "the whole proteome."
                    ),
                },
                {
                    "name": "DNA variant callset",
                    "value": "the portal's curated somatic calls, supplied",
                    "status": "deviation",
                    "canonical": (
                        "call-somatic-dna-vars on matched tumour/normal long-read WGS"
                    ),
                    "why": (
                        "Sid's WGS is short-read and Exacto's DNA callers want long "
                        "reads. The question asked here is whether Exacto finds these "
                        "mutations in the RNA and translates them, not whether it "
                        "rediscovers them in DNA."
                    ),
                },
            ],
        },
    ]


