"""Central configuration: where the data lives and how the test is parameterised.

Everything downstream reads from here so that a change of GENCODE release or a
moved bucket path is a one-line edit.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# Repository layout
# --------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
SITE_DIR = REPO_ROOT / "site"

# Big, regenerable intermediates. Overridable so CI can put them on a scratch
# disk instead of inside the checkout.
WORK_DIR = Path(os.environ.get("DEWY_WORK_DIR", REPO_ROOT / "work")).resolve()

def ensure_ca_bundle() -> str | None:
    """Point libcurl at a CA bundle before htslib opens an https:// file.

    pysam ships its own libcurl, and inside a conda environment that build often
    cannot find the system trust store — it fails with "Libcurl reported error 77
    (Problem with the SSL CA cert)" and htslib turns that into a bare I/O error.
    Conda's samtools is fine, which makes the failure look baffling: the same URL
    works from the shell and not from Python.

    Returns the bundle in use, or None if the environment already had one.
    """
    if os.environ.get("CURL_CA_BUNDLE"):
        return None

    candidates = []
    try:
        import certifi

        candidates.append(certifi.where())
    except ImportError:
        pass
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(str(Path(conda_prefix) / "ssl" / "cacert.pem"))
    candidates += [
        "/etc/ssl/certs/ca-certificates.crt",  # Debian, Ubuntu
        "/etc/pki/tls/certs/ca-bundle.crt",  # RHEL, Fedora
    ]

    for candidate in candidates:
        if candidate and Path(candidate).exists():
            os.environ["CURL_CA_BUNDLE"] = candidate
            os.environ.setdefault("SSL_CERT_FILE", candidate)
            return candidate
    return None


# --------------------------------------------------------------------------
# osteosarc.com — Sid Sijbrandij's open osteosarcoma data portal
# --------------------------------------------------------------------------

OSTEOSARC = "https://osteosarc.com"
B2 = "https://b2.osteosarc.com"

# Neoantigens selected into each of the five personalised vaccines, with
# ELISPOT status and a VAF trend per mutation. Source of the variant table.
VACCINE_OVERLAP_URL = f"{OSTEOSARC}/data/vaccine_overlap.json"

# Long-format per-assay read counts / VAFs for every curated somatic variant.
# Supplies ref/alt alleles, consequence and per-timepoint DNA + RNA support.
VARIANT_VAFS_URL = f"{OSTEOSARC}/variants/variant_vafs_long.tsv"
VARIANT_VAFS_COLUMNS_URL = f"{OSTEOSARC}/variants/variant_vafs_long.columns.tsv"

# The curated pVACtools neoantigen prediction the vaccine designs drew on. Its
# "MT Epitope Seq" column is the closest thing published to the actual peptide
# sequences that were manufactured, which lets the test ask whether Exacto's
# proteoform literally contains the vaccine epitope rather than merely carrying
# the right amino-acid change. It covers a subset of the vaccine variants.
PVACTOOLS_BASE = (
    f"{B2}/neoantigen_prediction/pvactools/"
    "2025.04.25.sg.curated.neoantigen.predictions"
)
PVACTOOLS_EPITOPE_URLS = {
    "I": f"{PVACTOOLS_BASE}/MHC_Class_I/"
    "SG.WGS_SG.WGS.UCLA.2025.01.tumor.all_epitopes.tsv",
    "II": f"{PVACTOOLS_BASE}/MHC_Class_II/"
    "SG.WGS_SG.WGS.UCLA.2025.01.tumor.all_epitopes.tsv",
}

# GATK hg38 as distributed by the portal. Served with byte-range support, so
# `samtools faidx <url> <region>` pulls only the slices we ask for.
HG38_FASTA_URL = f"{B2}/ref_genome/Homo_sapiens_assembly38.fasta"

# --------------------------------------------------------------------------
# The long-read RNA-seq samples Exacto is run on
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Sample:
    """One long-read RNA-seq dataset — the unit a run is keyed on.

    Deliberately *not* the biopsy. T1 was sequenced twice, on ONT and on PacBio,
    so a run keyed on "T1" has to pick one platform and silently drop the other;
    keying on the sample makes a second platform at the same timepoint a config
    entry rather than a special case. ``timepoint`` is kept alongside so the
    biopsies still group back together for the across-time question.

    ``name`` is the identifier everywhere else: FASTQ names, the work/exacto
    directory, the CI matrix leg, the scored JSON, the site's columns.
    """

    name: str
    timepoint: str
    platform: str
    # "long" or "short". Decides which preprocessing an arm can offer: a 150 bp
    # read is not a transcript, so the reads and corrected arms are meaningless
    # for it, while assembly is the only way it can reach Exacto at all.
    read_type: str
    # Key into assays.ASSAY_META, so a sample and a portal VAF column that mean
    # the same sequencing agree on their label rather than drifting apart.
    assay: str
    label: str
    biopsy_date: str
    biosample: str
    bam_url: str
    library: str
    # Whether the portal's own variant table genotyped these mutations against
    # this sample. False means the only numbers for it are the ones this test
    # derives itself.
    portal_genotyped: bool
    provenance: str

    @property
    def bai_url(self) -> str:
        return self.bam_url + ".bai"


def _ont_bam(timepoint: str) -> str:
    """UMI-deduplicated, genome-aligned reads from ONT's wf-single-cell.

    These are the same BAMs the portal's own variant table is genotyped from,
    so read support numbers are directly comparable.
    """
    stem = f"IPISRC044_{timepoint}_sclrs_ONT"
    return (
        f"{B2}/ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/"
        f"{stem}/{stem}/{stem}_dedup/{stem}_dedup.bam"
    )


_ONT_PROVENANCE = (
    "ONT wf-single-cell: barcode/UMI-tagged, deduplicated, aligned to "
    "refdata-gex-GRCh38-2024-A with minimap2 -ax splice."
)

SAMPLES = (
    Sample(
        name="T1-ONT", timepoint="T1", platform="ONT", assay="scRNA_ONT", read_type="long",
        label="T1 · ONT", biopsy_date="2024-06", biosample="IPISRC044_T1",
        bam_url=_ont_bam("T1"), library="single-cell, long read",
        portal_genotyped=True, provenance=_ONT_PROVENANCE,
    ),
    Sample(
        name="T2-ONT", timepoint="T2", platform="ONT", assay="scRNA_ONT", read_type="long",
        label="T2 · ONT", biopsy_date="2025-01", biosample="IPISRC044_T2",
        bam_url=_ont_bam("T2"), library="single-cell, long read",
        portal_genotyped=True, provenance=_ONT_PROVENANCE,
    ),
    Sample(
        name="T3-ONT", timepoint="T3", platform="ONT", assay="scRNA_ONT", read_type="long",
        label="T3 · ONT", biopsy_date="2025-04", biosample="IPISRC044_T3",
        bam_url=_ont_bam("T3"), library="single-cell, long read",
        portal_genotyped=True, provenance=_ONT_PROVENANCE,
    ),
    # The portal's variant table has no PacBio rows, so there is no published
    # VAF for these 37 mutations on this platform. That is a reason to run
    # Exacto on it, not a reason to leave it out: the BAM is right there, it is
    # the same T1 biopsy, and Exacto's own RNA caller can supply the numbers the
    # portal never published. It is also the more favourable input on paper —
    # Iso-Seq consensus transcripts are full-length and HiFi-accurate — which
    # makes it a real control on how much of a miss is ONT's error rate.
    Sample(
        name="T1-PacBio", timepoint="T1", platform="PacBio", assay="scRNA_PacBio",
        read_type="long",
        label="T1 · PacBio Iso-Seq", biopsy_date="2024-06",
        biosample="IPISRC044_T1_sclrs_live",
        bam_url=f"{B2}/pacbio/IPISRC044_T1_sclrs_live_pbmm2_mapped.bam",
        library="single-cell, long read",
        portal_genotyped=False,
        provenance=(
            "PacBio Iso-Seq v4: lima 5p--3p, refine, cluster, groupdedup "
            "(--keep-non-real-cells), aligned with pbmm2 1.14 --preset ISOSEQ "
            "to refdata-gex-GRCh38-2020-A. Records are deduplicated consensus "
            "transcripts, not raw subreads."
        ),
    ),
)

SAMPLES = SAMPLES + (
    # Illumina bulk RNA-seq. Exacto calls itself a long-read toolkit and it is
    # right to: a 150 bp read cannot carry a transcript structure, so raw short
    # reads have no meaningful path in. Assembled contigs do — they are
    # transcript-like sequences and Exacto does not ask where a sequence came
    # from. Including this is the only way to answer whether the long-read
    # requirement is about the reads or about the assembly they enable.
    #
    # 2025.01.06 is the UCLA resection matching the T2 biopsy; it is the one
    # bulk RNA alignment on the portal that is coordinate-sorted and indexed,
    # which byte-range access requires.
    Sample(
        name="T2-ILMN", timepoint="T2", platform="Illumina", assay="RNA",
        read_type="short",
        label="T2 · Illumina bulk", biopsy_date="2025-01",
        biosample="sj.rna.2025.01.resection.ucla",
        bam_url=f"{B2}/genomics/genomics-bulk/2025.01.06/RNA/2025.01.06.rna."
                "ucla-core/processed/STAR/25.03.23.rna.ucla.2025.01.resection."
                "tcga.d32.protocolAligned.sorted.bam",
        library="bulk, short read",
        portal_genotyped=True,
        provenance=(
            "STAR, TCGA protocol, GENCODE d32, UCLA Clinical Genomics. "
            "Coordinate-sorted and indexed, so the same byte-range access works."
        ),
    ),
)

SAMPLES_BY_NAME = {sample.name: sample for sample in SAMPLES}

# The biopsies, for grouping samples back into the across-time question.
TIMEPOINT_ORDER = ("T1", "T2", "T3")


def samples_named(names: list[str] | None) -> list[Sample]:
    """Resolve CLI/CI sample names, failing loudly on a typo.

    A silently-empty selection would make a three-hour CI leg finish in seconds
    and report nothing missing.
    """
    if not names:
        return list(SAMPLES)
    unknown = [name for name in names if name not in SAMPLES_BY_NAME]
    if unknown:
        raise SystemExit(
            f"unknown sample(s): {', '.join(unknown)} — "
            f"known: {', '.join(SAMPLES_BY_NAME)}"
        )
    return [SAMPLES_BY_NAME[name] for name in names]

# --------------------------------------------------------------------------
# Reference annotation
# --------------------------------------------------------------------------

# The source BAMs were aligned against 10x's refdata-gex-GRCh38-2024-A, which
# is built from GENCODE v44 — matching it keeps gene models consistent.
GENCODE_RELEASE = "44"
GENCODE_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    f"release_{GENCODE_RELEASE}"
)
GENCODE_GTF_URL = f"{GENCODE_BASE}/gencode.v{GENCODE_RELEASE}.annotation.gtf.gz"
GENCODE_PROTEINS_URL = (
    f"{GENCODE_BASE}/gencode.v{GENCODE_RELEASE}.pc_translations.fa.gz"
)

# Padding around each vaccine-mutation gene when carving the reference and
# pulling reads. build_reference then grows the windows further until every
# transcript that overlaps one sits entirely on real sequence.
GENE_FLANK_BP = 10_000

# Coverage across the vaccine loci spans five orders of magnitude — the
# mitochondrial window holds ~1.3M reads against VPS13B's 20 — and Exacto's
# transcript-model construction costs roughly a third of a CPU-second per read,
# so the deepest loci would otherwise set the runtime for everything.
#
# Reads that do not touch a variant only exist to give RNA-Bloom2 something to
# extend transcripts with, so they are capped per region.
CONTEXT_READS_PER_REGION = 5_000

# Reads that do cover a variant are capped per variant instead. 3,000 reads is
# far more than any caller needs to decide a locus — it still resolves a variant
# sitting at 1% VAF with 30 supporting reads — and it keeps the hyper-expressed
# loci from dominating. The uncapped count is recorded alongside so the site can
# show what was sampled from.
SPANNING_READS_PER_VARIANT = 3_000

# --------------------------------------------------------------------------
# Exacto
# --------------------------------------------------------------------------

# Which Exacto is under test is chosen by scripts/install_exacto.sh, from the
# EXACTO_VERSION and EXACTO_REPO environment variables. The resolved version
# lands in results/environment.json and is shown on the site.

# Two ways of feeding long reads to Exacto, both reported:
#   assembly — RNA-Bloom2 assembles reads into full-length transcripts first,
#              which is what the Exacto docs prescribe for polyA long reads.
#   reads    — the reads themselves are handed to Exacto as "transcripts",
#              skipping assembly. Cheaper, and a useful control.
ARMS = ("assembly", "reads", "corrected")

# Which arms make sense for which kind of read. Enforced rather than left to a
# CI matrix comment, so an invalid pairing fails loudly instead of producing an
# empty result that looks like a negative finding.
ARMS_BY_READ_TYPE = {
    "long": ("assembly", "reads", "corrected"),
    # No reads arm: a 150 bp read is not a transcript. No corrected arm:
    # isONclust/isONcorrect cluster and polish using reads that span shared
    # transcript structure, which short reads do not.
    "short": ("assembly",),
}


def arms_for(sample: Sample, requested: list[str] | None = None) -> list[str]:
    allowed = ARMS_BY_READ_TYPE.get(sample.read_type, ARMS)
    return [arm for arm in (requested or ARMS) if arm in allowed]

# Andy Lee's Nexus wraps Exacto in a canonical Nextflow subworkflow
# (PEPTIDE_PREDICTION_EXACTO). Step 6 of it filters RNA-Bloom2's output before
# anything is aligned, on these thresholds — without it Exacto is handed every
# contig the assembler emitted, single-read junk included.
NEXUS_VERSION = "v0.2.0a7"
RNABLOOM_FILTER = {
    "min-mapping-quality": "30",
    "min-read-support": "3",
    "min-fraction-match": "0.5",
    # The same trick this harness arrived at independently: assembled contigs
    # carry no per-base quality, and call-rna-vars panics without one.
    "base-quality": "30",
}

# The portal's BAMs were aligned with plain `minimap2 -ax splice --MD`, without
# --cs, and Exacto reads the CS tag to find variants — so everything gets
# realigned here regardless of arm.
#
# Preset differs by platform and arm. splice:hq is minimap2's high-accuracy
# spliced mode; it is right for assembled contigs whatever assembled them, and
# right for PacBio Iso-Seq reads, which are deduplicated HiFi consensus
# transcripts rather than raw subreads. Raw ONT reads are the one case that
# needs the tolerant preset — using splice:hq on them would cost real alignments.
MINIMAP2_PRESET = {
    ("ONT", "assembly"): "splice:hq",
    ("ONT", "reads"): "splice",
    # Corrected reads are consensus-polished, so they earn the high-quality
    # preset the raw ONT reads do not.
    ("ONT", "corrected"): "splice:hq",
    ("PacBio", "assembly"): "splice:hq",
    ("PacBio", "reads"): "splice:hq",
    ("PacBio", "corrected"): "splice:hq",
    # rnaSPAdes contigs are short-read assemblies: accurate, so splice:hq.
    ("Illumina", "assembly"): "splice:hq",
}
MINIMAP2_COMMON_FLAGS = ("-uf", "--cs", "--eqx", "-Y", "-L", "--secondary=no")


def minimap2_preset(platform: str, arm: str) -> str:
    """Fail rather than guess: a new platform must state its own preset."""
    try:
        return MINIMAP2_PRESET[(platform, arm)]
    except KeyError as error:
        raise SystemExit(
            f"no minimap2 preset for platform {platform!r} arm {arm!r}"
        ) from error

# GENCODE annotates the mitochondrial genes at level 3, so the Exacto defaults
# (levels 1 and 2) would silently drop MT-ND5 — one of the vaccine targets.
GENE_LEVELS = ("1", "2", "3")
GENE_TYPES = ("protein_coding",)
