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
# Long-read RNA-seq: ONT single-cell long-read sequencing of the tumour
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Timepoint:
    """One tumour biopsy that has ONT long-read RNA-seq."""

    name: str
    label: str
    biopsy_date: str
    sample: str

    @property
    def bam_url(self) -> str:
        """UMI-deduplicated, genome-aligned reads from ONT's wf-single-cell.

        These are the same BAMs the portal's own variant table is genotyped
        from, so read support numbers are directly comparable.
        """
        stem = f"IPISRC044_{self.name}_sclrs_ONT"
        return (
            f"{B2}/ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/"
            f"{stem}/{stem}/{stem}_dedup/{stem}_dedup.bam"
        )

    @property
    def bai_url(self) -> str:
        return self.bam_url + ".bai"


TIMEPOINTS = (
    Timepoint("T1", "T1 — first recurrence biopsy", "2024-06", "IPISRC044_T1"),
    Timepoint("T2", "T2 — second recurrence biopsy", "2025-01", "IPISRC044_T2"),
    Timepoint("T3", "T3 — third recurrence biopsy", "2025-04", "IPISRC044_T3"),
)

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
ARMS = ("assembly", "reads")

# The portal's BAMs were aligned with plain `minimap2 -ax splice --MD`, without
# --cs, and Exacto reads the CS tag to find variants — so everything gets
# realigned here regardless of arm.
#
# Preset differs by arm: assembled contigs are accurate enough for splice:hq,
# raw ONT reads are not.
MINIMAP2_PRESET = {"assembly": "splice:hq", "reads": "splice"}
MINIMAP2_COMMON_FLAGS = ("-uf", "--cs", "--eqx", "-Y", "-L", "--secondary=no")

# GENCODE annotates the mitochondrial genes at level 3, so the Exacto defaults
# (levels 1 and 2) would silently drop MT-ND5 — one of the vaccine targets.
GENE_LEVELS = ("1", "2", "3")
GENE_TYPES = ("protein_coding",)
