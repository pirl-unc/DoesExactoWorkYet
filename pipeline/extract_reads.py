"""Slice each long-read RNA-seq sample down to the vaccine-gene loci.

The three ONT dedup BAMs are 37, 67 and 53 GB and the PacBio Iso-Seq BAM is
1.8 GB. Backblaze serves them all with byte-range support and the BAM indexes
are small, so htslib can fetch just the blocks covering our regions — a couple
of minutes instead of 159 GB.

Volume is wildly uneven: the mitochondrial window alone holds ~1.3M reads, 88%
of everything in scope, while VPS13B's variant has 20. So reads are split in two
and each half is capped, with a fixed seed so runs stay reproducible:

  * **spanning reads** — their alignment covers a vaccine variant position.
    These are the only reads that can carry the mutation. Capped per variant,
    high enough that the cap cannot change whether a variant is callable.
  * **context reads** — elsewhere in the same gene. Interchangeable filler that
    helps RNA-Bloom2 extend transcripts. Capped per region.

Reads come out in their original (transcript-sense) orientation, which is what
both wf-single-cell and Iso-Seq put in and what minimap2's ``-uf`` expects on
the way back in.
"""

from __future__ import annotations

import gzip
import json
import random
import time
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pysam is only needed to actually read a BAM. Importing it
    import pysam  # eagerly would drag htslib into the site build and the unit
                   # tests, which have no business needing it.

from .build_reference import DOWNLOAD_DIR, REGIONS_JSON, download, load_variants
from .config import (
    CONTEXT_READS_PER_REGION,
    SPANNING_READS_PER_VARIANT,
    WORK_DIR,
    Sample,
    ensure_ca_bundle,
    samples_named,
)

READS_DIR = WORK_DIR / "reads"

# Streaming tens of thousands of BGZF blocks over HTTPS occasionally trips an
# HTTP/2 framing error or a connection reset. Seen once in testing, so it is
# worth surviving rather than failing a three-hour CI job.
REGION_FETCH_ATTEMPTS = 4
REGION_FETCH_BACKOFF_SECONDS = 5


def remote_size(url: str) -> int | None:
    """Byte size of a remote file, for the report. Best effort — never fatal."""
    import urllib.error
    import urllib.request

    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            length = response.headers.get("Content-Length")
            return int(length) if length else None
    except (urllib.error.URLError, ValueError, TimeoutError):
        return None


def load_regions() -> list[dict]:
    return json.loads(REGIONS_JSON.read_text())


def spanning_fastq(sample: Sample) -> Path:
    """Reads whose alignment covers a vaccine variant. The only ones that can
    carry a mutation, so the reads arm needs nothing else."""
    return READS_DIR / f"{sample.name}.spanning.fastq.gz"


def context_fastq(sample: Sample) -> Path:
    """Everything else in the gene — filler that lets RNA-Bloom2 extend
    transcripts. Capped per region."""
    return READS_DIR / f"{sample.name}.context.fastq.gz"


def assembly_inputs(sample: Sample) -> list[Path]:
    return [spanning_fastq(sample), context_fastq(sample)]


def stats_path(sample: Sample) -> Path:
    return READS_DIR / f"{sample.name}.extraction.json"


def variant_span(variant: dict) -> tuple[int, int]:
    """Reference bases a read must cover to be informative for this variant.

    For an SNV that is the single mutated base; for the VCF-style deletions on
    the portal it is the deleted bases, anchored at the shared leading base.
    """
    start = variant["pos"]
    end = variant["pos"] + max(len(variant["ref"]), 1) - 1
    return start, end


def _fastq_record(read: pysam.AlignedSegment) -> str | None:
    sequence = read.get_forward_sequence()
    qualities = read.get_forward_qualities()
    if not sequence or qualities is None:
        return None
    quality_string = "".join(chr(value + 33) for value in qualities)
    return f"@{read.query_name}\n{sequence}\n+\n{quality_string}\n"


def scan_region(
    bam,
    sample: Sample,
    region: dict,
    spans: list[tuple],
    already_seen: set[str],
) -> dict:
    """Sample one region's reads, retrying if the connection drops.

    Reading tens of thousands of BGZF blocks over HTTPS occasionally trips an
    HTTP/2 framing error or a reset, which htslib surfaces as an OSError
    mid-iteration. Everything is accumulated locally and only handed back once
    the region has been read end to end, so a retry starts from a clean slate
    rather than half a region's worth of state.
    """
    for attempt in range(1, REGION_FETCH_ATTEMPTS + 1):
        # Seeded per region, not per attempt, so a retry samples identically.
        rng = random.Random(f"{sample.name}:{region['chrom']}:{region['start']}")
        context_reservoir: list[str] = []
        spanning_reservoirs: dict[str, list[str]] = {
            variant["variant_id"]: [] for variant, _, _ in spans
        }
        spanning_seen: dict[str, int] = {key: 0 for key in spanning_reservoirs}
        names: set[str] = set()
        n_context_seen = 0

        try:
            for read in bam.fetch(region["chrom"], region["start"] - 1, region["end"]):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                name = read.query_name
                if name in already_seen or name in names:
                    # A read long enough to touch two neighbouring regions.
                    continue

                start = read.reference_start + 1
                end = read.reference_end or start
                covered = [
                    variant
                    for variant, span_start, span_end in spans
                    if start <= span_end and end >= span_start
                ]

                record = _fastq_record(read)
                if record is None:
                    continue
                names.add(name)

                if covered:
                    for variant in covered:
                        variant_id = variant["variant_id"]
                        spanning_seen[variant_id] += 1
                        reservoir = spanning_reservoirs[variant_id]
                        if len(reservoir) < SPANNING_READS_PER_VARIANT:
                            reservoir.append(record)
                        else:
                            slot = rng.randrange(spanning_seen[variant_id])
                            if slot < SPANNING_READS_PER_VARIANT:
                                reservoir[slot] = record
                else:
                    n_context_seen += 1
                    if len(context_reservoir) < CONTEXT_READS_PER_REGION:
                        context_reservoir.append(record)
                    else:
                        slot = rng.randrange(n_context_seen)
                        if slot < CONTEXT_READS_PER_REGION:
                            context_reservoir[slot] = record
        except (OSError, ValueError) as error:
            if attempt == REGION_FETCH_ATTEMPTS:
                raise
            wait = REGION_FETCH_BACKOFF_SECONDS * attempt
            print(
                f"    {region['chrom']}:{region['start']}-{region['end']} failed "
                f"({error}); retry {attempt}/{REGION_FETCH_ATTEMPTS - 1} in {wait}s"
            )
            time.sleep(wait)
            continue

        return {
            "spanning": spanning_reservoirs,
            "spanning_seen": spanning_seen,
            "context": context_reservoir,
            "context_seen": n_context_seen,
            "names": names,
        }

    raise RuntimeError("unreachable")


def extract(sample: Sample, regions: list[dict], variants: list[dict]) -> dict:
    """Write one sample's reads to two gzipped FASTQs, spanning and context."""
    bundle = ensure_ca_bundle()
    if bundle:
        print(f"  using CA bundle {bundle}")
    import pysam

    out_spanning = spanning_fastq(sample)
    out_context = context_fastq(sample)
    out_spanning.parent.mkdir(parents=True, exist_ok=True)

    # htslib wants the index beside the file; fetching it once locally avoids a
    # remote index read per region.
    index_path = download(sample.bai_url, DOWNLOAD_DIR / f"{sample.name}.bam.bai")
    bam_bytes = remote_size(sample.bam_url)

    variants_by_chrom: dict[str, list[dict]] = {}
    for variant in variants:
        variants_by_chrom.setdefault(variant["chrom"], []).append(variant)

    seen: set[str] = set()
    spanning_seen: dict[str, int] = {variant["variant_id"]: 0 for variant in variants}
    spanning_kept: dict[str, int] = {variant["variant_id"]: 0 for variant in variants}
    per_region: list[dict] = []
    n_spanning = 0
    n_context = 0
    n_bases = 0

    with pysam.AlignmentFile(
        sample.bam_url, "rb", index_filename=str(index_path)
    ) as bam, gzip.open(out_spanning, "wt") as sink, gzip.open(
        out_context, "wt"
    ) as context_sink:
        for region in regions:
            in_region = [
                variant
                for variant in variants_by_chrom.get(region["chrom"], [])
                if region["start"] <= variant["pos"] <= region["end"]
            ]
            spans = [(variant, *variant_span(variant)) for variant in in_region]

            scanned = scan_region(bam, sample, region, spans, seen)
            context_reservoir = scanned["context"]
            spanning_reservoirs = scanned["spanning"]
            n_context_seen = scanned["context_seen"]
            for variant_id, count in scanned["spanning_seen"].items():
                spanning_seen[variant_id] += count
            seen |= scanned["names"]

            # A read covering two variants sits in both reservoirs; write once.
            written: set[str] = set()
            region_spanning = 0
            for variant_id, reservoir in spanning_reservoirs.items():
                spanning_kept[variant_id] = len(reservoir)
                for record in reservoir:
                    if record in written:
                        continue
                    written.add(record)
                    sink.write(record)
                    n_bases += len(record.split("\n")[1])
                    region_spanning += 1
            n_spanning += region_spanning

            for record in context_reservoir:
                context_sink.write(record)
                n_bases += len(record.split("\n")[1])
            n_context += len(context_reservoir)

            per_region.append(
                {
                    "chrom": region["chrom"],
                    "start": region["start"],
                    "end": region["end"],
                    "genes": region["genes"],
                    "spanning_reads": region_spanning,
                    "context_reads_seen": n_context_seen,
                    "context_reads_kept": len(context_reservoir),
                }
            )
            print(
                f"  {sample.name} {','.join(region['genes']):<20} "
                f"spanning={region_spanning:>6,}  "
                f"context={len(context_reservoir):>6,}/{n_context_seen:,}"
            )

    n_reads = n_spanning + n_context
    stats = {
        "sample": sample.name,
        "timepoint": sample.timepoint,
        "platform": sample.platform,
        "bam_url": sample.bam_url,
        "bam_bytes": bam_bytes,
        "spanning_fastq": str(out_spanning),
        "context_fastq": str(out_context),
        "context_reads_per_region_cap": CONTEXT_READS_PER_REGION,
        "spanning_reads_per_variant_cap": SPANNING_READS_PER_VARIANT,
        "n_reads": n_reads,
        "n_spanning_reads": n_spanning,
        "n_context_reads": n_context,
        "n_bases": n_bases,
        "mean_read_length": round(n_bases / n_reads, 1) if n_reads else 0.0,
        "spanning_reads_by_variant": spanning_kept,
        "spanning_reads_seen_by_variant": spanning_seen,
        "regions": per_region,
    }
    stats_path(sample).write_text(json.dumps(stats, indent=2) + "\n")
    print(
        f"{sample.name}: {n_spanning:,} spanning + {n_context:,} context reads "
        f"-> {out_spanning.name}, {out_context.name}"
    )
    return stats


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        nargs="*",
        help=(
            "Which sequencing samples to extract, e.g. T1-ONT T1-PacBio "
            "(default: all)."
        ),
    )
    args = parser.parse_args()

    regions = load_regions()
    variants = load_variants()
    for sample in samples_named(args.samples):
        done = all(path.exists() for path in assembly_inputs(sample))
        if done and stats_path(sample).exists():
            print(f"{sample.name}: reusing {spanning_fastq(sample).parent}")
            continue
        extract(sample, regions, variants)


if __name__ == "__main__":
    main()
