"""Carve a small, coordinate-preserving reference around the vaccine genes.

The three ONT BAMs total ~157 GB and hg38 is another 3 GB, so nothing here
downloads a whole genome. Instead we:

  * find the GENCODE gene body for each vaccine mutation,
  * pull just those slices of hg38 over HTTP byte ranges,
  * write them back at their true hg38 offsets inside otherwise all-N
    chromosomes.

Coordinates therefore stay identical to hg38 — the GTF, the variant table and
every Exacto output line up without any offset arithmetic — while the sequence
that actually has to be indexed shrinks to a few megabases. minimap2 ignores
N runs when it collects minimizers, so the masked genome indexes in seconds.
"""

from __future__ import annotations

import gzip
import json
import subprocess
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .config import (
    GENCODE_GTF_URL,
    GENCODE_PROTEINS_URL,
    GENE_FLANK_BP,
    HG38_FASTA_URL,
    RESULTS_DIR,
    WORK_DIR,
)
from .fetch_osteosarc import USER_AGENT

REFERENCE_DIR = WORK_DIR / "reference"
DOWNLOAD_DIR = WORK_DIR / "downloads"

# Left uncompressed on purpose: Exacto queries the reference one base at a
# time, and on a bgzipped FASTA every one of those queries pays for a BGZF block
# decompression. Measured 110 s vs 74 s for the same 1,179 reads. The file is
# ~1.7 GB of mostly N, which is cheap next to that.
MASKED_FASTA = REFERENCE_DIR / "vaccine_genes.hg38.fa"
SUBSET_GTF = REFERENCE_DIR / "vaccine_genes.gencode.gtf.gz"
GENE_PROTEINS = REFERENCE_DIR / "vaccine_genes.proteins.fa"
REGIONS_JSON = REFERENCE_DIR / "regions.json"

FASTA_LINE_WIDTH = 60

# Remote reads of hg38 and of the ONT BAMs go over HTTPS to Backblaze, which
# occasionally drops a long connection with an HTTP/2 framing error.
FETCH_ATTEMPTS = 3
FETCH_BACKOFF_SECONDS = 5

# Windows are grown until every transcript overlapping them fits entirely on
# real sequence. Transcripts can chain one window into the next, so this is
# iterative; a runaway chain is a bug, not something to keep expanding through.
TRANSCRIPT_EXPANSION_ROUNDS = 6


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    genes: tuple[str, ...]

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


def download(url: str, dest: Path) -> Path:
    """Fetch ``url`` to ``dest`` unless it is already there."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"downloading {url}")
    partial = dest.with_suffix(dest.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=600) as response, open(
        partial, "wb"
    ) as handle:
        while chunk := response.read(1 << 20):
            handle.write(chunk)
    partial.rename(dest)
    return dest


def load_variants() -> list[dict]:
    payload = json.loads((RESULTS_DIR / "vaccine_variants.json").read_text())
    return payload["variants"]


def transcript_spans(gtf_path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Every GENCODE transcript's full genomic span, by chromosome."""
    spans: dict[str, list[tuple[int, int, str]]] = {}
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "transcript":
                continue
            attributes = fields[8]
            transcript_id = attributes.split('transcript_id "', 1)[1].split('"', 1)[0]
            spans.setdefault(fields[0], []).append(
                (int(fields[3]), int(fields[4]), transcript_id)
            )
    for chrom_spans in spans.values():
        chrom_spans.sort()
    return spans


def gene_spans(gene_names: set[str], gtf_path: Path) -> dict[str, tuple[str, int, int]]:
    """Gene-body span for each requested gene symbol."""
    spans: dict[str, tuple[str, int, int]] = {}
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "gene":
                continue
            attributes = fields[8]
            marker = 'gene_name "'
            if marker not in attributes:
                continue
            name = attributes.split(marker, 1)[1].split('"', 1)[0]
            if name not in gene_names:
                continue
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            if name in spans:
                # Same symbol on more than one contig (PAR genes, patches):
                # keep the copy on the contig we already chose.
                previous = spans[name]
                if previous[0] != chrom:
                    continue
                start = min(start, previous[1])
                end = max(end, previous[2])
            spans[name] = (chrom, start, end)
    return spans


def build_regions(variants: list[dict], gtf_path: Path) -> list[Region]:
    """One padded, merged interval per gene (or cluster of adjacent genes)."""
    wanted = {variant["gene"] for variant in variants}
    spans = gene_spans(wanted, gtf_path)

    missing = sorted(wanted - set(spans))
    if missing:
        raise SystemExit(f"genes absent from the GENCODE GTF: {', '.join(missing)}")

    intervals: list[tuple[str, int, int, str]] = []
    for gene, (chrom, start, end) in spans.items():
        intervals.append((chrom, max(1, start - GENE_FLANK_BP), end + GENE_FLANK_BP, gene))

    # Every variant must land inside its gene's padded interval; if the portal
    # and GENCODE disagree about a locus we want to hear about it now.
    for variant in variants:
        chrom, start, end, _ = next(
            item for item in intervals if item[3] == variant["gene"]
        )
        if chrom != variant["chrom"] or not start <= variant["pos"] <= end:
            raise SystemExit(
                f"{variant['gene']} {variant['chrom']}:{variant['pos']} falls outside "
                f"its GENCODE gene body {chrom}:{start}-{end}"
            )

    merged = _merge(intervals)

    # Any transcript Exacto considers must lie entirely on real sequence.
    # call-rna-vars rebuilds a candidate transcript's sequence base by base and
    # panics on anything that is not A/C/G/T (exacto-caller/src/structs/
    # reference_transcript_sequence.rs:176), so a transcript poking out of a
    # window into the N padding takes the read with it. Grow the windows until
    # every overlapping transcript fits; transcripts can chain, so repeat.
    transcripts = transcript_spans(gtf_path)
    for _ in range(TRANSCRIPT_EXPANSION_ROUNDS):
        grown: list[tuple[str, int, int, str]] = []
        changed = False
        for region in merged:
            start, end = region.start, region.end
            for tx_start, tx_end, _ in transcripts.get(region.chrom, []):
                if tx_start > end:
                    break
                if tx_end < start:
                    continue
                if tx_start < start or tx_end > end:
                    start, end = min(start, tx_start), max(end, tx_end)
                    changed = True
            grown.append((region.chrom, start, end, region.genes))
        merged = _merge(
            [
                (chrom, start, end, gene)
                for chrom, start, end, genes in grown
                for gene in genes
            ]
        )
        if not changed:
            break
    else:
        raise SystemExit(
            "transcript expansion did not settle — a transcript chain is walking "
            "the windows across the chromosome"
        )

    return merged


def _merge(intervals: list[tuple[str, int, int, str]]) -> list[Region]:
    """Collapse overlapping padded gene intervals into regions."""
    ordered = sorted(intervals, key=lambda item: (item[0], item[1]))
    merged: list[Region] = []
    for chrom, start, end, gene in ordered:
        if merged and merged[-1].chrom == chrom and start <= merged[-1].end:
            previous = merged[-1]
            merged[-1] = Region(
                chrom,
                previous.start,
                max(previous.end, end),
                tuple(sorted(set(previous.genes) | {gene})),
            )
        else:
            merged.append(Region(chrom, start, end, (gene,)))
    return merged


def contig_lengths(fai_path: Path) -> dict[str, int]:
    lengths = {}
    for line in fai_path.read_text().splitlines():
        name, length = line.split("\t")[:2]
        lengths[name] = int(length)
    return lengths


def _parse_faidx(payload: bytes) -> dict[str, bytes]:
    sequences: dict[str, bytes] = {}
    locus = None
    chunks: list[bytes] = []
    for line in payload.split(b"\n"):
        if line.startswith(b">"):
            if locus is not None:
                sequences[locus] = b"".join(chunks)
            locus = line[1:].decode().strip()
            chunks = []
        elif line:
            chunks.append(line.upper())
    if locus is not None:
        sequences[locus] = b"".join(chunks)
    return sequences


def _faidx(fasta_url: str, loci: list[str]) -> dict[str, bytes]:
    result = subprocess.run(
        ["samtools", "faidx", fasta_url, *loci], check=True, capture_output=True
    )
    return _parse_faidx(result.stdout)


def fetch_region_sequences(regions: list[Region], fasta_url: str) -> dict[str, bytes]:
    """Pull every hg38 slice we need over HTTP byte ranges.

    One batched call is far quicker than 37, but a long multi-region read of the
    remote FASTA occasionally trips an HTTP/2 framing error. Retry the batch a
    couple of times, then fall back to fetching each region on its own
    connection, which is slower but much harder to knock over.
    """
    loci = [region.locus for region in regions]

    for attempt in range(1, FETCH_ATTEMPTS + 1):
        try:
            sequences = _faidx(fasta_url, loci)
        except subprocess.CalledProcessError as error:
            print(f"  batched hg38 read failed (attempt {attempt}): {error}")
        else:
            if not [locus for locus in loci if locus not in sequences]:
                return sequences
            print(f"  batched hg38 read came back short (attempt {attempt})")
        time.sleep(FETCH_BACKOFF_SECONDS * attempt)

    print("  falling back to one request per region")
    sequences = {}
    for locus in loci:
        for attempt in range(1, FETCH_ATTEMPTS + 1):
            try:
                sequences.update(_faidx(fasta_url, [locus]))
                break
            except subprocess.CalledProcessError:
                if attempt == FETCH_ATTEMPTS:
                    raise
                time.sleep(FETCH_BACKOFF_SECONDS * attempt)

    missing = [locus for locus in loci if locus not in sequences]
    if missing:
        raise SystemExit(f"hg38 returned nothing for: {', '.join(missing)}")
    return sequences


def write_masked_fasta(
    regions: list[Region],
    fai_path: Path,
    fasta_url: str,
    out_path: Path,
) -> None:
    """Write chromosomes that are all N except for the regions of interest.

    Each chromosome is truncated just past its last region — everything before
    that keeps its real hg38 offset, which is all we need.
    """
    lengths = contig_lengths(fai_path)
    sequences = fetch_region_sequences(regions, fasta_url)

    by_chrom: dict[str, list[Region]] = {}
    for region in regions:
        by_chrom.setdefault(region.chrom, []).append(region)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    total_real_bp = 0
    with open(out_path, "wb") as sink:
        for chrom in sorted(by_chrom, key=lambda name: (len(name), name)):
            chrom_regions = sorted(by_chrom[chrom], key=lambda item: item.start)
            length = min(lengths[chrom], chrom_regions[-1].end)
            sequence = bytearray(b"N" * length)
            for region in chrom_regions:
                end = min(region.end, length)
                expected = end - region.start + 1
                payload = sequences[region.locus]
                if len(payload) < expected:
                    raise SystemExit(
                        f"short read from hg38 for {region.locus}: "
                        f"got {len(payload)} bases, expected {expected}"
                    )
                sequence[region.start - 1 : end] = payload[:expected]
                total_real_bp += expected

            sink.write(f">{chrom}\n".encode())
            for offset in range(0, length, FASTA_LINE_WIDTH):
                sink.write(sequence[offset : offset + FASTA_LINE_WIDTH])
                sink.write(b"\n")
            print(f"  {chrom}: {length:,} bp, {len(chrom_regions)} region(s)")

    subprocess.run(["samtools", "faidx", str(out_path)], check=True)
    print(f"masked reference: {out_path} ({total_real_bp:,} bp of real sequence)")


def verify_transcripts_contained(regions: list[Region], gtf_path: Path) -> int:
    """No transcript in the subset may run off the end of a region.

    If one does, Exacto will query masked bases while rebuilding its sequence
    and drop every read that matched it — silently losing variants rather than
    failing loudly. Cheap to check, so check it.
    """
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for region in regions:
        by_chrom.setdefault(region.chrom, []).append((region.start, region.end))

    total = 0
    escaping = []
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "transcript":
                continue
            total += 1
            start, end = int(fields[3]), int(fields[4])
            spans = by_chrom.get(fields[0], [])
            if not any(
                span_start <= start and end <= span_end for span_start, span_end in spans
            ):
                escaping.append(f"{fields[0]}:{start}-{end}")

    if escaping:
        raise SystemExit(
            f"{len(escaping)} transcript(s) extend past the unmasked windows, e.g. "
            f"{', '.join(escaping[:3])} — Exacto would hit N and drop reads"
        )
    return total


def write_subset_gtf(regions: list[Region], gtf_path: Path, out_path: Path) -> None:
    """Keep only GTF features overlapping the regions of interest."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for region in regions:
        by_chrom.setdefault(region.chrom, []).append((region.start, region.end))

    kept = 0
    with gzip.open(gtf_path, "rt") as source, gzip.open(out_path, "wt") as sink:
        for line in source:
            if line.startswith("#"):
                sink.write(line)
                continue
            fields = line.split("\t", 5)
            spans = by_chrom.get(fields[0])
            if not spans:
                continue
            start, end = int(fields[3]), int(fields[4])
            if any(start <= span_end and end >= span_start for span_start, span_end in spans):
                sink.write(line)
                kept += 1
    print(f"subset GTF: {out_path} ({kept:,} features)")


def write_gene_proteins(
    variants: list[dict], proteome_path: Path, out_path: Path
) -> None:
    """GENCODE translations for the genes under test.

    Used as the wild-type background for ``call-peptide-vars``; keeping it to
    the tested genes means a "novel" peptide is one absent from the gene's own
    reference proteins, which is the comparison this test cares about.
    """
    wanted = {variant["gene"] for variant in variants}
    written = 0
    with gzip.open(proteome_path, "rt") as source, open(out_path, "w") as sink:
        keep = False
        for line in source:
            if line.startswith(">"):
                # GENCODE translation headers are pipe-delimited with the gene
                # symbol in field 6.
                parts = line[1:].split("|")
                keep = len(parts) > 6 and parts[6] in wanted
                if keep:
                    written += 1
            if keep:
                sink.write(line)
    print(f"reference proteins: {out_path} ({written} sequences)")
    if written == 0:
        raise SystemExit("no reference proteins matched the vaccine genes")


def main() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    gtf_path = download(GENCODE_GTF_URL, DOWNLOAD_DIR / "gencode.annotation.gtf.gz")
    proteome_path = download(
        GENCODE_PROTEINS_URL, DOWNLOAD_DIR / "gencode.pc_translations.fa.gz"
    )
    fai_path = download(HG38_FASTA_URL + ".fai", DOWNLOAD_DIR / "hg38.fasta.fai")

    variants = load_variants()
    regions = build_regions(variants, gtf_path)
    print(f"{len(regions)} regions covering {len(variants)} vaccine variants")

    write_masked_fasta(regions, fai_path, HG38_FASTA_URL, MASKED_FASTA)
    write_subset_gtf(regions, gtf_path, SUBSET_GTF)
    n_transcripts = verify_transcripts_contained(regions, SUBSET_GTF)
    print(f"all {n_transcripts:,} subset transcripts sit on real sequence")
    write_gene_proteins(variants, proteome_path, GENE_PROTEINS)

    REGIONS_JSON.write_text(
        json.dumps(
            [
                {
                    "chrom": region.chrom,
                    "start": region.start,
                    "end": region.end,
                    "genes": list(region.genes),
                }
                for region in regions
            ],
            indent=2,
        )
        + "\n"
    )
    print(f"regions: {REGIONS_JSON}")


if __name__ == "__main__":
    main()
