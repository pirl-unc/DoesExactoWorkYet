"""Run the Exacto mutant-proteoform pipeline on each timepoint.

Two arms per timepoint:

``assembly``
    The pipeline as Exacto documents it — RNA-Bloom2 stitches the long reads
    into full-length transcripts, those are realigned, unspliced RNAs are
    dropped, and Exacto calls variants on what is left.

``reads``
    The same thing minus the assembler: reads go straight in as transcripts.
    Cheaper, and it isolates whether a miss is Exacto's or the assembler's.

Both arms then annotate the known vaccine mutations as the somatic DNA callset,
integrate them with Exacto's own RNA calls, and translate.

Every step's command line, exit status, duration and log lives in the run JSON,
so a failure shows up on the site as a failure rather than as a silent zero.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path

from .build_reference import (
    GENE_PROTEINS,
    MASKED_FASTA,
    SUBSET_GTF,
    load_variants,
)
from .config import (
    ARMS,
    GENCODE_RELEASE,
    GENE_LEVELS,
    GENE_TYPES,
    MINIMAP2_COMMON_FLAGS,
    MINIMAP2_PRESET,
    TIMEPOINTS,
    WORK_DIR,
    Timepoint,
)
from .extract_reads import assembly_inputs, spanning_fastq, stats_path

EXACTO_DIR = WORK_DIR / "exacto"

# Phred character stamped on every base of an assembled contig. 'I' is Q40,
# comfortably over call-rna-vars' default --min-average-base-quality of 30.
ASSEMBLY_BASE_QUALITY = "I"

ANNOTATION_ARGS = [
    "--reference-gene-annotation-file", str(SUBSET_GTF),
    "--reference-gene-annotation-source", "gencode",
    "--reference-gene-annotation-assembly", "hg38",
    "--reference-gene-annotation-version", f"v{GENCODE_RELEASE}",
]

LEVEL_ARGS = [
    "--gene-types", *GENE_TYPES,
    "--gene-levels", *GENE_LEVELS,
    "--transcript-types", *GENE_TYPES,
    "--transcript-levels", *GENE_LEVELS,
]


class StepFailed(RuntimeError):
    """A pipeline step exited non-zero; the arm cannot continue."""


class Runner:
    """Runs shell steps, keeping a structured record of each one."""

    def __init__(self, log_dir: Path, threads: int) -> None:
        self.log_dir = log_dir
        self.threads = threads
        self.steps: list[dict] = []
        log_dir.mkdir(parents=True, exist_ok=True)

    def run(
        self,
        name: str,
        command: list[str],
        stdout_path: Path | None = None,
    ) -> dict:
        log_path = self.log_dir / f"{len(self.steps):02d}_{name}.log"
        print(f"    [{name}] {' '.join(command[:6])} ...")
        started = time.monotonic()
        with open(log_path, "wb") as log:
            if stdout_path is not None:
                with open(stdout_path, "wb") as sink:
                    completed = subprocess.run(command, stdout=sink, stderr=log)
            else:
                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        elapsed = round(time.monotonic() - started, 1)

        step = {
            "name": name,
            "command": command,
            "returncode": completed.returncode,
            "seconds": elapsed,
            "log": str(log_path),
        }
        if completed.returncode != 0:
            step["log_tail"] = log_path.read_text(errors="replace")[-4000:]
        self.steps.append(step)

        if completed.returncode != 0:
            raise StepFailed(f"{name} exited {completed.returncode}; see {log_path}")
        return step


# --------------------------------------------------------------------------
# The known vaccine mutations, expressed as an Exacto somatic DNA callset
# --------------------------------------------------------------------------

SOMATIC_TSV_COLUMNS = [
    "variant_call_id",
    "chromosome_1",
    "position_1",
    "strand_1",
    "operation_1",
    "chromosome_2",
    "position_2",
    "strand_2",
    "operation_2",
    "variant_size",
    "variant_type",
    "sequence",
]


def as_graph_operation(variant: dict) -> tuple[int, int, int, str, str]:
    """Translate a VCF-style ref/alt into Exacto's breakpoint encoding.

    Exacto brackets a variant with the untouched bases either side: position_1
    is the last reference base kept before the edit ("D", downstream of it) and
    position_2 the first kept after ("U"). ``sequence`` is what goes between
    them, so a deletion carries an empty sequence.
    """
    pos, ref, alt = variant["pos"], variant["ref"], variant["alt"]

    if len(ref) == 1 and len(alt) == 1:
        return pos - 1, pos + 1, 1, "SNV", alt

    if len(ref) > len(alt) and ref.startswith(alt):
        # VCF-style deletion anchored on a shared leading base.
        deleted = len(ref) - len(alt)
        position_1 = pos + len(alt) - 1
        return position_1, position_1 + deleted + 1, deleted, "DEL", ""

    if len(alt) > len(ref) and alt.startswith(ref):
        inserted = alt[len(ref) :]
        position_1 = pos + len(ref) - 1
        return position_1, position_1 + 1, len(inserted), "INS", inserted

    if len(ref) == len(alt):
        return pos - 1, pos + len(ref), len(alt), "MNV", alt

    raise SystemExit(
        f"cannot encode {variant['gene']} {variant['chrom']}:{pos} {ref}>{alt}"
    )


def write_somatic_tsv(variants: list[dict], out_path: Path) -> dict[int, dict]:
    """Write the vaccine mutations as Exacto's somatic-DNA-variant TSV.

    Exacto normally gets this from ``call-somatic-dna-vars`` on tumour/normal
    long-read WGS. Sid's WGS is short-read, so the portal's curated calls stand
    in: this test asks whether Exacto can find these mutations *in the RNA* and
    build a mutant protein from them, not whether it can rediscover them in DNA.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    by_call_id: dict[int, dict] = {}

    with open(out_path, "w") as sink:
        sink.write("\t".join(SOMATIC_TSV_COLUMNS) + "\n")
        for index, variant in enumerate(variants, start=1):
            position_1, position_2, size, variant_type, sequence = as_graph_operation(
                variant
            )
            sink.write(
                "\t".join(
                    [
                        str(index),
                        variant["chrom"],
                        str(position_1),
                        "+",
                        "D",
                        variant["chrom"],
                        str(position_2),
                        "+",
                        "U",
                        str(size),
                        variant_type,
                        sequence,
                    ]
                )
                + "\n"
            )
            by_call_id[index] = variant

    return by_call_id


# --------------------------------------------------------------------------
# Arms
# --------------------------------------------------------------------------


def fasta_to_fastq(fasta: Path, fastq: Path, quality: str = ASSEMBLY_BASE_QUALITY) -> int:
    """Give assembled contigs a flat base quality so Exacto will read them.

    ``call-rna-vars`` indexes into the per-base quality array unconditionally
    (exacto-caller/src/structs/alignment.rs:229), so a BAM built from a FASTA —
    which is exactly what RNA-Bloom2 emits and what Exacto's own docs align —
    panics on every record. Exacto's bundled test BAMs carry a flat quality
    string, so this matches what the tool is built against. Consensus contigs
    have no meaningful per-base quality anyway.
    """
    records = 0
    with open(fasta) as source, open(fastq, "w") as sink:
        name = None
        chunks: list[str] = []

        def flush() -> None:
            nonlocal records
            if name is None:
                return
            sequence = "".join(chunks)
            sink.write(f"@{name}\n{sequence}\n+\n{quality * len(sequence)}\n")
            records += 1

        for line in source:
            if line.startswith(">"):
                flush()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        flush()
    return records


def count_fasta_records(fasta: Path) -> int:
    """Count records without pulling a few hundred MB of sequence into memory."""
    with open(fasta, "rb") as handle:
        return sum(1 for line in handle if line.startswith(b">"))


def find_assembly(outdir: Path) -> Path:
    """Locate RNA-Bloom2's transcript FASTA, whose name varies by version."""
    candidates = [
        "rnabloom.transcripts.fa",
        "rnabloom.longreads.assembly4.pol.fa",
        "rnabloom.longreads.assembly3.fa",
    ]
    for name in candidates:
        path = outdir / name
        if path.exists() and path.stat().st_size > 0:
            return path
    found = sorted(outdir.glob("rnabloom*.fa"))
    if found:
        return max(found, key=lambda item: item.stat().st_size)
    raise StepFailed(f"RNA-Bloom2 produced no transcripts in {outdir}")


def align(runner: Runner, query: Path, arm: str, out_bam: Path) -> None:
    sam = out_bam.with_suffix(".sam")
    mapped = out_bam.with_suffix(".mapped.bam")
    runner.run(
        "minimap2",
        [
            "minimap2",
            "-ax",
            MINIMAP2_PRESET[arm],
            *MINIMAP2_COMMON_FLAGS,
            "-t",
            str(runner.threads),
            str(MASKED_FASTA),
            str(query),
        ],
        stdout_path=sam,
    )
    # minimap2 emits unmapped records, and Exacto 0.4.6a1 unwraps the reference
    # sequence id without checking (exacto-core/src/common/bam.rs:892), so
    # remove-unspliced-rnas panics the moment it meets one. Dropping them is
    # standard hygiene anyway — they carry no variant information.
    runner.run(
        "samtools_drop_unmapped",
        ["samtools", "view", "-b", "-F", "4", "-o", str(mapped), str(sam)],
    )
    runner.run(
        "samtools_sort",
        ["samtools", "sort", "-@", str(runner.threads), "-o", str(out_bam), str(mapped)],
    )
    runner.run("samtools_index", ["samtools", "index", str(out_bam)])
    sam.unlink(missing_ok=True)
    mapped.unlink(missing_ok=True)


def drop_transcriptless_rna_calls(source: Path, dest: Path) -> int:
    """Copy the RNA callset without the rows that crash ``integrate-vars``.

    Exacto parses ``reference_transcript_ids`` with ``field.split(",")``, and in
    Rust ``"".split(",")`` yields one empty string rather than nothing. A call
    that matched no reference transcript therefore looks like a call against a
    transcript named "", takes the branch meant for annotated transcripts, and
    dies on ``get_transcript("").unwrap()``
    (exacto-integrator/src/algorithms/variant_integration.rs:98).

    Only the copy handed to ``integrate-vars`` is filtered — ``translate-structs``
    still sees every call. Returns how many rows were dropped so the loss shows
    up in the run record instead of vanishing.

    Lines are copied verbatim rather than re-serialised: Exacto writes an empty
    field as a quoted ``""`` and reads a bare empty cell back as null, which
    panics elsewhere (rna_variant_call_set.rs:236). Round-tripping through a CSV
    writer silently drops those quotes.
    """
    dropped = 0
    with open(source) as handle, open(dest, "w") as sink:
        header = handle.readline()
        sink.write(header)
        try:
            column = header.rstrip("\n").split("\t").index("reference_transcript_ids")
        except ValueError as error:
            raise SystemExit(f"{source} has no reference_transcript_ids column") from error

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= column or not fields[column].strip('"'):
                dropped += 1
                continue
            sink.write(line)
    return dropped


def count_alignments(bam: Path) -> int:
    result = subprocess.run(
        ["samtools", "view", "-c", "-F", "0x900", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout.strip() or 0)


def run_arm(
    timepoint: Timepoint,
    arm: str,
    variants: list[dict],
    threads: int,
) -> dict:
    """One timepoint through one arm of the pipeline."""
    out_dir = EXACTO_DIR / timepoint.name / arm
    out_dir.mkdir(parents=True, exist_ok=True)
    runner = Runner(out_dir / "logs", threads)
    prefix = f"{timepoint.name}_{arm}"

    result: dict = {
        "timepoint": timepoint.name,
        "arm": arm,
        "status": "ok",
        "outputs": {},
        "counts": {},
    }

    stats = json.loads(stats_path(timepoint).read_text())
    query_fasta: Path | None = None

    try:
        if arm == "assembly":
            assembly_dir = out_dir / "rnabloom"
            if assembly_dir.exists():
                shutil.rmtree(assembly_dir)
            runner.run(
                "rnabloom",
                [
                    "rnabloom",
                    "-long", *(str(path) for path in assembly_inputs(timepoint)),
                    "-o", str(assembly_dir),
                    "-t", str(threads),
                    "-chimera",
                    "-f",
                ],
            )
            query_fasta = out_dir / f"{prefix}.transcripts.fa"
            shutil.copyfile(find_assembly(assembly_dir), query_fasta)
            query = out_dir / f"{prefix}.transcripts.fq"
            result["counts"]["query_sequences"] = fasta_to_fastq(query_fasta, query)
        else:
            # Only the variant-spanning reads. Without an assembler each read is
            # its own transcript, so a read that touches no vaccine variant can
            # never produce one of the mutant proteins under test — it would just
            # cost Exacto time. minimap2 carries their real base qualities into
            # the BAM, which Exacto needs.
            query = spanning_fastq(timepoint)
            result["counts"]["query_sequences"] = stats["n_spanning_reads"]

        aligned_bam = out_dir / f"{prefix}.aligned.bam"
        align(runner, query, arm, aligned_bam)
        result["counts"]["aligned"] = count_alignments(aligned_bam)

        if arm == "assembly":
            # Only the documented assembly pipeline filters unspliced RNAs;
            # applied to raw reads it would throw away most of the data.
            filtered_bam = out_dir / f"{prefix}.filtered.bam"
            try:
                runner.run(
                    "remove_unspliced_rnas",
                    [
                        "exacto", "remove-unspliced-rnas",
                        "--bam-file", str(aligned_bam),
                        "--bam-bai-file", str(aligned_bam) + ".bai",
                        "--fasta-file", str(query_fasta),  # the contigs, not the genome
                        *ANNOTATION_ARGS,
                        *LEVEL_ARGS,
                        "--output-bam-file", str(filtered_bam),
                        "--output-bam-bai-file", str(filtered_bam) + ".bai",
                        "--output-fasta-file", str(out_dir / f"{prefix}.filtered.fa"),
                        "--num-threads", str(threads),
                    ],
                )
            except StepFailed:
                # Exacto 0.4.6a1 writes the kept records in hash order while
                # stamping SO:coordinate on the header, so its own indexing step
                # rejects the file it just wrote. The filtering itself finished —
                # sort what it produced and carry on.
                if not filtered_bam.exists() or not count_alignments(filtered_bam):
                    raise
                result["workarounds"] = result.get("workarounds", [])
                result["workarounds"].append("sorted the unspliced-filter output")
                sorted_bam = out_dir / f"{prefix}.filtered.sorted.bam"
                runner.run(
                    "sort_filter_output",
                    ["samtools", "sort", "-@", str(threads),
                     "-o", str(sorted_bam), str(filtered_bam)],
                )
                sorted_bam.replace(filtered_bam)
                runner.run("index_filter_output", ["samtools", "index", str(filtered_bam)])

            variant_calling_bam = filtered_bam
            result["counts"]["after_unspliced_filter"] = count_alignments(filtered_bam)
        else:
            variant_calling_bam = aligned_bam

        rna_dir = out_dir / "rna_vars"
        runner.run(
            "call_rna_vars",
            [
                "exacto", "call-rna-vars",
                "--bam-file", str(variant_calling_bam),
                "--bam-bai-file", str(variant_calling_bam) + ".bai",
                "--reference-genome-fasta-file", str(MASKED_FASTA),
                *ANNOTATION_ARGS,
                "--output-dir", str(rna_dir),
                "--output-prefix", prefix,
                "--num-threads", str(threads),
            ],
        )
        rna_calls = rna_dir / f"{prefix}_exacto_rna_variant_calls.tsv"
        structures = rna_dir / f"{prefix}_exacto_transcript_structures.tsv"

        somatic_tsv = out_dir / f"{prefix}.vaccine_variants.tsv"
        write_somatic_tsv(variants, somatic_tsv)

        annotated_tsv = out_dir / f"{prefix}.vaccine_variants.annotated.tsv"
        runner.run(
            "annotate_vars",
            [
                "exacto", "annotate-vars",
                "--tsv-file", str(somatic_tsv),
                *ANNOTATION_ARGS,
                *LEVEL_ARGS,
                "--output-tsv-file", str(annotated_tsv),
                "--num-threads", str(threads),
            ],
        )

        integrable_calls = out_dir / f"{prefix}.rna_variant_calls.integrable.tsv"
        dropped = drop_transcriptless_rna_calls(rna_calls, integrable_calls)
        if dropped:
            result["counts"]["rna_calls_without_reference_transcript"] = dropped
            result.setdefault("workarounds", []).append(
                f"withheld {dropped:,} RNA calls with no reference transcript from "
                "integrate-vars"
            )

        integrated_tsv = out_dir / f"{prefix}.integrated.tsv"
        runner.run(
            "integrate_vars",
            [
                "exacto", "integrate-vars",
                "--annotated-dna-vars-tsv-file", str(annotated_tsv),
                "--rna-vars-tsv-file", str(integrable_calls),
                *ANNOTATION_ARGS,
                # Exacto defaults to linking a DNA variant to any RNA variant
                # within 10 kb of a transcript boundary, or 100 kb intergenically,
                # which produces mostly-spurious pairings — 19 of 3,359 were exact
                # at the defaults. Keep the small exon offset, which exists for
                # indel placement wobble, and drop the rest.
                "--max-exon-offset", "2",
                "--max-transcript-boundary-offset", "0",
                "--max-intergenic-distance", "0",
                "--output-tsv-file", str(integrated_tsv),
                "--num-threads", str(threads),
            ],
        )

        primary_tsv = out_dir / f"{prefix}.primary_structures.tsv"
        primary_fasta = out_dir / f"{prefix}.primary_structures.fasta"
        runner.run(
            "translate_structs",
            [
                "exacto", "translate-structs",
                "--transcript-structures-tsv-file", str(structures),
                "--rna-variant-calls-tsv-file", str(rna_calls),
                "--integrated-variants-tsv-file", str(integrated_tsv),
                "--strategy", "longest_orf",
                "--output-tsv-file", str(primary_tsv),
                "--output-fasta-file", str(primary_fasta),
                "--num-threads", str(threads),
            ],
        )

        result["outputs"] = {
            "rna_variant_calls": str(rna_calls),
            "transcript_structures": str(structures),
            "somatic_variants": str(somatic_tsv),
            "annotated_variants": str(annotated_tsv),
            "integrated_variants": str(integrated_tsv),
            "primary_structures_tsv": str(primary_tsv),
            "primary_structures_fasta": str(primary_fasta),
        }

        # Neoantigen candidates. Downstream of the question this repo asks, so a
        # failure here is recorded but does not sink the arm.
        peptides_tsv = out_dir / f"{prefix}.peptide_variants.tsv"
        peptides_fasta = out_dir / f"{prefix}.peptide_variants.fasta"
        try:
            runner.run(
                "call_peptide_vars",
                [
                    "exacto", "call-peptide-vars",
                    "--primary-structures-tsv-file", str(primary_tsv),
                    "--reference-fasta-file", str(GENE_PROTEINS),
                    "--output-tsv-file", str(peptides_tsv),
                    "--output-fasta-file", str(peptides_fasta),
                    "--min-k", "8",
                    "--max-k", "11",
                    "--num-threads", str(threads),
                ],
            )
            result["outputs"]["peptide_variants"] = str(peptides_tsv)
            result["outputs"]["peptide_variants_fasta"] = str(peptides_fasta)
        except StepFailed as error:
            result["peptide_error"] = str(error)

    except StepFailed as error:
        result["status"] = "failed"
        result["error"] = str(error)
        print(f"    !! {timepoint.name}/{arm}: {error}")

    result["steps"] = runner.steps
    result["seconds"] = round(sum(step["seconds"] for step in runner.steps), 1)
    (out_dir / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    import argparse
    import os

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 2)
    parser.add_argument("--timepoints", nargs="*", default=[tp.name for tp in TIMEPOINTS])
    parser.add_argument("--arms", nargs="*", default=list(ARMS))
    args = parser.parse_args()

    variants = load_variants()
    runs = []
    for timepoint in TIMEPOINTS:
        if timepoint.name not in args.timepoints:
            continue
        for arm in args.arms:
            print(f"== {timepoint.name} / {arm} ==")
            runs.append(run_arm(timepoint, arm, variants, args.threads))

    failures = [run for run in runs if run["status"] != "ok"]
    print(f"\n{len(runs) - len(failures)}/{len(runs)} arms completed")
    for failure in failures:
        print(f"  FAILED {failure['timepoint']}/{failure['arm']}: {failure['error']}")


def collect_runs(timepoints: list[str] | None = None) -> list[dict]:
    """Gather the per-arm run.json files on disk.

    CI runs the timepoints as separate jobs, so the runs are assembled from what
    each one left behind rather than from a single in-process list. Only the
    canonical ``<timepoint>/<arm>/`` paths count — anything else under work/ is
    somebody's leftovers, not a result.
    """
    known_timepoints = {timepoint.name for timepoint in TIMEPOINTS}
    runs = []
    for run_path in sorted(EXACTO_DIR.glob("*/*/run.json")):
        arm_dir = run_path.parent
        if arm_dir.parent.name not in known_timepoints or arm_dir.name not in ARMS:
            continue
        run = json.loads(run_path.read_text())
        if timepoints and run["timepoint"] not in timepoints:
            continue
        runs.append(run)
    return runs


if __name__ == "__main__":
    main()
