"""Run the Exacto mutant-proteoform pipeline on each sample.

Two arms per sample:

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
    SAMPLES,
    WORK_DIR,
    Sample,
    arms_for,
    minimap2_preset,
    samples_named,
)
from .extract_reads import (
    SYNTHETIC_BASE_QUALITY,
    assembly_inputs,
    reads_arm_fastq,
    stats_path,
)
from .methods import METHODS_BY_NAME

EXACTO_DIR = WORK_DIR / "exacto"

# RNA-Bloom2 names its polished assembly and the read-to-contig alignment it
# is filtered against.
ASSEMBLY_NAME = "rnabloom.longreads.assembly4.pol.fa"
PAF_NAME = "rnabloom.longreads.assembly3.map.paf.gz"

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
                    completed = subprocess.run(
                        command, stdout=sink, stderr=log, check=False
                    )
            else:
                completed = subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT, check=False
                )
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


def count_fasta_records(fasta: Path) -> int:
    """Count records without pulling a few hundred MB of sequence into memory."""
    with open(fasta, "rb") as handle:
        return sum(1 for line in handle if line.startswith(b">"))


def find_assembly(outdir: Path) -> Path:
    """Locate RNA-Bloom2's transcript FASTA, whose name varies by version."""
    candidates = [ASSEMBLY_NAME, "rnabloom.transcripts.fa"]
    for name in candidates:
        path = outdir / name
        if path.exists() and path.stat().st_size > 0:
            return path
    found = sorted(outdir.glob("rnabloom*.fa"))
    if found:
        return max(found, key=lambda item: item.stat().st_size)
    raise StepFailed(f"RNA-Bloom2 produced no transcripts in {outdir}")


def align(
    runner: Runner, query: Path, sample: Sample, arm: str, out_bam: Path
) -> None:
    sam = out_bam.with_suffix(".sam")
    mapped = out_bam.with_suffix(".mapped.bam")
    runner.run(
        "minimap2",
        [
            "minimap2",
            "-ax",
            minimap2_preset(sample.platform, arm),
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


def fasta_to_fastq(source: Path, dest: Path, quality: int) -> int:
    """Give assembled contigs a flat quality so call-rna-vars will read them.

    Nexus does this for RNA-Bloom2 output inside its own filter. rnaSPAdes has
    no such wrapper, so the same trick is applied here — and it is the same
    fabrication, disclosed the same way.
    """
    written = 0
    with open(source) as handle, open(dest, "w") as sink:
        name, chunks = None, []

        def flush() -> None:
            nonlocal written
            if name is None:
                return
            sequence = "".join(chunks)
            if not sequence:
                return
            sink.write(f"@{name}\n{sequence}\n+\n{chr(quality + 33) * len(sequence)}\n")
            written += 1

        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                name, chunks = line[1:].split()[0], []
            elif name is not None:
                chunks.append(line.strip())
        flush()
    return written


def run_rnaspades(runner: Runner, sample, out_dir: Path, threads: int) -> Path:
    """Assemble short reads into transcript-like contigs.

    Exacto is a long-read tool and says so, but it never asks where a sequence
    came from — it wants transcripts. rnaSPAdes turns short reads into exactly
    that, which is the only honest way to ask whether the long-read requirement
    is about the reads themselves or about the assembly they make possible.

    Run in single-end mode. The mates are extracted independently by region, so
    a pair is frequently split across the spanning/context files or has one mate
    outside the window entirely; feeding that to -1/-2 as though it were an
    intact library would be worse than not claiming pairing at all.
    """
    assembly_dir = out_dir / "rnaspades"
    if assembly_dir.exists():
        shutil.rmtree(assembly_dir)
    runner.run(
        "rnaspades",
        [
            "rnaspades.py",
            *(arg for path in assembly_inputs(sample) for arg in ("-s", str(path))),
            "-o", str(assembly_dir),
            "-t", str(threads),
            "-m", "12",
        ],
    )
    contigs = assembly_dir / "transcripts.fasta"
    if not contigs.exists() or not contigs.stat().st_size:
        raise StepFailed(f"rnaSPAdes produced no transcripts in {assembly_dir}")
    return contigs


def run_isoncorrect(runner: Runner, sample, out_dir: Path, threads: int) -> Path:
    """Cluster long reads by shared structure, then error-correct within cluster.

    The point is to remove basecalling indels without removing the variant.
    Assembly does both: it averages over every read at a locus, so a subclonal
    allele is outvoted and disappears. isONclust groups reads that share
    transcript structure and isONcorrect polishes each read against the others
    in its own group, emitting one corrected read per input read — so a minority
    allele keeps its own read and nothing outvotes it.

    Reference-free throughout, which is why this is preferred to anchoring on an
    annotated transcript: a novel junction defines its own cluster instead of
    being measured against a reference that does not contain it.
    """
    clusters = out_dir / "isonclust"
    if clusters.exists():
        shutil.rmtree(clusters)
    mode = "--isoseq" if sample.platform == "PacBio" else "--ont"
    runner.run(
        "isonclust",
        [
            "isONclust", mode,
            "--fastq", str(reads_arm_fastq(sample)),
            "--outfolder", str(clusters),
            "--t", str(threads),
        ],
    )
    per_cluster = clusters / "fastq_files"
    runner.run(
        "isonclust_write_fastq",
        [
            "isONclust", "write_fastq",
            "--clusters", str(clusters / "final_clusters.tsv"),
            "--fastq", str(reads_arm_fastq(sample)),
            "--outfolder", str(per_cluster),
            "--N", "1",
        ],
    )
    corrected_dir = out_dir / "isoncorrect"
    if corrected_dir.exists():
        shutil.rmtree(corrected_dir)
    runner.run(
        "isoncorrect",
        [
            "run_isoncorrect",
            "--fastq_folder", str(per_cluster),
            "--outfolder", str(corrected_dir),
            "--t", str(threads),
        ],
    )
    merged = out_dir / f"{sample.name}_corrected.fastq"
    parts = sorted(corrected_dir.glob("*/corrected_reads.fastq"))
    if not parts:
        raise StepFailed(f"isONcorrect produced no corrected reads in {corrected_dir}")
    merged.write_text("".join(part.read_text() for part in parts))
    return merged


def run_isonform(runner: Runner, sample, out_dir: Path, threads: int,
                 corrected: Path) -> Path:
    """Assemble each corrected cluster into isoforms.

    The distinction from RNA-Bloom2 is where the collapsing happens. A global
    assembler pools every read at a locus, so a subclonal allele is a minority
    inside its own contig and is averaged out. isONform assembles within a
    cluster that isONclust already separated by transcript structure, so the
    averaging is over reads that agree, not over reads that differ.
    """
    per_cluster = out_dir / "isonclust" / "fastq_files"
    forms = out_dir / "isonform"
    if forms.exists():
        shutil.rmtree(forms)
    runner.run(
        "isonform",
        [
            "isONform_parallel",
            "--fastq_folder", str(per_cluster),
            "--outfolder", str(forms),
            "--t", str(threads),
            "--exact_instance_limit", "50",
            "--split_wrt_batches",
        ],
    )
    merged = out_dir / f"{sample.name}_isonform.fastq"
    parts = sorted(forms.rglob("transcripts.fastq"))
    if not parts:
        raise StepFailed(f"isONform produced no transcripts in {forms}")
    merged.write_text("".join(part.read_text() for part in parts))
    return merged


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
    sample: Sample,
    arm: str,
    variants: list[dict],
    threads: int,
) -> dict:
    """One sample through one method."""
    method = METHODS_BY_NAME[arm]
    out_dir = EXACTO_DIR / sample.name / arm
    out_dir.mkdir(parents=True, exist_ok=True)
    runner = Runner(out_dir / "logs", threads)
    prefix = f"{sample.name}_{arm}"

    result: dict = {
        "sample": sample.name,
        "timepoint": sample.timepoint,
        "platform": sample.platform,
        "assay": sample.assay,
        "label": sample.label,
        "arm": arm,
        "method": {
            "name": method.name,
            "family": method.family,
            "label": method.label,
            "tool": method.tool,
            "params": method.params,
        },
        "status": "ok",
        "outputs": {},
        "counts": {},
    }

    stats = json.loads(stats_path(sample).read_text())
    query_fasta: Path | None = None

    try:
        if method.family == "assembly" and sample.read_type == "short":
            # Short reads reach Exacto only as contigs. rnaSPAdes emits FASTA
            # with no per-base quality, and call-rna-vars panics without one,
            # so the same flat score Nexus writes for RNA-Bloom2 output is
            # written here.
            query_fasta = run_rnaspades(runner, sample, out_dir, threads)
            query = out_dir / f"{prefix}.transcripts.fastq"
            n = fasta_to_fastq(query_fasta, query, SYNTHETIC_BASE_QUALITY)
            result["counts"]["assembled_transcripts"] = n
            result["counts"]["query_sequences"] = n
            result.setdefault("workarounds", []).append(
                f"wrote a flat Q{SYNTHETIC_BASE_QUALITY} for {n:,} rnaSPAdes "
                "contigs, which carry no quality"
            )
        elif method.family == "corrected":
            query = run_isoncorrect(runner, sample, out_dir, threads)
            if method.params.get("assemble") == "isonform":
                # Assemble inside each cluster instead of leaving the corrected
                # reads as they are. The cluster is already allele-separated, so
                # this should not average away the minority allele the way a
                # global assembler does — which is the whole hypothesis.
                query = run_isonform(runner, sample, out_dir, threads, query)
            with open(query) as handle:
                result["counts"]["query_sequences"] = sum(
                    1 for index, _ in enumerate(handle) if index % 4 == 0
                )
            result["counts"]["spanning_reads_available"] = stats["n_spanning_reads"]
        elif method.family == "assembly":
            assembly_dir = out_dir / "rnabloom"
            if assembly_dir.exists():
                shutil.rmtree(assembly_dir)
            runner.run(
                "rnabloom",
                [
                    "rnabloom",
                    # -long is RNA-Bloom2's long-read mode and is not
                    # platform-specific; the same flag covers ONT and PacBio.
                    # What does differ by platform is the minimap2 preset the
                    # polished contigs are then realigned with.
                    "-long", *(str(path) for path in assembly_inputs(sample)),
                    "-o", str(assembly_dir),
                    "-t", str(threads),
                    "-chimera",
                    "-f",
                ],
            )
            # Step 6 of Andy Lee's canonical Nexus subworkflow, which the
            # Exacto docs point at but do not spell out: drop assembled
            # transcripts without enough read support before anything is
            # aligned. It also emits the FASTQ — assembled contigs have no
            # per-base quality, and call-rna-vars panics on a BAM without one.
            raw_assembly = find_assembly(assembly_dir)
            query_fasta = out_dir / f"{prefix}.transcripts.fa"
            query = out_dir / f"{prefix}.transcripts.fastq.gz"
            runner.run(
                "nexus_filter_rnabloom2_transcripts",
                [
                    "nexus_filter_rnabloom2_transcripts",
                    "--assembly4-pol-fasta-file", str(raw_assembly),
                    "--assembly3-map-paf-file", str(assembly_dir / PAF_NAME),
                    *(
                        argument
                        for key, value in method.params.items()
                        for argument in (f"--{key}", value)
                    ),
                    "--output-reads-tsv-file", str(out_dir / f"{prefix}.reads.tsv"),
                    "--output-transcripts-tsv-file",
                    str(out_dir / f"{prefix}.transcripts.tsv"),
                    "--output-fasta-file", str(query_fasta),
                    "--output-fastq-file", str(query),
                ],
            )
            result["counts"]["assembled_transcripts"] = count_fasta_records(raw_assembly)
            result["counts"]["query_sequences"] = count_fasta_records(query_fasta)
        else:
            # Only the variant-spanning reads, and capped harder than the
            # assembly arm gets them. Without an assembler each read is its own
            # transcript, so a read that touches no vaccine variant can never
            # produce one of the mutant proteins under test — it would just cost
            # Exacto time — and nothing collapses the depth before call-rna-vars,
            # which is what exhausts the runner. See READS_ARM_READS_PER_VARIANT.
            # minimap2 carries their real base qualities into the BAM, which
            # Exacto needs.
            query = reads_arm_fastq(sample)
            result["counts"]["query_sequences"] = stats.get(
                "n_reads_arm_reads", stats["n_spanning_reads"]
            )
            result["counts"]["spanning_reads_available"] = stats["n_spanning_reads"]

        aligned_bam = out_dir / f"{prefix}.aligned.bam"
        align(runner, query, sample, arm, aligned_bam)
        result["counts"]["aligned"] = count_alignments(aligned_bam)

        if (
            method.family == "assembly"
            and query_fasta is not None
            and method.params.get("unspliced-filter") != "off"
        ):
            # Only the documented assembly pipeline filters unspliced RNAs;
            # applied to raw or corrected reads it would throw away most of the
            # data, since an individual read need not look spliced.
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
            if method.family == "assembly":
                # Recorded so the funnel does not silently show a stage that was
                # never run as a stage that lost nothing.
                result["counts"]["unspliced_filter"] = "skipped"

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
                # Exacto's defaults, matching Andy's Nexus subworkflow, which
                # passes no extra arguments here. They are permissive — a DNA
                # variant links to any RNA variant within 10 kb of a transcript
                # boundary or 100 kb intergenically, and only 19 of 3,359
                # integrations were exact in one measured arm — but tightening
                # them would mean testing something other than Exacto as shipped.
                # The verdict does not read this table; see evaluate.py.
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
        print(f"    !! {sample.name}/{arm}: {error}")

    result["steps"] = runner.steps
    result["seconds"] = round(sum(step["seconds"] for step in runner.steps), 1)
    (out_dir / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    import argparse
    import os

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 2)
    parser.add_argument(
        "--samples",
        nargs="*",
        help="Sequencing samples to run, e.g. T1-ONT T1-PacBio (default: all).",
    )
    parser.add_argument("--arms", nargs="*", default=list(ARMS))
    args = parser.parse_args()

    variants = load_variants()
    runs = []
    for sample in samples_named(args.samples):
        skipped = [a for a in args.arms if a not in arms_for(sample, args.arms)]
        for name in skipped:
            print(f"== {sample.name} / {name}: not applicable to "
                  f"{sample.read_type} reads, skipping ==")
        for arm in arms_for(sample, args.arms):
            print(f"== {sample.name} / {arm} ==")
            runs.append(run_arm(sample, arm, variants, args.threads))

    failures = [run for run in runs if run["status"] != "ok"]
    print(f"\n{len(runs) - len(failures)}/{len(runs)} arms completed")
    for failure in failures:
        print(f"  FAILED {failure['sample']}/{failure['arm']}: {failure['error']}")


def collect_runs(samples: list[str] | None = None) -> list[dict]:
    """Gather the per-arm run.json files on disk.

    CI runs the samples as separate jobs, so the runs are assembled from what
    each one left behind rather than from a single in-process list. Only the
    canonical ``<sample>/<arm>/`` paths count — anything else under work/ is
    somebody's leftovers, not a result.
    """
    known_samples = {sample.name for sample in SAMPLES}
    runs = []
    for run_path in sorted(EXACTO_DIR.glob("*/*/run.json")):
        arm_dir = run_path.parent
        if arm_dir.parent.name not in known_samples or arm_dir.name not in ARMS:
            continue
        run = json.loads(run_path.read_text())
        if samples and run["sample"] not in samples:
            continue
        runs.append(run)
    return runs


if __name__ == "__main__":
    main()
