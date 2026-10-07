# DoesExactoWorkYet

An automated, end-to-end test of [Exacto](https://github.com/pirl-unc/exacto) on real
data: **can it recover the mutant proteins that went into Sid Sijbrandij's personalised
cancer vaccines, from the long-read RNA-seq of his own tumour?**

The benchmark uses [osteosarc](https://github.com/iskandr/osteosarc) **0.15.0** and
its frozen **2026-09-28** corrected catalogue. The inputs and their checksums live
in `data/osteosarc/`; weekly runs and site rebuilds never refresh the variant panel.
Reads are acquired through osteosarc into verified, cached regional BAMs. Each
sample's prepared FASTQs and reference are shared by all its methods.

**Results: <https://pirl-unc.github.io/DoesExactoWorkYet/>** — a summary page, plus
[data sources &amp; method](https://pirl-unc.github.io/DoesExactoWorkYet/sources.html)
listing every file and parameter, and
[bug reports](https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html) with the failing
commands and their output, ready to file upstream.

---

## The question

The panel contains **51 vaccine-associated variants**, using the union of the
catalogue's vaccine-count assertions, source-variant vaccine flags, and vaccine-overlap
memberships. The overlap table alone contains 38 peptide entries for 37 mutations;
the additional 14 variants are retained with their source assertions. Unknown vaccine
names or protein effects are not inferred. Unresolved alleles remain visible but are
excluded from Exacto inputs and the measured recovery denominator.

The earlier 37-mutation scores remain in history. They are not measurements of this
expanded input set; results carry a catalogue identity, and the site refuses to attach
an older evaluation to the new panel.

A sequence-backed subset is available as
[`results/vaccine_peptide_subset.fasta`](results/vaccine_peptide_subset.fasta)
and [`results/vaccine_peptide_subset.json`](results/vaccine_peptide_subset.json).
It contains **36 variants with 78 recorded vaccine-peptide entries (77 distinct
sequences)**. Inclusion requires both a sequence and an explicit `in_vaccines`
label in `published_vaccine_peptides`. The export retains the recorded long
peptides and mRNA minimal epitopes, their vaccine labels, and the original
mutation annotations. Identical sequences retain separate records where the
source does. NME1 remains included because a vaccine peptide is recorded, despite
its conflicting consequence annotation; ABCF2 has no recorded vaccine peptide.
The subset excludes pVACtools-only predictions and peptides without vaccine
membership. Its denominator describes available reference sequences, not verified
expression in every sample. The full benchmark still evaluates the 51-variant
panel. Rebuild both exports offline with `python -m pipeline.fetch_osteosarc`.

The subset reanalysis is saved in
[`results/vaccine_peptide_analysis.json`](results/vaccine_peptide_analysis.json).
The website's [per-sequence report](https://pirl-unc.github.io/DoesExactoWorkYet/peptides.html)
shows all reference peptides, filters by gene/sequence, vaccine, sample and method,
and highlights each contained peptide inside its reconstructed protein. Missing
sample/method outputs are shown as not evaluated. The report is tied to the exact
scored result so a new run cannot inherit stale containment results.
Its metric is **vaccine sequence contained in reconstruction**: the recorded
vaccine peptide must appear as a contiguous substring of a reconstructed protein,
which can extend on either side. Every counted occurrence is tied to an RNA call
at the target's exact locus and allele, and overlaps its translated codon/junction
or frameshifted tail. The analysis reads the full exported protein FASTAs, since
the main results JSON retains only three candidate examples per target/method.

Re-scoring run 51 (Exacto 0.5.0a1) found **15/36 targets** with at least one
recorded vaccine sequence contained in a reconstruction: **37/78 peptide entries**.
All 15 also matched a non-minimal vaccine peptide. The subset's broader candidate
recovery was 22/36, with the expected residue in 17/19 recovered missense targets.
These are unions over 22 completed methods; T2-ONT input preparation failed and
its methods are unavailable. This reanalysis uses the saved translations and
does not rerun alignment, assembly, or Exacto.

The full benchmark workflow regenerates this report from the current run's
archived proteins and RNA calls before publishing. If those outputs cannot be
verified, the site marks the sequence report unavailable while still publishing
the main benchmark diagnostics.

To reproduce it while the GitHub artifacts remain available:

```bash
mkdir -p work/run-51
git show 01475a2:results/exacto_results.json > work/run-51/exacto_results.json
gh run download 37536966270 --repo pirl-unc/DoesExactoWorkYet \
  --pattern 'exacto-outputs-*' --dir work/run-51/outputs
python -m pipeline.score_vaccine_peptides \
  --results work/run-51/exacto_results.json --outputs-dir work/run-51/outputs \
  --run-url https://github.com/pirl-unc/DoesExactoWorkYet/actions/runs/37536966270 \
  --exacto-version 0.5.0a1
```

The verdict is graded, not binary:

| Outcome | Meaning |
|---|---|
| `no_reads` | the ONT data doesn't cover the locus — not counted against Exacto |
| `no_call` | reads cover it, Exacto called no RNA variant there |
| `rna_only` | Exacto called the variant at that exact locus and allele in the RNA, but translated no mutant protein carrying it |
| `proteoform` | a translated primary structure carries the mutation |
| `peptide` | ...and `call-peptide-vars` emitted novel mutant peptides from it |

Failed or missing runs leave mutations **unscored**, rather than claiming no
reads cover them. If no mutations can be evaluated, the site says **Run failed**
and shows the last evaluated result separately. Completed runs with no covered
mutations say **No covered mutations**. Neither produces a recovery fraction or
replaces the last measurement in the history.

Every rung is keyed on the RNA variant call Exacto made **at the mutation's exact locus
with its exact allele**, and never on `integrate-vars` output. That distinction is not
pedantic: `integrate-vars` links a DNA variant to any RNA variant within 10 kb of a
transcript boundary, or 100 kb intergenically, and at those defaults only 19 of 3,359
integrations in one T1 arm were exact. Scoring off the integration table turned 5 recovered
mutations into 28. `integrate-vars` still runs at Exacto's defaults, as Nexus runs it —
the point is to test Exacto as shipped, not a tuned version of it — but nothing downstream
of the verdict depends on its output.

Two further checks run on top of the ladder:

- **Right residue.** For missense mutations the amino acid Exacto produced is compared
  against the one the portal's HGVS annotation predicts. A change at the right codon but
  the wrong residue is reported, not counted as a win.
- **Predicted epitope match.** For the 15 variants linked to the frozen pVACtools
  predictions, the current scorer asks whether Exacto's translated proteoform
  contains a predicted mutant epitope as a substring. This score does not yet use
  the published vaccine-peptide sequences in the 36-variant subset above, and
  should not be interpreted as recovery of all recorded vaccine peptides.

## What actually runs

```
fetch_osteosarc  →  build_reference  →  extract_reads  →  run_exacto  →  evaluate  →  build_site
```

**`pipeline/fetch_osteosarc.py`** verifies the frozen catalogue and epitope checksums
and uses osteosarc's allele-aware membership join. It builds all 51 variant records
from the corrected alleles, preserving correction IDs, original membership assertions,
ELISPOT annotations and per-assay support. It makes no network requests.

**`pipeline/build_reference.py`** builds a reference that is small but still in hg38
coordinates. Each mutation's GENCODE v44 gene body (±10 kb) is fetched from hg38 over HTTP
byte ranges and written back at its true offset inside otherwise all-N chromosomes — about
8 Mb of real sequence. minimap2 skips N runs when it collects minimizers, so it indexes in
seconds. The windows are then grown until every transcript that overlaps one lies entirely
on real sequence, because Exacto drops any read whose candidate transcript touches an N,
and the build fails loudly if one still escapes. The GTF is subset to the same windows.

**`pipeline/extract_reads.py`** asks osteosarc for verified whole-gene regional BAMs
from the pinned ONT, PacBio and Illumina source files, then samples locally. Coverage is wildly uneven: the
mitochondrial window alone holds ~1.3M reads, 88% of everything in scope, while VPS13B's
variant has 20. So reads land in two files:

- **spanning** — the read's alignment covers a vaccine variant. The only reads that can
  carry a mutation, capped per variant at a depth no caller needs to exceed.
- **context** — anything else in the gene. Interchangeable filler that helps RNA-Bloom2
  extend transcripts, capped per region.

Both are sampled by seeded reservoir; gzip timestamps are fixed, so rebuilding produces identical FASTQ bytes.
`pipeline.prepare_inputs` pins checksums for those files and the reference, and each method
verifies them before running. Changed preparation code or input checksums invalidate the cache.
The uncapped counts are recorded alongside, and shown per variant on the site.

**`pipeline/run_exacto.py`** runs each timepoint through two arms:

- **`assembly`** — the pipeline as Exacto documents it, by way of Andy Lee's canonical
  `PEPTIDE_PREDICTION_EXACTO` subworkflow in [Nexus](https://github.com/pirl-unc/nexus).
  RNA-Bloom2 assembles spanning + context reads, `nexus_filter_rnabloom2_transcripts`
  drops contigs without enough read support (min MAPQ 30, min 3 reads, min 0.5 fraction
  match — it halved a T1 assembly, 5,365 contigs to 2,656), minimap2 realigns them
  (`splice:hq`), `remove-unspliced-rnas` filters, then `call-rna-transcript-vars`.

  That filter also writes the FASTQ, not just the FASTA. Nexus's own comment explains
  why — *"so the downstream BAM consumed by call-rna-vars carries QUAL fields"* — which
  is the same workaround this harness arrived at independently before finding it, and
  independent confirmation that the crash below is real.
- **`reads`** — the same without the assembler; the spanning reads go straight in as
  transcripts (`splice` preset, no unspliced filter). Without an assembler each read *is*
  a transcript, so a read touching no variant cannot produce one of the mutant proteins
  under test. Cheaper, and it separates an Exacto miss from an assembler miss.

Both arms then feed the known vaccine mutations in as the somatic DNA callset and run
`annotate-vars` → `integrate-vars` → `translate-transcripts` → `call-peptide-vars`.
These are the Exacto 0.5 commands; explicit 0.4 installations retain the older
`call-rna-vars`/`translate-structs` path. The 0.5 translator also receives the
original transcript sequences in a separate support TSV. Scoring reads its
per-nucleotide RNA provenance, including deletion events, and requires a novel
peptide to overlap that call's translated residue or frameshifted tail.
Read-support counts come from the input reads or Nexus's contig-to-read table;
they are unavailable for rnaSPAdes/isONform contigs, whose names are not counted
as supporting reads.

The corrected methods require the `spoa` executable in `environment.yml`.
isONform receives corrected clusters and emits `transcriptome.fastq`; worker
stderr is retained when either correction or isoform assembly fails.

The **Exacto compatibility** PR check builds the current release and runs a
small synthetic input through raw reads, a preassembled contig with splice
filtering, corrected reads, and isONform, all the way through peptide scoring.
It exercises the installed CLIs and table contracts without downloading patient
reads. Run it locally in the installed environment with
`DEWY_EXACTO_INTEGRATION=1 python -m pytest tests/test_exacto_integration.py -v`.

Realignment is not optional: the portal's BAMs were produced with
`minimap2 -ax splice --MD`, without the `--cs` tag that Exacto reads variants from.

**`pipeline/evaluate.py`** scores every mutation against every run. It scores per timepoint
into `results/scored/`, then merges into `results/exacto_results.json` — scoring has to
happen next to the run because the primary-structures TSVs it reads are far too big to
move between CI jobs. **`pipeline/build_site.py`** renders the GitHub Pages site and
appends the run to `results/history.json`, so the answer to "does it work yet" has a track
record rather than just a current value.

Anything Exacto did that looked like a bug rather than a result is written up by hand in
[`results/findings.json`](results/findings.json) and shown on the site.

## Data volume

The corrected catalogue, source inventory and epitope tables are vendored with
checksums. Large alignments stay remote. Preparation runs once per sample, and
every method downloads the same verified reference and capped FASTQs:

| Source | Transferred | |
|---|---|---|
| GENCODE v44 GTF + protein translations | 58 MB | fetched during preparation |
| hg38 `.fai` | 160 KB | fetched during preparation |
| hg38 sequence for the vaccine gene bodies | a few MB | HTTP byte ranges, not the 3 GB FASTA |
| Exacto release tarball | 67 MB | |
| ONT BAM index | 12–17 MB | tells htslib which blocks to ask for |
| ONT BAM reads over the gene windows | several GB | acquired once per sample; full 51-locus volume pending measurement |

The three ONT BAMs total 157 GB and are never downloaded whole. Osteosarc acquires
indexed gene windows from the official public S3 URLs and retains a verified local
BAM and source receipt. The earlier 37-locus T2 extraction transferred about 4.6 GB;
the expanded panel has not yet been measured. Whole-gene context can include many
records even when only a capped fraction is retained. Prepared references and FASTQs
are cached by sample, frozen data, environment and preparation code.

Runtime is dominated by `call-rna-vars`, which rebuilds each candidate reference
transcript's sequence one base at a time for every read and caches nothing (see
[`results/findings.json`](results/findings.json)). Budget roughly two to three hours per
timepoint on a four-vCPU runner; the job timeout is set accordingly. That cost is exactly
why the read caps exist.

Long HTTPS reads can fail. An unsuccessful osteosarc acquisition stops preparation;
it cannot become a result reporting absent read coverage. Rebuilding FASTQs from the
same verified BAM produces byte-identical output.

Preparation needs space for regional BAMs, capped FASTQs and the masked reference.
Method jobs receive only the reference, FASTQs and receipts, plus their own Exacto
intermediates. The masked reference is mostly N and deliberately uncompressed.
The `jlumbroso/free-disk-space` step frees disk space while preserving Clang/LLVM;
the full expanded matrix still needs runtime and disk validation.

## Running it yourself

```bash
micromamba env create -f environment.yml     # or conda/mamba
micromamba activate does-exacto-work-yet
bash scripts/install_exacto.sh               # EXACTO_VERSION=latest-release by default

export DEWY_WORK_DIR=$PWD/work               # big intermediates live here
python -m pipeline.prepare_inputs
python -m pipeline.run_exacto --threads "$(nproc)"
python -m pipeline.evaluate
python -m pipeline.build_site
python -m http.server -d site 8000
```

`run_exacto` takes `--samples T1-ONT` and `--arms reads` if you want a quick single pass.

To test an unreleased Exacto, set `EXACTO_VERSION` to any git ref:

```bash
EXACTO_VERSION=dev bash scripts/install_exacto.sh
```

## Automation

| Workflow | Trigger | Does |
|---|---|---|
| `.github/workflows/exacto-test.yml` | weekly cron, manual dispatch, pushes to `pipeline/` | the full run, commits `results/`, publishes the site |
| `.github/workflows/site.yml` | pushes to `web/` or `results/`, manual dispatch | rebuilds and publishes the site only (~2 min) |
| `.github/workflows/ci.yml` | every push and PR | unit tests, site build and reproducible frozen catalogue rebuild |
| `.github/workflows/install-exacto.yml` | PR changes to pipeline/toolchain, manual dispatch | builds Exacto and runs the synthetic compatibility tests |

Five preparation jobs create one verified input artifact per sample. The 29
sample/method jobs each score their own output and upload a compact JSON;
a final job merges, publishes and commits. Both publishing workflows call
`actions/configure-pages` with `enablement: true`, so the first run turns Pages on without
anyone touching repo settings.

The manual dispatch takes an `exacto_version` input, so testing a candidate release is a
one-click job.

## Layout

```
pipeline/          the six steps, each runnable on its own
data/osteosarc/    frozen corrected catalogue, epitopes and source manifest
web/               three pages — summary, data sources & method, bug reports;
                   build_site.py copies these to site/ alongside data.json
scripts/           Exacto installer, records the exact build under test
tests/             the fiddly bits — variant encoding, region sampling and retry,
                   the streaming proteoform reader, epitope matching
results/           committed outputs — variant table, findings, scored run, history
environment.yml    conda environment (samtools, minimap2, RNA-Bloom2, Exacto's stack)
```

The fast tests run without pysam or samtools installed; installed-toolchain tests
are skipped unless explicitly enabled in the compatibility environment.

## Audited against the author's own pipeline

Checked line by line against Andy Lee's `PEPTIDE_PREDICTION_EXACTO` subworkflow in
[Nexus](https://github.com/pirl-unc/nexus) and against Exacto's source:

| | |
|---|---|
| minimap2 RNA args | identical: `-ax splice:hq -uf --cs --eqx -Y -L --secondary=no` |
| RNA-Bloom2 | `-chimera`, matching; `-lrpb` deliberately dropped — it means PacBio, and this is ONT |
| `nexus_filter_rnabloom2_transcripts` | run, at its defaults |
| `call-rna-vars`, `translate-structs`, `integrate-vars` | no extra args, as Nexus |
| `call-peptide-vars --min-k 8 --max-k 11` | equals Exacto's own defaults |
| `--preset ont` | not applicable — only the DNA callers take a preset, and this pipeline does not run them |
| `samtools calmd` | Nexus pipes through it; safe to omit, because Exacto reads the `cs` tag and never reads `MD` |
| GENCODE v44 vs Nexus's placeholder v45 | deliberate: v44 is what the source BAMs were aligned against |
| Gene/transcript levels 1–3 vs the default 1–2 | deliberate: MT-ND5 is level 3 and would otherwise be dropped silently |

## Caveats

This measures one thing well and several things not at all.

- **DNA variants are supplied, not discovered.** Sid's WGS/WES is short-read; Exacto's DNA
  callers want long reads. The portal's curated somatic calls stand in as the DNA callset.
  The question asked is whether Exacto finds them *in the RNA* and translates them.
- **Only the vaccine genes are analysed.** This is sensitivity at known loci, not
  genome-wide precision. Masking the rest of the genome also removes paralogues that would
  otherwise compete for alignments, which makes alignment easier than it would be
  genome-wide.
- **Contextual reads are capped**, so the assembly arm sees less depth than a whole-sample
  run would in the most highly expressed windows. Variant-spanning reads are never capped.
- **`call-peptide-vars` uses the tested genes' own reference proteins** as the wild-type
  background, so "novel" means absent from that gene's reference isoforms rather than from
  the whole human proteome.
- **GENCODE levels 1–3** are allowed, rather than Exacto's default 1–2, because the
  mitochondrial genes are annotated at level 3 and MT-ND5 is one of the vaccine targets.
- **Exacto has no documented single-cell path.** Not one mention of single-cell,
  barcode, UMI or 10x in its docs, README or Python — yet every long-read dataset here is
  single-cell. Nexus has unwired scRNA tooling (`find_scrna_barcode_knee`,
  `convert_scrna_bam2fastq`, `count_scrna_assembly_support`, `filter_scrna_assembly`) that
  sketches an intended path this test does not yet follow.
- **There is no pooled mode.** Only `call-somatic-dna-vars --control-bam-files` takes more
  than one sample, and those are the *normals*. Everything else is one BAM at a time, and
  Nexus runs Exacto per sample. Pooling T1+T2+T3 would mean merging BAMs by hand.
- **Four Exacto crashes are worked around** rather than reported as failures, or the run
  would stop at the first step every time. Each is written up in
  [`results/findings.json`](results/findings.json) and shown on the site, and any
  workaround applied to a run is recorded in that run's JSON. Without them, Exacto 0.4.6a1
  cannot complete its own documented pipeline on this data.
- **A recovered proteoform is not automatically the neoantigen.** `translate-structs
  --strategy longest_orf` has no reference CDS to anchor on, and on short or truncated
  transcript models it regularly picks the wrong frame — three of the first five
  proteoforms recovered carried the wrong residue at the right codon. That is why the
  residue check exists and why the headline reports it separately.

## Credit

The data is Sid Sijbrandij's, published openly at [osteosarc.com](https://osteosarc.com)
alongside the rest of his osteosarcoma research. Exacto is from
[PIRL at UNC](https://github.com/pirl-unc).
