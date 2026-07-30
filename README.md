# DoesExactoWorkYet

An automated, end-to-end test of [Exacto](https://github.com/pirl-unc/exacto) on real
data: **can it recover the mutant proteins that went into Sid Sijbrandij's personalised
cancer vaccines, from the long-read RNA-seq of his own tumour?**

Everything comes from the open data portal at [osteosarc.com](https://osteosarc.com) —
the vaccine neoantigen list, the somatic variant calls, and the ONT single-cell long-read
RNA-seq from three recurrence biopsies (T1, T2, T3). Nothing is vendored; each run fetches
the current data.

**Results: <https://pirl-unc.github.io/DoesExactoWorkYet/>** — a summary page, plus
[data sources &amp; method](https://pirl-unc.github.io/DoesExactoWorkYet/sources.html)
listing every file and parameter, and
[bug reports](https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html) with the failing
commands and their output, ready to file upstream.

---

## The question

Sid's five personalised vaccines (an mRNA vaccine, three JLF peptide versions, and a
CeGaT vaccine) between them encode 38 neoantigen peptides drawn from 37 somatic mutations.
Each of those is a mutant protein that somebody committed to manufacturing — and 14 of
them have a positive ELISPOT, so we know his T cells saw them.

That makes them an unusually good benchmark. Exacto's job is to go from long reads to
mutant proteoforms. So: point it at the tumour RNA and see how many of those 37 known
mutant proteins it puts back together.

The verdict is graded, not binary:

| Outcome | Meaning |
|---|---|
| `no_reads` | the ONT data doesn't cover the locus — not counted against Exacto |
| `no_call` | reads cover it, Exacto called no RNA variant there |
| `rna_only` | Exacto called the variant at that exact locus and allele in the RNA, but translated no mutant protein carrying it |
| `proteoform` | a translated primary structure carries the mutation |
| `peptide` | ...and `call-peptide-vars` emitted novel mutant peptides from it |

Every rung is keyed on the RNA variant call Exacto made **at the mutation's exact locus
with its exact allele**, and never on `integrate-vars` output. That distinction is not
pedantic: `integrate-vars` links a DNA variant to any RNA variant within 10 kb of a
transcript boundary, or 100 kb intergenically, and at those defaults only 19 of 3,359
integrations in one T1 arm were exact. Scoring off the integration table turned 5 recovered
mutations into 28. The pipeline also runs `integrate-vars` with those tolerances closed
down, so its output means something.

Two further checks run on top of the ladder:

- **Right residue.** For missense mutations the amino acid Exacto produced is compared
  against the one the portal's HGVS annotation predicts. A change at the right codon but
  the wrong residue is reported, not counted as a win.
- **Right peptide.** The portal publishes the curated pVACtools run the vaccine designs
  were picked from, whose `MT Epitope Seq` column is the closest thing available to the
  peptides that were actually manufactured. For the 10 mutations it covers, the test asks
  whether Exacto's translated proteoform *literally contains* that epitope as a substring.
  That is the strictest available form of "the mutant proteoform matches what was in the
  vaccine".

## What actually runs

```
fetch_osteosarc  →  build_reference  →  extract_reads  →  run_exacto  →  evaluate  →  build_site
```

**`pipeline/fetch_osteosarc.py`** pulls `vaccine_overlap.json` and
`variant_vafs_long.tsv` from osteosarc.com and merges them on locus, giving each vaccine
mutation its ref/alt alleles, consequence, ELISPOT status, VAF trend, and — usefully — the
portal's own genotyping of the same ONT BAMs, which is the yardstick for what is even
recoverable.

**`pipeline/build_reference.py`** builds a reference that is small but still in hg38
coordinates. Each mutation's GENCODE v44 gene body (±10 kb) is fetched from hg38 over HTTP
byte ranges and written back at its true offset inside otherwise all-N chromosomes — about
8 Mb of real sequence. minimap2 skips N runs when it collects minimizers, so it indexes in
seconds. The windows are then grown until every transcript that overlaps one lies entirely
on real sequence, because Exacto drops any read whose candidate transcript touches an N,
and the build fails loudly if one still escapes. The GTF is subset to the same windows.

**`pipeline/extract_reads.py`** range-reads the three dedup ONT BAMs — 37, 67 and 53 GB,
which stay on Backblaze — for reads over those windows. Coverage is wildly uneven: the
mitochondrial window alone holds ~1.3M reads, 88% of everything in scope, while VPS13B's
variant has 20. So reads land in two files:

- **spanning** — the read's alignment covers a vaccine variant. The only reads that can
  carry a mutation, capped per variant at a depth no caller needs to exceed.
- **context** — anything else in the gene. Interchangeable filler that helps RNA-Bloom2
  extend transcripts, capped per region.

Both are sampled by seeded reservoir, so re-running an extraction reproduces it exactly.
The uncapped counts are recorded alongside, and shown per variant on the site.

**`pipeline/run_exacto.py`** runs each timepoint through two arms:

- **`assembly`** — the pipeline as Exacto documents it. RNA-Bloom2 assembles spanning +
  context reads into full-length transcripts, minimap2 realigns them (`splice:hq`),
  `remove-unspliced-rnas` filters, then `call-rna-vars`.
- **`reads`** — the same without the assembler; the spanning reads go straight in as
  transcripts (`splice` preset, no unspliced filter). Without an assembler each read *is*
  a transcript, so a read touching no variant cannot produce one of the mutant proteins
  under test. Cheaper, and it separates an Exacto miss from an assembler miss.

Both arms then feed the known vaccine mutations in as the somatic DNA callset and run
`annotate-vars` → `integrate-vars` → `translate-structs` → `call-peptide-vars`.

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

Nothing is vendored — every run pulls the data fresh — and nothing large is downloaded
whole. Per CI job (one timepoint):

| Source | Transferred | |
|---|---|---|
| GENCODE v44 GTF + protein translations | 58 MB | `actions/cache`d across runs |
| hg38 `.fai` | 160 KB | cached |
| hg38 sequence for 37 gene bodies | ~6 MB | HTTP byte ranges, not the 3 GB FASTA |
| Exacto release tarball | 67 MB | |
| ONT BAM index | 12–17 MB | tells htslib which blocks to ask for |
| ONT BAM reads over the vaccine loci | **~4.6 GB** | measured on T2; 7% of that 67 GB file |

The three ONT BAMs total 157 GB and are never downloaded whole. Backblaze serves them
with `accept-ranges: bytes`, so htslib fetches only the BGZF blocks covering the gene
windows. That still comes to ~4.6 GB per timepoint (measured: 396 s wall for T2) because
htslib has to *read* every record in those windows even though the caps mean only a
fraction is kept — the mitochondrial window alone holds over a million reads. The
reference is built the same way: `samtools faidx <url> chr9:70248978-70364873` pulls one
gene, not one genome, for about 6 MB in total.

Since each timepoint is its own CI job, that is ~4.7 GB per job rather than 14 GB in one.

Runtime is dominated by `call-rna-vars`, which rebuilds each candidate reference
transcript's sequence one base at a time for every read and caches nothing (see
[`results/findings.json`](results/findings.json)). Budget roughly two to three hours per
timepoint on a four-vCPU runner; the job timeout is set accordingly. That cost is exactly
why the read caps exist.

Long HTTPS reads against the bucket occasionally die with an HTTP/2 framing error — seen
twice in testing — so both the reference build and the read extraction retry with backoff,
and the extraction restarts a region from scratch rather than resuming half-populated
state. Re-running an extraction produces byte-identical output.

Working set on disk stays well inside a GitHub runner: ~150 MB of extracted FASTQ, a 1.7 GB
reference (mostly N, deliberately uncompressed — see the runtime note above) and the Exacto
intermediates. The `jlumbroso/free-disk-space` step at the top of the workflow clears the
runner's preinstalled toolchains for headroom.

## Running it yourself

```bash
micromamba env create -f environment.yml     # or conda/mamba
micromamba activate does-exacto-work-yet
bash scripts/install_exacto.sh               # EXACTO_VERSION=latest-release by default

export DEWY_WORK_DIR=$PWD/work               # big intermediates live here
python -m pipeline.fetch_osteosarc
python -m pipeline.build_reference
python -m pipeline.extract_reads
python -m pipeline.run_exacto --threads "$(nproc)"
python -m pipeline.evaluate
python -m pipeline.build_site
python -m http.server -d site 8000
```

`run_exacto` takes `--timepoints T1` and `--arms reads` if you want a quick single pass.

To test an unreleased Exacto, set `EXACTO_VERSION` to any git ref:

```bash
EXACTO_VERSION=dev bash scripts/install_exacto.sh
```

## Automation

| Workflow | Trigger | Does |
|---|---|---|
| `.github/workflows/exacto-test.yml` | weekly cron, manual dispatch, pushes to `pipeline/` | the full run, commits `results/`, publishes the site |
| `.github/workflows/site.yml` | pushes to `web/` or `results/`, manual dispatch | rebuilds and publishes the site only (~2 min) |
| `.github/workflows/ci.yml` | every push and PR | unit tests plus a smoke test that osteosarc.com still serves the tables |

The timepoints run as a three-way matrix; one job doing all three would eat most of
GitHub's six-hour ceiling. Each job scores its own timepoint and uploads a compact JSON,
and a final job merges, publishes and commits. Both publishing workflows call
`actions/configure-pages` with `enablement: true`, so the first run turns Pages on without
anyone touching repo settings.

The manual dispatch takes an `exacto_version` input, so testing a candidate release is a
one-click job.

## Layout

```
pipeline/          the six steps, each runnable on its own
web/               three pages — summary, data sources & method, bug reports;
                   build_site.py copies these to site/ alongside data.json
scripts/           Exacto installer, records the exact build under test
tests/             the fiddly bits — variant encoding, region sampling and retry,
                   the streaming proteoform reader, epitope matching
results/           committed outputs — variant table, findings, scored run, history
environment.yml    conda environment (samtools, minimap2, RNA-Bloom2, Exacto's stack)
```

The tests run without pysam or samtools installed, which is what lets the site and CI
workflows stay lightweight.

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
