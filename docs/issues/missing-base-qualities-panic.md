# call-rna-vars panics on every read when the BAM has no base qualities

**Severity:** crash
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto-caller/src/structs/alignment.rs:229`

## What happens

base_quality_scores[i] is indexed once per base with no length check, so a BAM whose QUAL column is '*' panics on every record ('the len is 0 but the index is 0'). Zero transcript models then get built and the run dies in transcript_model_set.rs:337. This is reachable straight from Exacto's own documented pipeline, which aligns RNA-Bloom2's FASTA output — FASTA carries no qualities. Exacto's bundled test BAMs all have a flat quality string, which is why the tests do not catch it.

## Workaround currently in use

Reads are aligned from FASTQ; assembled contigs get a flat Q40 quality string before alignment.

## Suggested fix

Treat a missing QUAL as a uniform default rather than indexing into an empty slice.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
