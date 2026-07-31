# remove-unspliced-rnas cannot index the BAM it just wrote

**Severity:** crash
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto-core/src/common/bam.rs:1315`

## What happens

The kept records are written in hash-set iteration order while the header still carries SO:coordinate, so the bam::fs::index call on the next line fails with 'invalid reference sequence ID'. The filtering work itself completes — the output BAM is on disk and correct — but the command exits non-zero and the pipeline stops. Observed with 55 out-of-order records in a 1,755-record output.

## Workaround currently in use

The failure is caught, the output BAM is sorted with samtools and reindexed, and the arm continues.

## Suggested fix

Sort the kept records by (reference sequence id, alignment start) before writing, or drop the SO:coordinate claim and index with a sort first.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
