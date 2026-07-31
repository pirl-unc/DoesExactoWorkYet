# call-rna-vars drops any read whose candidate transcript touches a non-ACGT base

**Severity:** silent data loss
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto-caller/src/structs/reference_transcript_sequence.rs:176`

## What happens

Rebuilding a reference transcript's sequence calls Nucleotide::from_str(...).unwrap() on each base, which panics on N. The panic is caught per read ('Panic while processing read name ...') and the read is skipped, so the run completes and simply reports fewer variants — the failure mode is a quiet loss of sensitivity rather than an error. hg38 has N runs in real assembly gaps, so this is reachable on an unmodified reference too.

## Workaround currently in use

The masked reference is grown until every transcript in the subset annotation lies entirely on real sequence, and the build fails if one does not.

## Suggested fix

Map unknown bases to N rather than unwrapping, and surface skipped reads as a count instead of a stderr line per read.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
