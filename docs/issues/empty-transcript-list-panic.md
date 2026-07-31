# integrate-vars panics on RNA calls that matched no reference transcript

**Severity:** crash
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto-integrator/src/algorithms/variant_integration.rs:98`

## What happens

reference_transcript_ids is parsed with field.split(","), and in Rust "".split(",") yields one empty string rather than nothing. A call against a novel transcript model therefore looks like a call against a transcript named "", fails the is_empty() check that guards the intergenic path, takes the annotated-transcript branch, and dies on get_transcript("").unwrap(). call-rna-vars produces these routinely — 4,039 of 30,697 calls in one T1 arm.

## Workaround currently in use

integrate-vars is handed a copy of the callset with those rows removed, and the number withheld is recorded in the run. translate-structs still receives every call.

## Suggested fix

Filter empty strings out after the split, so a call with no reference transcript takes the intergenic branch it was written for.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
