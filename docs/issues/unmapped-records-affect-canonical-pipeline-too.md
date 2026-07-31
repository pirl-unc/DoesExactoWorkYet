# The unmapped-record crash is reachable from Nexus's own subworkflow

**Severity:** crash
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto-core/src/common/bam.rs:892`

## What happens

Noted while auditing against Andy Lee's PEPTIDE_PREDICTION_EXACTO subworkflow in Nexus. It pipes minimap2 through 'samtools view -bS | samtools calmd | samtools sort' with no -F 4, so unmapped records reach remove-unspliced-rnas exactly as they did here. This is the same defect as unmapped-records-panic, recorded separately because it means the canonical pipeline is affected, not just this harness's variation on it.

## Workaround currently in use

This harness drops unmapped records with samtools view -F 4 before sorting.

## Suggested fix

Same as unmapped-records-panic: skip records whose reference sequence id is None.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
