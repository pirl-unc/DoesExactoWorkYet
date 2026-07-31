# translate-structs emits nothing for deletion RNA calls

**Severity:** silent data loss
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto translate-structs`

## What happens

Every deletion in the test set is called in the RNA by call-rna-vars and none reaches a proteoform, on any sample or any arm: 22 of 22 translated variants are SNVs, 0 of 6 are deletions.

This is not coverage and not a marginal call. H1-2 carries 2,998 alt reads of 18,173 at T2, and call-rna-vars emits 150 separate RNA variant calls for it on T1-PacBio and 114 on T1-ONT.

Instrumenting the primary-structures table settles where it is lost. Counting every row that names each RNA call id, split by what translate-structs said the amino acid did, on T1-ONT/reads:

    SNVs        17 of 18 called variants are referenced, with mutant rows
    deletions    0 of 4 called variants are referenced at all

So the deletion calls are not referenced and judged unchanged, and not dropped by the consumer -- translate-structs emits no primary-structure row naming them. The same pattern holds on all four samples. The SNV column is the control: the same instrument, the same file, the same run.

## Workaround currently in use

None. Deletions cannot be recovered as mutant proteoforms in this version, which for a neoantigen application removes frameshift and inframe-deletion neoantigens entirely -- 6 of the 37 mutations here, 5 of them frameshifts.

## Suggested fix

A test that translates the same transcript with an SNV call and with a deletion call at the same position and asserts both produce mutant proteoforms would catch this; see docs/reproductions.md. Fixing it takes recovery in this test from 22 of 37 to 28 of 37 -- 28 of the 29 mutations whose allele is present in the RNA.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
