# longest_orf picks the wrong frame on short or partial transcripts

**Severity:** limitation
**Observed in:** Exacto 0.4.6a1
**Where:** `translate-structs --strategy longest_orf`

## What happens

Three of the five mutant proteoforms recovered at T1 carried the wrong amino acid at the right codon — DYNC1H1 gave H where the annotation predicts I, VPS13B N for F, MT-ND5 V for T. All three were short ORFs (92-129 aa) on novel or truncated transcript models, one terminating in a stop immediately after the mutation. With no reference CDS to anchor on, the longest open reading frame in a partial transcript is often not the real one, and the peptide that comes out is not the neoantigen.

## Workaround currently in use

Not worked around. Every recovered proteoform's residue is checked against the portal's HGVS annotation and a mismatch is reported rather than counted, so the headline is not inflated by these.

## Suggested fix

Where a transcript model matches a reference transcript, prefer that transcript's CDS frame over the longest ORF, or expose the reference frame as a strategy.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
