# longest_orf picks a frame from the read alone, and read error moves it

**Severity:** silent data loss
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto translate-structs --strategy longest_orf`

## What happens

Across the four reads-arm runs, 145 proteoforms were translated and 32 are labelled frameshifted. Every one is attributed to a variant declared an SNV of size 1 in the somatic DNA callset. A substitution cannot itself change a reading frame, so in each case something else moved it.

Two things can: a real frame-shifting variant upstream on the same molecule, or a miscalled indel in the read. The first is legitimate and is exactly what long reads are for — call-rna-vars works de novo, so it can and should call an upstream somatic indel this test never supplied, and a downstream missense codon genuinely does sit in a shifted frame when one is phased with it. None of the 37 vaccine variants shares a gene with another, so none can explain these through the supplied callset, but an unsupplied somatic indel could.

The platform gradient is what separates the two at the population level. 26% of proteoforms frameshift on raw ONT reads, 9% on PacBio Iso-Seq consensus reads, and 0% on RNA-Bloom2 contigs, which are a consensus over many reads. A real upstream somatic frameshift is a property of the tumour and would appear at the same rate on every platform; an error rate that tracks basecalling accuracy is basecalling. So most of these are read error, though this evidence is distributional and does not adjudicate any single case.

It is silent rather than a crash: the proteoform is emitted, carries the right variant call id, and looks like a result, so call-peptide-vars derives neoantigen candidates from sequence that is wrong downstream of the mutation.

## Workaround currently in use

None that is oracle-free at the harness level. This test can detect the bad frames only because it holds the portal's HGVS annotation to compare against — a real discovery run, where Exacto is the only tool, has no such check and would carry the wrong peptides forward.

## Suggested fix

Take the frame from the reference CDS of the transcript the model was matched to, and let only called variants move it — then read error cannot shift the frame at all, while a genuine phased upstream indel still does, which is the distinction that matters. Reporting how many independent molecules support each frame would let a caller separate a recurrent upstream indel from a one-read basecalling artefact without any outside knowledge. The variant class is a weaker but free check: an SNV cannot be the cause of its own frameshift, so a shifted frame implies an upstream cause that should be nameable.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
