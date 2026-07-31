# call-rna-vars memory grows ~0.55 MB per read, then spikes ~6.8 GB at the end

**Severity:** limitation
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto call-rna-vars`

## What happens

Handed long reads as transcripts (the reads arm, where nothing assembles them into contigs first), call-rna-vars accumulates resident memory linearly at roughly 0.55 MB per read, then allocates a further ~6.8 GB in a final burst. Sampled every 60 s on a 16 GB GitHub runner for T2-ONT: 14,919 MB available at 00:58 falling to 4,668 MB by 01:48 — 10.2 GB consumed, about 2.5 GB per 10 minutes — then a 57-minute plateau, then a single 60-second window in which available memory fell to 799 MB and 3,071 MB of swap was consumed, a further 6.8 GB. Peak demand ~17 GB. The runner stopped responding and was killed (exit 143) after 1h45m in the step, with no error of its own. The terminal burst after a long plateau suggests a final aggregation or serialisation pass holding a second copy of what was built.

Peak is therefore about (0.55 MB x reads) + 6.8 GB, putting the ceiling near 15,000 reads on a 16 GB machine — below what a single long-read sample produces. Confirmed against three legs of one run rather than modelled from one: T2-ONT (18,818 reads, predicted 16.6 GB) and T3-ONT (19,858, predicted 17.2 GB) both died; T1-PacBio (7,692, predicted 10.7 GB) completed in 127 minutes on the same runner image.

## Workaround currently in use

Cap the reads arm at 600 reads per variant — 42% of the reads, ~8,000 per sample, predicted peak ~10.8 GB, which is the depth the PacBio leg completed at. The assembly arm needs no cap because RNA-Bloom2 collapses the reads into a few hundred contigs before Exacto sees them, and peaks near 2 GB.

## Suggested fix

Stream transcript models to disk rather than holding all of them, or process reads in bounded batches. Failing that, documenting the memory cost per read would let callers size the machine instead of discovering the limit as an unexplained runner death three hours in.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
