# Reference transcript sequences are rebuilt one base at a time, per read

**Severity:** performance
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto-caller/src/structs/reference_transcript_sequence.rs:150-190`

## What happens

from_reference_transcript issues a separate fasta_reader.query() per base per exon, and is called again for every candidate transcript of every read with no cache. A CPU sample of call-rna-vars found 1,599 of 1,618 samples inside that one function. Cost therefore scales as reads x isoforms x transcript length, which is what makes read depth the dominant runtime factor here. A bgzipped reference makes it worse still, since each query then also pays for a BGZF block decompression: 110 s versus 74 s for the same 1,179 reads.

## Workaround currently in use

Read counts are capped per variant and per region, and the masked reference is written uncompressed.

## Suggested fix

Fetch each exon's sequence in one query, and memoise per transcript id for the lifetime of the run.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
