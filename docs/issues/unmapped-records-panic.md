# remove-unspliced-rnas panics on unmapped BAM records

**Severity:** crash
**Observed in:** Exacto 0.4.6a1
**Where:** `exacto-core/src/common/bam.rs:892`

## What happens

record.reference_sequence_id().unwrap().unwrap() is called for every record that is not spliced. minimap2 emits unmapped reads by default, and samtools sort keeps them, so any BAM built straight from the documented alignment command takes the tool down. Observed with 123 unmapped records out of 5,560.

## Workaround currently in use

The pipeline runs `samtools view -b -F 4` between alignment and sorting.

## Suggested fix

Skip records whose reference sequence id is None instead of unwrapping.

## Reproduction and a suggested test

See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/docs/reproductions.md) for a minimal reproduction and a suggested unit test in the crate where this lives.

---

Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), an automated end-to-end test of Exacto on open data from [osteosarc.com](https://osteosarc.com). Full detail, including the failing command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html
