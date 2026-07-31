# Ready-to-file issue bodies

Generated from `results/findings.json` by `scripts/build_issues.py`.
Nothing here has been filed — these are drafts.

| file | severity | title |
|---|---|---|
| [`unmapped-records-panic.md`](unmapped-records-panic.md) | crash | remove-unspliced-rnas panics on unmapped BAM records |
| [`unsorted-filter-output.md`](unsorted-filter-output.md) | crash | remove-unspliced-rnas cannot index the BAM it just wrote |
| [`empty-transcript-list-panic.md`](empty-transcript-list-panic.md) | crash | integrate-vars panics on RNA calls that matched no reference transcript |
| [`non-acgt-reference-panic.md`](non-acgt-reference-panic.md) | silent data loss | call-rna-vars drops any read whose candidate transcript touches a non-ACGT base |
| [`per-base-reference-queries.md`](per-base-reference-queries.md) | performance | Reference transcript sequences are rebuilt one base at a time, per read |
| [`missing-base-qualities-panic.md`](missing-base-qualities-panic.md) | crash | call-rna-vars panics on every read when the BAM has no base qualities |
| [`longest-orf-frame-on-partial-transcripts.md`](longest-orf-frame-on-partial-transcripts.md) | limitation | longest_orf picks the wrong frame on short or partial transcripts |
| [`unmapped-records-affect-canonical-pipeline-too.md`](unmapped-records-affect-canonical-pipeline-too.md) | crash | The unmapped-record crash is reachable from Nexus's own subworkflow |
| [`call-rna-vars-memory-scales-with-read-count.md`](call-rna-vars-memory-scales-with-read-count.md) | limitation | call-rna-vars memory grows ~0.55 MB per read, then spikes ~6.8 GB at the end |
| [`longest-orf-ignores-the-variant-class-it-was-given.md`](longest-orf-ignores-the-variant-class-it-was-given.md) | silent data loss | longest_orf picks a frame from the read alone, and read error moves it |
