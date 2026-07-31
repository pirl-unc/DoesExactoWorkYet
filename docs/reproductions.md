# Reproducing the Exacto findings, with suggested tests

Ten findings from an end-to-end run of Exacto 0.4.6a1. Each has a minimal
reproduction that needs no patient data, and a suggested unit test in the crate
where the fault lives — the point being that every one of these is reachable
from a handful of synthetic records, so none of them need this harness to catch.

Filed by [pirl-unc/DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet).
Live detail, including the failing command and its stderr for each:
[the bug reports page](https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html).

Severity here means: **crash** stops the pipeline; **silent data loss** produces
a plausible wrong answer, which is worse.

---

## 1. `remove-unspliced-rnas` panics on unmapped BAM records

**Severity** crash · `exacto-core/src/common/bam.rs:892`

The reference sequence id of every record is unwrapped without checking, so the
first unmapped record ends the run. minimap2 emits unmapped records by default,
so this is reachable from any alignment that was not pre-filtered.

**Reproduce**

```bash
# Any BAM containing at least one unmapped record.
samtools view -H aligned.bam > tiny.sam
samtools view aligned.bam | head -5 >> tiny.sam        # mapped
samtools view -f 4 aligned.bam | head -1 >> tiny.sam   # one unmapped
samtools view -b tiny.sam | samtools sort -o tiny.bam && samtools index tiny.bam

exacto remove-unspliced-rnas --bam-file tiny.bam --bam-bai-file tiny.bam.bai \
  --fasta-file transcripts.fa --reference-gene-annotation-file genes.gtf \
  --reference-gene-annotation-source gencode \
  --reference-gene-annotation-assembly hg38 \
  --reference-gene-annotation-version v44 \
  --output-bam-file out.bam --output-bam-bai-file out.bam.bai \
  --output-fasta-file out.fa
```

**Expected** unmapped records skipped. **Actual** panic on unwrap.

**Suggested test** — `exacto-core/src/common/bam.rs`

```rust
#[test]
fn unmapped_records_are_skipped_not_unwrapped() {
    // A record with reference_sequence_id = None must not end the run: minimap2
    // emits these by default, so any caller who has not pre-filtered hits it.
    let records = vec![mapped_record("chr1", 1000), unmapped_record()];
    let kept = collect_mapped(&records);
    assert_eq!(kept.len(), 1);
}
```

**Workaround in use** `samtools view -b -F 4` before every Exacto step.

---

## 2. `remove-unspliced-rnas` cannot index the BAM it just wrote

**Severity** crash · `exacto-core/src/common/bam.rs:1315`

Records are written in hash-map iteration order under a header still declaring
`SO:coordinate`, so the indexing step rejects the file the same function
produced. The filtering itself completes — only the index fails.

**Reproduce** run the command in §1 on a BAM with no unmapped records and more
than one reference sequence. Output BAM is written; indexing fails.

**Suggested test** — `exacto-core/src/common/bam.rs`

```rust
#[test]
fn written_records_are_sorted_when_header_claims_coordinate_order() {
    // Records arrive from a HashMap, whose iteration order is unspecified.
    // Writing them unsorted under SO:coordinate makes the file unindexable by
    // the very next call.
    let out = write_filtered(unsorted_records_across_two_contigs());
    assert!(is_coordinate_sorted(&out),
            "header declares SO:coordinate; records must match");
}
```

**Workaround in use** catch the failure, `samtools sort`, reindex, continue.

---

## 3. `integrate-vars` panics on RNA calls that matched no reference transcript

**Severity** crash · `exacto-integrator/src/algorithms/variant_integration.rs:98`

`reference_transcript_ids` is parsed with `field.split(",")`. In Rust
`"".split(",")` yields **one empty string**, not nothing — so a call that matched
no transcript looks like a call against a transcript named `""`, takes the
annotated-transcript branch, and dies on `get_transcript("").unwrap()`.

This is not a rare input: measured across four samples, **969–2,029 RNA calls
per run** matched no reference transcript.

**Reproduce**

```bash
printf 'variant_call_id\ttranscript_model_id\treference_transcript_ids\t...\n1\tTM1\t""\t...\n' > rna.tsv
exacto integrate-vars --annotated-dna-vars-tsv-file dna.tsv \
  --rna-vars-tsv-file rna.tsv --output-tsv-file out.tsv ...
```

**Suggested test** — `exacto-integrator/src/algorithms/variant_integration.rs`

```rust
#[test]
fn empty_reference_transcript_list_is_empty_not_one_empty_string() {
    // "".split(",") yields [""] in Rust. Anything treating that as a transcript
    // id reaches get_transcript("").unwrap().
    assert_eq!(parse_transcript_ids(""), Vec::<String>::new());
    assert_eq!(parse_transcript_ids("\"\""), Vec::<String>::new());
    assert_eq!(parse_transcript_ids("ENST1,ENST2"), vec!["ENST1", "ENST2"]);
}

#[test]
fn rna_call_without_reference_transcript_integrates_as_intergenic() {
    let call = rna_call_with_reference_transcripts("");
    assert!(integrate(&call, &dna_call()).is_ok());
}
```

**Workaround in use** withhold those rows from `integrate-vars` only, copying
lines verbatim — re-serialising through a CSV writer drops Exacto's own `""`
quoting and triggers a *different* panic at `rna_variant_call_set.rs:236`.

---

## 4. `call-rna-vars` silently drops any read whose transcript touches a non-ACGT base

**Severity** silent data loss · `exacto-caller/src/structs/reference_transcript_sequence.rs:176`

No warning, no count. A masked or padded reference makes reads disappear with no
indication that anything was lost — the run reports success on a subset.

**Reproduce** build a reference with an `N` inside a transcript's span, align
reads across it, run `call-rna-vars`. Compare input read count to the number
appearing in the output.

**Suggested test** — `exacto-caller/src/structs/reference_transcript_sequence.rs`

```rust
#[test]
fn reads_touching_ambiguous_bases_are_reported_not_silently_dropped() {
    let reference = reference_with_n_at(1_000);
    let result = build_transcript_sequence(&reference, span(950, 1_050));
    // Whatever the policy, the caller must be able to tell it happened.
    assert!(result.is_err() || result.unwrap().dropped_bases > 0);
}
```

**Workaround in use** grow every masked window until no overlapping transcript
touches an `N`, and fail the build loudly if one still does.

---

## 5. Reference transcript sequences are rebuilt one base at a time, per read

**Severity** performance · `exacto-caller/src/structs/reference_transcript_sequence.rs:150-190`

Each read triggers a fresh per-base walk of its transcript's reference sequence.
Measured cost: **110 s versus 74 s** for the same 1,179 reads, purely from
whether the reference FASTA was bgzipped — because every base pays a BGZF block
decompression.

**Suggested test** — a benchmark rather than a unit test

```rust
#[bench]
fn transcript_sequence_is_not_rebuilt_per_read(b: &mut Bencher) {
    // 1,000 reads over one transcript should cost about one transcript
    // extraction, not one thousand.
    let reference = test_reference();
    b.iter(|| build_sequences_for(&reference, &thousand_reads_one_transcript()));
}
```

**Workaround in use** keep the reference FASTA uncompressed.

---

## 6. `call-rna-vars` panics on every read when the BAM has no base qualities

**Severity** crash · `exacto-caller/src/structs/alignment.rs:229`

`len is 0 but index is 0` on the first read. Reached by two ordinary inputs:
assembled contigs, which carry no `QUAL`, and PacBio Iso-Seq, where
`isoseq groupdedup` writes `QUAL` as `*` on **every** record — measured at 206
of 207 records across three loci.

**Reproduce**

```bash
# A FASTA-derived or Iso-Seq BAM: QUAL is "*" for every record.
samtools view aligned.bam | awk '$11 == "*"' | head   # confirm
exacto call-rna-vars --bam-file aligned.bam ...
```

**Suggested test** — `exacto-caller/src/structs/alignment.rs`

```rust
#[test]
fn records_without_base_qualities_do_not_panic() {
    // QUAL "*" is valid SAM and is what isoseq groupdedup and every
    // FASTA-derived alignment produce.
    let record = aligned_record_without_qualities();
    let result = Alignment::from_record(&record);
    assert!(result.is_ok(), "QUAL '*' is valid SAM, not a malformed record");
}
```

**Workaround in use** write a flat Q30 for records with no `QUAL`, matching what
Nexus does for RNA-Bloom2 contigs.

---

## 7. `longest_orf` picks a frame from the read alone

**Severity** silent data loss · `translate-structs --strategy longest_orf`

A single miscalled indel can make a different frame the longest, and nothing
constrains the choice. Measured: **all 32 frameshifted proteoforms** are
attributed to variants Exacto was told are SNVs of size 1 — a substitution
cannot shift its own frame.

The frameshift rate tracks basecalling accuracy rather than biology: **26%** on
raw ONT reads, **9%** on PacBio consensus reads, **0%** on RNA-Bloom2 contigs. A
real upstream phased indel would not care which sequencer saw it.

Silent rather than a crash: the proteoform is emitted, carries the right variant
call id, and feeds `call-peptide-vars` neoantigen candidates built on wrong
sequence downstream of the mutation.

**Suggested tests** — `exacto-translator`

```rust
#[test]
fn an_snv_cannot_produce_a_frameshifted_proteoform() {
    // The variant class is Exacto's own input, not outside knowledge.
    let structures = translate(&transcript_with_indel_error(), &snv_call());
    for proteoform in structures.iter().filter(|p| p.carries(&snv_call())) {
        assert_eq!(proteoform.frameshift_state, FrameshiftState::InFrame,
                   "a size-1 substitution cannot shift the frame");
    }
}

#[test]
fn frame_follows_the_annotated_cds_when_a_reference_transcript_matched() {
    // One synthetic indel in a homopolymer must not move the frame when the
    // model matched an annotated transcript whose CDS phase is known.
    let clean = translate_with_reference_frame(&transcript(), &reference_cds());
    let noisy = translate_with_reference_frame(&transcript_with_homopolymer_indel(),
                                               &reference_cds());
    assert_eq!(clean.peptide_at(&mutation()), noisy.peptide_at(&mutation()));
}
```

**Suggested fix** a reference-CDS-anchored strategy alongside `longest_orf` and
`all_orfs`, applied where an annotated transcript matched; plus the variant-class
consistency check, which is free.

---

## 8. `call-rna-vars` memory grows ~0.55 MB per read, then spikes ~6.8 GB

**Severity** limitation · `exacto call-rna-vars`

Sampled every 60 s on a 16 GB runner: 14,919 MB available falling to 4,668 MB
over 50 minutes, a 57-minute plateau, then one 60-second window in which
available memory fell to 799 MB and 3,071 MB of swap was consumed.

Peak ≈ `(0.55 MB × reads) + 6.8 GB`, so the ceiling is near **15,000 reads on a
16 GB machine** — below what one long-read sample produces. Confirmed three
ways: 18,818 and 19,858 reads both killed the runner; 7,692 completed in 127
minutes.

The failure mode is worse than the limit: the runner stops responding, so the
job dies with no error of its own and nothing is uploaded.

**Suggested test** — a bounded-memory assertion

```rust
#[test]
fn memory_does_not_scale_without_bound_in_read_count() {
    // 20k reads must not need proportionally more resident memory than 2k.
    let small = peak_rss(|| call_rna_vars(&reads(2_000)));
    let large = peak_rss(|| call_rna_vars(&reads(20_000)));
    assert!(large < small * 4, "10x the reads should not cost ~10x the memory");
}
```

**Suggested fix** stream transcript models to disk or process in bounded
batches. Failing that, documenting the per-read cost turns an unexplained
three-hour death into a sizing decision.

**Workaround in use** cap the reads arm at 600 reads per variant.

---

## 9. Candidate proteoforms are emitted unranked

**Severity** limitation · `translate-structs`

The reads route returns a **median of 5.5 candidate proteins per recovered
mutation, up to 107**, roughly a quarter of them frameshifted, with nothing to
choose between them. A downstream consumer must invent a rule.

**Suggested fix** report per-proteoform **molecule support** — how many
independent reads or UMIs produced each translation — and two flags for frame
provenance: did the transcript model match a reference transcript, and did it
reach that transcript's annotated start codon. Both are information Exacto
already has.

```rust
#[test]
fn proteoforms_report_the_molecules_that_support_them() {
    let structures = translate(&three_identical_reads_and_one_divergent());
    let modal = structures.iter().max_by_key(|p| p.supporting_reads).unwrap();
    assert_eq!(modal.supporting_reads, 3);
}
```

---

## 10. Deletions are called in the RNA and never translated

**Severity** open — attribution not yet established · `translate-structs`

Clean split by variant type: **22 of 22 translated variants are SNVs, 0 of 6
deletions.** Not coverage — H1-2 has 2,998 alt reads of 18,173 and 150 separate
RNA calls from Exacto's own caller, and yields no proteoform on any sample or
any arm.

**Not yet filed as an Exacto bug.** It could equally be this harness's evaluator,
which keys proteoforms on RNA call ids and only keeps rows where
`amino_acid_change == "mutant"`. A diagnostic is instrumented and running: it
counts every primary-structures row naming each RNA call, split by what
`translate-structs` said the amino acid did, which distinguishes *never
referenced the call* from *referenced it and judged the protein unchanged* from
*we dropped it*.

**Suggested test regardless** — `exacto-translator`

```rust
#[test]
fn a_deletion_call_produces_a_mutant_proteoform() {
    // Same coverage, same transcript, only the variant type differs.
    let snv = translate(&transcript(), &snv_call_at(500));
    let del = translate(&transcript(), &deletion_call_at(500, 3));
    assert!(!snv.mutant_proteoforms().is_empty());
    assert!(!del.mutant_proteoforms().is_empty(),
            "an in-frame deletion should translate as readily as a substitution");
}
```

---

## Filing these

Every finding on the [bug reports page](https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html)
has a **Copy as issue** button producing a GitHub-ready body with the version it
was seen on, the source location, the failing command and its stderr tail. The
machine-readable source is
[`results/findings.json`](https://github.com/pirl-unc/DoesExactoWorkYet/blob/main/results/findings.json).

Two caveats on all of it. Every reproduction above is written from what the
harness observed rather than from a minimised case actually re-run in isolation
— the line numbers and messages are real, the synthetic inputs are proposed. And
finding 10 has no attribution yet; it may be ours.
