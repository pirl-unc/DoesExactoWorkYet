const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const test = require("node:test");
const { sequenceRows, highlightedPieces, sharedOffset, referenceGroups,
  sequenceComparison, alignmentSlice, reconstructionWindows, windowSlice, rnaSupport, sequenceStatus, sourceRnaState } = require("../web/peptides.js");

const reportPath = path.join(__dirname, "../results/vaccine_peptide_analysis.json");
const report = fs.existsSync(reportPath) ? JSON.parse(fs.readFileSync(reportPath, "utf8")) : null;
const fixture = {
  analysis: { artifacts: [
    { sample: "T1-ONT", arm: "reads" },
    { sample: "T1-PacBio", arm: "corrected" },
    { sample: "T2-ILMN", arm: "spades" },
  ] },
  variants: [
    { variant_id: "SPG11", gene: "SPG11", protein_change: "p.Ile147Val", published_vaccine_peptides: [
      { sequence: "ATILYSCSREALQKLIDDQDVSISLLSLRIL", in_vaccines: ["mRNA"], matches: [] },
      { sequence: "KLIDDQDVSI", in_vaccines: ["mRNA"], matches: [] },
    ] },
    { variant_id: "ZNF436", gene: "ZNF436", protein_change: "p.Ser207Cys", published_vaccine_peptides: [
      { sequence: "KSFGRSCHL", in_vaccines: ["mRNA"], matches: [
        { sample: "T1-ONT", arm: "reads" }, { sample: "T1-PacBio", arm: "corrected" },
      ] },
    ] },
  ],
};

test("the per-sequence report reproduces the recorded containment totals", { skip: !report }, () => {
  const rows = sequenceRows(report);
  const found = rows.filter((row) => row.status === "contained");
  assert.equal(rows.length, report.summary.n_peptide_entries);
  assert.equal(found.length, report.summary.n_peptide_entries_matched);
  assert.equal(new Set(found.map((row) => row.variant.variant_id)).size,
    report.summary.n_variants_any_peptide_matched);
});

test("T2-ONT absence is not displayed as a sequence miss", () => {
  const rows = sequenceRows(fixture, { sample: "T2-ONT" });
  assert.equal(rows.length, 3);
  assert.ok(rows.every((row) => row.status === "not_evaluated"));
  assert.equal(sequenceRows(fixture, { sample: "T2-ONT", status: "not_found" }).length, 0);
});

test("an unavailable sample/method combination is unscored", () => {
  const rows = sequenceRows(fixture, { sample: "T2-ILMN", method: "reads" });
  assert.ok(rows.every((row) => row.status === "not_evaluated"));
});

test("method and sample filters cannot inherit matches from other runs", () => {
  const rows = sequenceRows(fixture, { sample: "T1-PacBio", method: "corrected", status: "contained" });
  assert.ok(rows.length > 0);
  assert.ok(rows.every((row) => row.matches.every((hit) => hit.sample === "T1-PacBio" && hit.arm === "corrected")));
});

test("gene, amino acid, vaccine, and status filters combine", () => {
  const rows = sequenceRows(fixture, { query: "spg11", vaccine: "mRNA", status: "not_found" });
  assert.equal(rows.length, 2);
  assert.ok(rows.every((row) => row.variant.gene === "SPG11"));
  const peptide = rows[0].peptide.sequence;
  assert.ok(sequenceRows(fixture, { query: peptide.toLowerCase() }).some((row) => row.peptide.sequence === peptide));
  assert.equal(sequenceRows(fixture, { query: "not-a-gene" }).length, 0);
});

test("highlight preserves longer reconstructions and combines overlapping occurrences", () => {
  const pieces = highlightedPieces("AAACGKSFGRSCGGG", [4], 9);
  assert.equal(pieces.map((piece) => piece.text).join(""), "AAACGKSFGRSCGGG");
  assert.equal(pieces.filter((piece) => piece.matched).map((piece) => piece.text).join(""), "CGKSFGRSC");
  assert.deepEqual(highlightedPieces("AAAAA", [1, 2], 3).filter((piece) => piece.matched),
    [{ text: "AAAA", matched: true }]);
});

test("every displayed reconstruction contains the highlighted reference peptide", { skip: !report }, () => {
  for (const { peptide, matches } of sequenceRows(report)) {
    for (const match of matches) {
      const protein = report.reconstructions[match.reconstruction_id];
      assert.ok(protein);
      for (const start of match.amino_acid_starts) {
        assert.equal(protein.slice(start - 1, start - 1 + peptide.sequence.length), peptide.sequence);
      }
    }
  }
});

test("reference offsets align nested peptides and overlapping extensions", () => {
  const variant = { published_vaccine_peptides: [
    { sequence: "KSLRPRKVNTPAGSSQKAREERALLPLELQD" },
    { sequence: "ERALLPLEL" },
    { sequence: "GSSQKAREERALLPLELQDDGSDSRKS" },
  ] };
  const groups = referenceGroups(variant.published_vaccine_peptides.map((peptide) => ({ variant, peptide })));
  assert.equal(groups.length, 1);
  const offsets = new Map(groups[0].map((entry) => [entry.row.peptide.sequence, entry.offset]));
  assert.equal(offsets.get("ERALLPLEL"), 20);
  assert.equal(offsets.get("GSSQKAREERALLPLELQDDGSDSRKS"), 12);
  assert.equal(sharedOffset("AAAAAAXAAAAAA", "AAAAAA"), null, "repeated anchors are ambiguous");
  assert.equal(sharedOffset("KELPLYLWQPSTSEIAVIRDWKK", "KKSVIRTLSTIDDVEDRENEKGR"), null);
});

test("trimmed regions retain shared columns, highlights, and full-protein coordinates", () => {
  const slice = alignmentSlice("ABCDEFGHIJKLMN", -5, -2, 6, [[5, 8], [7, 9]]);
  assert.equal(slice.from, 3);
  assert.equal(slice.to, 11);
  assert.equal(slice.pieces.map((piece) => piece.text).join(""), "DEFGHIJK");
  assert.equal(slice.pieces.filter((piece) => piece.matched).map((piece) => piece.text).join(""), "FGHI");
  assert.ok(slice.clippedLeft && slice.clippedRight);
  const short = alignmentSlice("FGHI", 0, -2, 6);
  assert.equal(short.left, "  ");
  assert.equal(short.right, "  ");
  assert.equal(alignmentSlice("FGHI", 20, 0, 10), null);
});

test("identical window sequences combine full proteins and preserve their provenance", () => {
  const proteins = [
    { sequence: "MPEPTIDEKK*", offset: -1, ranges: [[1, 8]] },
    { sequence: "MQQPEPTIDEK*", offset: -3, ranges: [] },
    { sequence: "MPEPTVDEKK*", offset: -1, ranges: [] },
  ];
  const windows = reconstructionWindows({ reconstructions: proteins }, 0, 7, 4);
  assert.equal(windows.length, 2);
  assert.deepEqual(windows[0].members.map((m) => m.protein), proteins.slice(0, 2));
  assert.equal(windows[1].members[0].protein, proteins[2]);
  assert.deepEqual(windows.map((w) => w.number), [4, 5]);
  assert.equal(windows[0].sequence, "PEPTIDE");
  assert.equal(windowSlice(windows[0], 0, 7).pieces.filter((p) => p.matched).map((p) => p.text).join(""), "PEPTIDE");
});

test("visible missing ends, different placement, and stops are not labeled identical", () => {
  const proteins = [
    { sequence: "PEPTIDE", offset: 0, ranges: [] },
    { sequence: "PEPTIDE*", offset: 0, ranges: [] },
    { sequence: "PEPTIDE", offset: 1, ranges: [] },
    { sequence: "PEPTID", offset: 0, ranges: [] },
  ];
  const windows = reconstructionWindows({ reconstructions: proteins }, 0, 8);
  assert.equal(windows.length, 4);
  assert.ok(windows.every((w) => w.members.length === 1));
});

test("window equivalence is computed before screen wrapping and keeps differences in later blocks", () => {
  const proteins = [
    { sequence: "MPEPTIDEKKK", offset: -1, ranges: [] },
    { sequence: "MPEPTIDEKKV", offset: -1, ranges: [] },
  ];
  const windows = reconstructionWindows({ reconstructions: proteins }, 0, 10);
  assert.equal(windows.length, 2);
  assert.deepEqual(windowSlice(windows[0], 0, 7).pieces, windowSlice(windows[1], 0, 7).pieces);
  assert.notDeepEqual(windowSlice(windows[0], 7, 10).pieces, windowSlice(windows[1], 7, 10).pieces);
});

test("comparison windows are 45 aa or wide enough to preserve every reference", () => {
  for (const length of [9, 31, 45, 80]) {
    const peptide = { sequence: "A".repeat(length), matches: [] };
    const [group] = sequenceComparison({ reconstructions: {} }, [{ peptide, matches: [] }]);
    assert.equal(group.end - group.start, Math.max(45, length));
    assert.ok(group.start <= 0 && group.end >= length);
  }
});

test("equivalent window support deduplicates RNA inputs while keeping samples and methods separate", () => {
  const runData = { analysis: { artifacts: [
    { sample: "T1-ONT", arm: "reads" }, { sample: "T1-ONT", arm: "corrected" },
    { sample: "T3-ONT", arm: "reads" },
  ] } };
  const hits = [
    { sample: "T1-ONT", arm: "reads", protein_id: "r1|orf_0-29" },
    { sample: "T1-ONT", arm: "reads", protein_id: "r1|orf_3-32" },
    { sample: "T1-ONT", arm: "corrected", protein_id: "r1|orf_0-29" },
    { sample: "T3-ONT", arm: "reads", protein_id: "r1|orf_0-29" },
    { sample: "T3-ONT", arm: "reads", protein_id: "r2|orf_0-29" },
  ];
  const proteins = hits.map((hit, i) => ({
    sequence: "MPEPTIDE" + "K".repeat(i), offset: -1, ranges: [],
    observations: new Map([[String(i), { ...hit, reconstruction_id: String(i) }]]),
  }));
  const [window] = reconstructionWindows({ reconstructions: proteins }, 0, 7);
  assert.equal(window.members.length, 5);
  assert.equal(window.observations.size, 5);
  assert.deepEqual(rnaSupport(runData, [...window.observations.values()]), [
    { sample: "T1-ONT", evaluated: true, methods: [{ method: "corrected", count: 1 }, { method: "reads", count: 1 }] },
    { sample: "T3-ONT", evaluated: true, methods: [{ method: "reads", count: 2 }] },
  ]);
});

test("all supporting sequences and anchored occurrences survive inline grouping", { skip: !report }, () => {
  for (const variant of report.variants) {
    const rows = sequenceRows({ ...report, variants: [variant] });
    const groups = sequenceComparison(report, rows);
    const expected = new Set(rows.flatMap((row) => row.matches.map((match) => match.reconstruction_id)));
    (variant.reconstructed_candidates || []).forEach((hit) => expected.add(hit.reconstruction_id));
    const displayed = new Set(groups.flatMap((group) => group.reconstructions.flatMap((protein) =>
      [...protein.observations.values()].map((hit) => hit.reconstruction_id))));
    assert.deepEqual(displayed, expected, variant.gene);
    assert.equal(groups.reduce((n, group) => n + group.references.length, 0), rows.length);
    for (const group of groups) {
      const windows = reconstructionWindows(group, group.start, group.end);
      const groupedIds = new Set(windows.flatMap((window) => [...window.observations.values()].map((hit) => hit.reconstruction_id)));
      const groupIds = new Set(group.reconstructions.flatMap((protein) => [...protein.observations.values()].map((hit) => hit.reconstruction_id)));
      assert.deepEqual(groupedIds, groupIds, `${variant.gene}: window grouping preserves every output`);
      for (const { row, offset } of group.references) {
        for (const match of row.matches) {
          for (const start of match.amino_acid_starts) {
            assert.ok(group.reconstructions.some((protein) =>
              protein.sequence === report.reconstructions[match.reconstruction_id] &&
              protein.offset + start - 1 === offset &&
              protein.ranges.some(([a, b]) => a === start - 1 && b === start - 1 + row.peptide.sequence.length)));
          }
        }
      }
    }
  }
});

test("protein disagreement is distinct from no translation and unavailable runs", () => {
  const variant = { reconstructed_candidates: [{ sample: "T1-ONT", arm: "reads" }] };
  assert.equal(sequenceStatus(fixture, variant, [], {}), "sequence_disagreement");
  assert.equal(sequenceStatus(fixture, variant, [], { sample: "T1-PacBio" }), "no_sequence");
  assert.equal(sequenceStatus(fixture, variant, [], { sample: "T2-ONT" }), "not_evaluated");
  assert.equal(sequenceStatus(fixture, variant, [{ sample: "T1-ONT", arm: "reads" }]), "contained");
});

test("source RNA distinguishes mutant support, measured zero, no coverage and missing counts", () => {
  assert.equal(sourceRnaState([{ alt_reads: 4, total_reads: 20 }]), "supported");
  assert.equal(sourceRnaState([{ alt_reads: 0, total_reads: 20 }]), "no_alt_reads");
  assert.equal(sourceRnaState([{ alt_reads: 0, total_reads: 0 }]), "no_coverage");
  assert.equal(sourceRnaState([{ alt_reads: null, total_reads: 20 }]), "not_reported");
  assert.equal(sourceRnaState([]), "not_reported");
});

test("source RNA evidence can find missing proteins without crossing selected samples", () => {
  const data = structuredClone(fixture);
  data.variants[0].source_rna_support = [
    { benchmark_sample: "T1-ONT", alt_reads: 3, total_reads: 20 },
    { benchmark_sample: "T2-ILMN", alt_reads: 0, total_reads: 25 },
  ];
  const rows = sequenceRows(data, { status: "no_sequence", sourceRna: "supported" });
  assert.equal(rows.length, 2);
  assert.ok(rows.every((row) => row.variant.gene === "SPG11"));
  assert.equal(sequenceRows(data, { sample: "T2-ILMN", sourceRna: "supported" }).length, 0);
  assert.equal(sequenceRows(data, { sample: "T2-ILMN", sourceRna: "no_alt_reads" }).length, 2);
});

test("sample and method filters also exclude non-matching reconstructed proteins", { skip: !report }, () => {
  const variant = report.variants.find((v) => v.gene === "SPG11");
  const data = { ...report, variants: [variant] };
  for (const filters of [{ sample: "T3-ONT" }, { method: "corrected" }, { sample: "T2-ONT" }]) {
    const groups = sequenceComparison(data, sequenceRows(data, filters), filters);
    assert.ok(groups.every((group) => group.reconstructions.length === 0));
  }
});

test("identical full proteins combine support but differences outside the crop remain separate", () => {
  const peptide = { peptide_id: "P1", sequence: "PEPTIDE", matches: [
    { reconstruction_id: "a", protein_id: "a|orf_0-29", amino_acid_starts: [2], sample: "T1-ONT", arm: "reads", rna_call_id: "1" },
    { reconstruction_id: "b", protein_id: "b|orf_0-29", amino_acid_starts: [2], sample: "T1-ONT", arm: "reads", rna_call_id: "2" },
    { reconstruction_id: "c", protein_id: "c|orf_0-32", amino_acid_starts: [2], sample: "T1-PacBio", arm: "corrected", rna_call_id: "1" },
  ] };
  const report = { reconstructions: { a: "MPEPTIDEK*", b: "MPEPTIDEK*", c: "MPEPTIDEKK*" } };
  const [group] = sequenceComparison(report, [{ peptide, matches: peptide.matches }]);
  assert.equal(group.reconstructions.length, 2);
  assert.equal(group.reconstructions[0].observations.size, 2);
  assert.deepEqual(rnaSupport(fixture, [...group.reconstructions[0].observations.values()])[0], {
    sample: "T1-ONT", evaluated: true, methods: [{ method: "reads", count: 2 }],
  });
});

test("RNA support counts input RNAs once per sample and method, without counting repeated ORFs or calls", () => {
  const data = { analysis: { samples_unavailable: ["T2-ONT"], artifacts: [
    { sample: "T1-ONT", arm: "reads" }, { sample: "T1-ONT", arm: "corrected" },
    { sample: "T2-ILMN", arm: "spades" },
  ] } };
  const hits = [
    { sample: "T1-ONT", arm: "reads", protein_id: "read-1|orf_0-29", rna_call_id: "1" },
    { sample: "T1-ONT", arm: "reads", protein_id: "read-1|orf_3-29", rna_call_id: "2" },
    { sample: "T1-ONT", arm: "reads", protein_id: "read-2|orf_0-29", rna_call_id: "3" },
    { sample: "T1-ONT", arm: "corrected", protein_id: "read-1|orf_0-29", rna_call_id: "1" },
  ];
  const support = rnaSupport(data, hits);
  assert.deepEqual(support, [
    { sample: "T1-ONT", evaluated: true, methods: [{ method: "corrected", count: 1 }, { method: "reads", count: 2 }] },
    { sample: "T2-ILMN", evaluated: true, methods: [{ method: "spades", count: 0 }] },
    { sample: "T2-ONT", evaluated: false, methods: [] },
  ]);
  assert.deepEqual(rnaSupport(data, hits, { sample: "T2-ILMN", method: "reads" }), [
    { sample: "T2-ILMN", evaluated: false, methods: [] },
  ]);
  assert.deepEqual(rnaSupport(data, hits, { sample: "T1-ONT", method: "reads" })[0].methods,
    [{ method: "reads", count: 2 }]);
});
