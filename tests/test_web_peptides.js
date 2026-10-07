const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const test = require("node:test");
const { sequenceRows, highlightedPieces } = require("../web/peptides.js");

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
