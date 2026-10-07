"use strict";

function sequenceRows(report, filters = {}) {
  const available = report.analysis.artifacts.some((run) =>
    (!filters.sample || run.sample === filters.sample) &&
    (!filters.method || run.arm === filters.method));
  const query = (filters.query || "").trim().toLowerCase();
  return report.variants.flatMap((variant) => variant.published_vaccine_peptides.map((peptide) => {
    const matches = peptide.matches.filter((hit) =>
      (!filters.sample || hit.sample === filters.sample) &&
      (!filters.method || hit.arm === filters.method));
    return { variant, peptide, matches,
      status: !available ? "not_evaluated" : matches.length ? "contained" : "not_found" };
  })).filter(({ variant, peptide, status }) =>
    (!filters.status || filters.status === "all" || filters.status === status) &&
    (!filters.vaccine || peptide.in_vaccines.includes(filters.vaccine)) &&
    (!query || [variant.gene, variant.protein_change, variant.variant_id,
      peptide.sequence, ...peptide.in_vaccines].join(" ").toLowerCase().includes(query)));
}

function highlightedPieces(sequence, starts, length) {
  const ranges = starts.map((start) => [start - 1, start - 1 + length]).sort((a, b) => a[0] - b[0]);
  const merged = [];
  for (const range of ranges) {
    if (merged.length && range[0] <= merged[merged.length - 1][1]) {
      merged[merged.length - 1][1] = Math.max(merged[merged.length - 1][1], range[1]);
    } else merged.push([...range]);
  }
  const pieces = [];
  let cursor = 0;
  for (const [start, end] of merged) {
    pieces.push({ text: sequence.slice(cursor, start), matched: false });
    pieces.push({ text: sequence.slice(start, end), matched: true });
    cursor = end;
  }
  pieces.push({ text: sequence.slice(cursor), matched: false });
  return pieces;
}

const node = (tag, className, text) => {
  const item = document.createElement(tag);
  if (className) item.className = className;
  if (text !== undefined) item.textContent = text;
  return item;
};
const queryNode = (selector) => document.querySelector(selector);

function sequenceEvidence(report, row) {
  const details = node("details", "sequence-evidence");
  details.append(node("summary", "", `Inspect ${row.matches.length} supporting reconstruction${row.matches.length === 1 ? "" : "s"}`));
  const label = node("label", "", "Choose a reconstruction");
  const select = node("select");
  select.setAttribute("aria-label", `Reconstruction for ${row.peptide.peptide_id}`);
  row.matches.forEach((match, index) => {
    const option = node("option", "", `${match.sample} / ${match.arm} — ${index + 1}`);
    option.value = String(index);
    select.append(option);
  });
  label.append(select);
  const evidence = node("div", "sequence-reconstruction");
  const render = () => {
    evidence.replaceChildren();
    const match = row.matches[Number(select.value) || 0];
    const protein = report.reconstructions[match.reconstruction_id];
    const start = match.amino_acid_starts[0] - 1;
    const end = start + row.peptide.sequence.length;
    evidence.append(node("p", "section-note", `Vaccine sequence at amino acids ${start + 1}–${end} of a ${protein.replace(/\*$/, "").length}-residue reconstruction.`));
    const context = node("code", "sequence-context");
    context.append(
      node("span", "sequence-flank", (start > 18 ? "…" : "") + protein.slice(Math.max(0, start - 18), start)),
      node("mark", "", protein.slice(start, end)),
      node("span", "sequence-flank", protein.slice(end, end + 18) + (protein.length > end + 18 ? "…" : "")),
    );
    evidence.append(context);
    const full = node("details", "sequence-full");
    full.append(node("summary", "", "Full reconstructed sequence"));
    const sequence = node("code", "sequence-context");
    for (const piece of highlightedPieces(protein, match.amino_acid_starts, row.peptide.sequence.length)) {
      sequence.append(node(piece.matched ? "mark" : "span", "", piece.text));
    }
    full.append(sequence);
    evidence.append(full, node("p", "sequence-identifiers", `Protein: ${match.protein_id}\nRNA call: ${match.rna_call_id}`));
  };
  select.addEventListener("change", render);
  details.append(label, evidence);
  details.addEventListener("toggle", () => { if (details.open && !evidence.childElementCount) render(); });
  return details;
}

function sequenceRow(report, row, filters) {
  const item = node("div", "sequence-row");
  item.id = row.peptide.peptide_id;
  const top = node("div", "sequence-row-top");
  const statusText = { contained: "Contained", not_found: "Not found", not_evaluated: "Not evaluated" };
  top.append(node("span", `sequence-status ${row.status}`, statusText[row.status]),
    node("span", "sequence-meta", `${row.peptide.is_mrna_minimal_epitope ? "mRNA minimal epitope" : "Vaccine peptide"} · ${row.peptide.sequence.length} amino acids · ${row.peptide.in_vaccines.join(", ")}`));
  item.append(top, node("code", "vaccine-sequence", row.peptide.sequence));
  if (row.status === "contained") {
    const groups = [...new Set(row.matches.map((hit) => `${hit.sample} / ${hit.arm}`))];
    item.append(node("p", "sequence-methods", groups.join(" · ")), sequenceEvidence(report, row));
  } else {
    const explanations = {
      peptide: "Translated candidates were produced, but none contains this vaccine sequence.",
      proteoform: "Translated candidates were produced, but none contains this vaccine sequence.",
      rna_only: "The target RNA variant was called, but no mutant protein carrying it was translated.",
      no_call: "No RNA call for the target allele was made in the available outputs.",
      no_reads: "No reads covered this target in the available outputs.",
    };
    const reason = row.status === "not_evaluated" ? "No completed output for this sample and method selection."
      : filters.sample || filters.method ? "Not found in the selected reconstructions."
      : explanations[row.variant.candidate_outcome] || "Not found in the available reconstructions.";
    item.append(node("p", "sequence-methods", reason));
  }
  return item;
}

function renderSequenceRows(report) {
  const filters = {
    query: queryNode("#sequence-search").value,
    status: queryNode("#sequence-status").value,
    sample: queryNode("#sequence-sample").value,
    method: queryNode("#sequence-method").value,
    vaccine: queryNode("#sequence-vaccine").value,
  };
  const rows = sequenceRows(report, filters);
  const groups = new Map();
  rows.forEach((row) => {
    if (!groups.has(row.variant.variant_id)) groups.set(row.variant.variant_id, []);
    groups.get(row.variant.variant_id).push(row);
  });
  queryNode("#sequence-count").textContent = `${rows.length} peptide entries across ${groups.size} target${groups.size === 1 ? "" : "s"} shown. ` +
    `${rows.filter((row) => row.status === "contained").length} contained; ` +
    `${rows.filter((row) => row.status === "not_found").length} not found; ` +
    `${rows.filter((row) => row.status === "not_evaluated").length} not evaluated in this selection.`;
  const container = queryNode("#sequence-targets");
  container.replaceChildren();
  for (const entries of groups.values()) {
    const variant = entries[0].variant;
    const target = node("article", "sequence-target");
    const title = node("h3", "sequence-target-title");
    title.append(node("span", "", variant.gene), node("span", "mono", variant.protein_change || "Protein change not annotated"));
    target.append(title);
    entries.forEach((row) => target.append(sequenceRow(report, row, filters)));
    container.append(target);
  }
  if (!rows.length) container.append(node("p", "empty", "No sequences match these filters. Try another gene, sequence, or recovery status."));
}

function populateSequenceReport(report) {
  if (report.status === "unavailable") {
    queryNode("#sequence-summary").textContent = report.message;
    return;
  }
  const { summary, analysis } = report;
  queryNode("#sequence-summary").textContent = `${summary.n_variants_any_peptide_matched}/${summary.n_variants} targets have a vaccine sequence contained in a reconstruction. ` +
    `${summary.n_peptide_entries_matched}/${summary.n_peptide_entries} recorded peptide entries recovered across all available samples and methods.`;
  const link = node("a", "", `Exacto ${analysis.exacto_version} · benchmark run`);
  link.href = analysis.run_url;
  queryNode("#sequence-provenance").append(link, node("span", "", ` · ${analysis.n_completed_methods} completed methods. Saved outputs were re-scored; a target counts if any candidate contains a recorded vaccine peptide.`));
  const unavailable = analysis.samples_unavailable || [];
  if (unavailable.length) {
    const notice = queryNode("#sequence-availability");
    notice.hidden = false;
    notice.textContent = `Not evaluated: ${unavailable.join(", ")}. No completed Exacto outputs are available for these samples; they contribute neither matches nor sequence-level misses.`;
  }
  const options = {
    "#sequence-sample": [...analysis.samples_available, ...unavailable].sort(),
    "#sequence-method": Object.keys(report.by_method).sort(),
    "#sequence-vaccine": [...new Set(report.variants.flatMap((v) => v.published_vaccine_peptides.flatMap((p) => p.in_vaccines)))].sort(),
  };
  for (const [selector, values] of Object.entries(options)) {
    const select = queryNode(selector);
    for (const value of values) {
      const option = node("option", "", value + (unavailable.includes(value) ? " (not evaluated)" : ""));
      option.value = value;
      select.append(option);
    }
  }
  queryNode("#sequence-controls").hidden = false;
  for (const selector of ["#sequence-search", "#sequence-status", ...Object.keys(options)]) {
    queryNode(selector).addEventListener("input", () => renderSequenceRows(report));
  }
  renderSequenceRows(report);
}

async function main() {
  const response = await fetch("vaccine_peptide_analysis.json");
  if (!response.ok) throw new Error(`HTTP ${response.status}`);
  populateSequenceReport(await response.json());
}

if (typeof module !== "undefined") module.exports = { sequenceRows, highlightedPieces };
if (typeof document !== "undefined") main().catch((error) => {
  queryNode("#sequence-summary").textContent = `The sequence report could not be loaded (${error.message}). Reload the page to try again.`;
});
