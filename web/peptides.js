"use strict";

function sequenceRows(report, filters = {}) {
  const query = (filters.query || "").trim().toLowerCase();
  return report.variants.flatMap((variant) => variant.published_vaccine_peptides.map((peptide) => {
    const matches = peptide.matches.filter((hit) =>
      (!filters.sample || hit.sample === filters.sample) &&
      (!filters.method || hit.arm === filters.method));
    return { variant, peptide, matches, status: sequenceStatus(report, variant, matches, filters) };
  })).filter(({ variant, peptide, status }) =>
    (!filters.status || filters.status === "all" || filters.status === status ||
      (filters.status === "not_found" && ["sequence_disagreement", "no_sequence"].includes(status))) &&
    (!filters.vaccine || peptide.in_vaccines.includes(filters.vaccine)) &&
    (!filters.sourceRna || sourceRnaState((variant.source_rna_support || []).filter((entry) =>
      !filters.sample || entry.benchmark_sample === filters.sample)) === filters.sourceRna) &&
    (!query || [variant.gene, variant.protein_change, variant.variant_id,
      peptide.sequence, ...peptide.in_vaccines].join(" ").toLowerCase().includes(query)));
}

function selectedRun(run, filters) {
  return (!filters.sample || run.sample === filters.sample) && (!filters.method || run.arm === filters.method);
}

function sequenceStatus(report, variant, matches, filters = {}) {
  if (!report.analysis.artifacts.some((run) => selectedRun(run, filters))) return "not_evaluated";
  if (matches.some((hit) => selectedRun(hit, filters))) return "contained";
  const hasProtein = (variant.reconstructed_candidates || []).some((run) => selectedRun(run, filters)) ||
    (variant.exacto_runs || []).some((run) => selectedRun(run, filters) && run.n_proteoforms > 0);
  return hasProtein ? "sequence_disagreement" : "no_sequence";
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

// Place references by their longest unambiguous shared stretch. This is a
// display offset, not an inferred gapped alignment or an additional match.
function sharedOffset(reference, sequence) {
  const offsets = [];
  let best = 5;
  for (let offset = 1 - sequence.length; offset < reference.length; offset++) {
    let run = 0, longest = 0;
    for (let i = Math.max(0, offset); i < Math.min(reference.length, offset + sequence.length); i++) {
      run = reference[i] === sequence[i - offset] ? run + 1 : 0;
      longest = Math.max(longest, run);
    }
    if (longest > best) { best = longest; offsets.length = 0; }
    if (longest === best && best >= 6) offsets.push(offset);
  }
  return offsets.length === 1 ? { offset: offsets[0], shared: best } : null;
}

function referenceGroups(rows) {
  const remaining = [...rows].sort((a, b) => b.peptide.sequence.length - a.peptide.sequence.length);
  const groups = [];
  while (remaining.length) {
    const entries = [{ row: remaining.shift(), offset: 0 }];
    while (remaining.length) {
      let best = null;
      remaining.forEach((row, index) => entries.forEach((entry) => {
        const placement = sharedOffset(entry.row.peptide.sequence, row.peptide.sequence);
        if (placement && (!best || placement.shared > best.shared)) {
          best = { index, offset: entry.offset + placement.offset, shared: placement.shared };
        }
      }));
      if (!best) break;
      entries.push({ row: remaining.splice(best.index, 1)[0], offset: best.offset });
    }
    const origin = Math.min(...entries.map((entry) => entry.offset));
    entries.forEach((entry) => { entry.offset -= origin; });
    groups.push(entries);
  }
  return groups;
}

function sequenceComparison(report, rows, filters = {}, windowSize = 45) {
  const groups = referenceGroups(rows).map((references) => {
    const proteins = new Map();
    for (const { row, offset } of references) {
      for (const match of row.matches) {
        for (const start of match.amino_acid_starts) {
          const proteinOffset = offset - (start - 1);
          // A reconstruction matching several vaccine peptides appears once at
          // each supported placement. Repeated occurrences retain separate rows.
          const sequence = report.reconstructions[match.reconstruction_id];
          const key = JSON.stringify([sequence, proteinOffset]);
          if (!proteins.has(key)) proteins.set(key, {
            sequence, offset: proteinOffset, observations: new Map(),
            ranges: [], peptideIds: new Set(),
          });
          const protein = proteins.get(key);
          protein.observations.set(JSON.stringify([match.reconstruction_id, match.rna_call_id]), match);
          protein.peptideIds.add(row.peptide.peptide_id);
          protein.ranges.push([start - 1, start - 1 + row.peptide.sequence.length]);
        }
      }
    }
    const reconstructions = [...proteins.values()].sort((a, b) =>
      b.observations.size - a.observations.size || a.sequence.localeCompare(b.sequence) || a.offset - b.offset);
    return { references, reconstructions, start: 0,
      end: Math.max(...references.map(({ row, offset }) => offset + row.peptide.sequence.length)) };
  });
  const shown = new Set(groups.flatMap((group) => group.reconstructions.flatMap((protein) =>
    [...protein.observations.values()].map((hit) => hit.reconstruction_id))));
  const unaligned = { references: [], reconstructions: [],
    start: -Math.floor(windowSize / 2), end: Math.ceil(windowSize / 2), unaligned: true };
  const placements = new Map();
  for (const candidate of rows[0].variant?.reconstructed_candidates || []) {
    if (!selectedRun(candidate, filters) || shown.has(candidate.reconstruction_id)) continue;
    const sequence = report.reconstructions[candidate.reconstruction_id];
    if (!placements.has(sequence)) {
      let best = null;
      for (const group of groups) for (const entry of group.references) {
        const position = sharedOffset(entry.row.peptide.sequence, sequence);
        if (position && (!best || position.shared > best.shared))
          best = { group, offset: entry.offset + position.offset, shared: position.shared };
      }
      placements.set(sequence, best);
    }
    const placement = placements.get(sequence);
    const group = placement?.group || unaligned;
    const offset = placement ? placement.offset : 1 - candidate.variant_amino_acid_start;
    let protein = group.reconstructions.find((p) => p.sequence === sequence && p.offset === offset);
    if (!protein) {
      protein = { sequence, offset, observations: new Map(), ranges: [], peptideIds: new Set() };
      group.reconstructions.push(protein);
    }
    protein.observations.set(JSON.stringify([candidate.reconstruction_id, candidate.rna_call_id]), candidate);
  }
  for (const group of groups) {
    // The comparison window is independent of screen width. Expand rather than
    // truncate when the aligned vaccine references span more than 45 residues.
    const padding = Math.max(0, windowSize - (group.end - group.start));
    group.start -= Math.floor(padding / 2);
    group.end += Math.ceil(padding / 2);
  }
  if (unaligned.reconstructions.length) groups.push(unaligned);
  return groups;
}

function alignmentSlice(sequence, offset, start, end, ranges = []) {
  const from = Math.max(0, start - offset);
  const to = Math.min(sequence.length, end - offset);
  if (to <= from) return null;
  const pieces = [];
  let text = "", previous = null;
  for (let i = from; i < to; i++) {
    const matched = ranges.some(([a, b]) => i >= a && i < b);
    if (previous !== null && previous !== matched) { pieces.push({ text, matched: previous }); text = ""; }
    text += sequence[i];
    previous = matched;
  }
  pieces.push({ text, matched: previous });
  return { from, to, left: " ".repeat(Math.max(0, offset - start)),
    right: " ".repeat(Math.max(0, end - offset - sequence.length)), pieces,
    clippedLeft: from > 0, clippedRight: to < sequence.length };
}

function reconstructionWindows(group, start, end, firstNumber = 1) {
  const identical = new Map();
  for (const protein of group.reconstructions) {
    const slice = alignmentSlice(protein.sequence, protein.offset, start, end, protein.ranges);
    if (!slice) continue;
    // Compare the complete displayed fragment and its columns, including stops
    // and missing ends. Highlighting is evidence, not a sequence difference.
    const sequence = slice.pieces.map((piece) => piece.text).join("");
    const key = JSON.stringify([slice.left, sequence, slice.right]);
    if (!identical.has(key)) identical.set(key, {
      sequence, offset: start + slice.left.length, number: firstNumber + identical.size,
      members: [], observations: new Map(), ranges: [],
    });
    const window = identical.get(key);
    window.members.push({ protein, slice });
    for (const [key, hit] of protein.observations || []) window.observations.set(key, hit);
    for (const [a, b] of protein.ranges) {
      const from = Math.max(a, slice.from), to = Math.min(b, slice.to);
      if (to > from) window.ranges.push([from - slice.from, to - slice.from]);
    }
  }
  return [...identical.values()];
}

function windowSlice(window, start, end) {
  const slice = alignmentSlice(window.sequence, window.offset, start, end, window.ranges);
  if (slice) {
    slice.clippedLeft ||= window.members.some((member) => member.slice.clippedLeft);
    slice.clippedRight ||= window.members.some((member) => member.slice.clippedRight);
  }
  return slice;
}

const statusText = { contained: "Contained", sequence_disagreement: "Sequence disagreement",
  no_sequence: "No reconstructed sequence", not_evaluated: "Not evaluated" };

function rnaSupport(report, matches, filters = {}) {
  const samples = [...new Set(report.analysis.artifacts.map((run) => run.sample)
    .concat(report.analysis.samples_unavailable || []))].sort()
    .filter((sample) => !filters.sample || sample === filters.sample);
  return samples.map((sample) => {
    const methods = [...new Set(report.analysis.artifacts.filter((run) => run.sample === sample &&
      (!filters.method || run.arm === filters.method)).map((run) => run.arm))].sort();
    return { sample, evaluated: methods.length > 0, methods: methods.map((method) => {
      // Different ORFs, calls, and peptide matches from the same input RNA
      // must not inflate support. Methods reuse inputs, so never sum them.
      const transcripts = new Set(matches.filter((match) => match.sample === sample && match.arm === method)
        .map((match) => match.protein_id.replace(/\|orf_\d+-\d+$/, "")));
      return { method, count: transcripts.size };
    }) };
  });
}

function supportGrid(report, matches, filters, variant = null) {
  const grid = node("dl", "sequence-support");
  for (const support of rnaSupport(report, matches, filters)) {
    const cell = node("div", "sequence-support-sample");
    cell.append(node("dt", "", support.sample));
    if (variant) cell.append(node("dd", "sequence-support-verdict", statusText[sequenceStatus(report, variant, matches, { ...filters, sample: support.sample })]));
    if (!support.evaluated && !variant) cell.append(node("dd", "sequence-support-unavailable", "Not evaluated"));
    const displayed = support.methods.filter((entry) => entry.count || entry.method === "reads" || support.methods.length === 1);
    for (const { method, count } of displayed) {
      const label = method === "reads" ? `raw read${count === 1 ? "" : "s"}`
        : method === "corrected" ? `corrected RNA${count === 1 ? "" : "s"}` : `${method} RNA${count === 1 ? "" : "s"}`;
      const value = node("dd", count ? "sequence-support-positive" : "sequence-support-zero");
      value.append(node("strong", "", String(count)), document.createTextNode(` ${label}`));
      cell.append(value);
    }
    const otherZeros = support.methods.length - displayed.length;
    if (otherZeros) cell.append(node("dd", "sequence-support-zero", `0 in ${otherZeros} other tested method${otherZeros === 1 ? "" : "s"}`));
    grid.append(cell);
  }
  return grid;
}

function sourceRnaState(entries) {
  if (entries.some((entry) => entry.alt_reads > 0)) return "supported";
  if (!entries.length || entries.some((entry) => entry.alt_reads === null || entry.alt_reads === undefined ||
    entry.total_reads === null || entry.total_reads === undefined)) return "not_reported";
  return entries.some((entry) => entry.total_reads > 0) ? "no_alt_reads" : "no_coverage";
}

function sourceRnaEvidence(report, variant, filters) {
  const section = node("div", "sequence-source-evidence");
  section.append(node("h4", "", "Sid dataset: mutant RNA at this locus"), node("p", "sequence-methods",
    "Published mutant / total RNA reads. This evidence is independent of Exacto and does not establish that a read spans the entire vaccine peptide. Libraries are shown separately, without pooling."));
  const states = { supported: "Mutant RNA observed", no_alt_reads: "No mutant reads observed",
    no_coverage: "No coverage", not_reported: "Not reported" };
  const entries = variant.source_rna_support || [];
  const samples = rnaSupport(report, [], filters);
  const grid = node("div", "sequence-source-samples");
  for (const { sample, evaluated } of samples) {
    const rows = entries.filter((entry) => entry.benchmark_sample === sample);
    const state = sourceRnaState(rows);
    const item = node("div", `sequence-source-sample ${state}`);
    item.append(node("strong", "", sample), node("span", "", states[state]));
    for (const row of rows) item.append(node("span", "sequence-source-count", `${row.alt_reads ?? "?"} / ${row.total_reads ?? "?"} reads`));
    const selected = { ...filters, sample };
    const runs = (variant.exacto_runs || []).filter((run) => selectedRun(run, selected));
    const hasProtein = (variant.reconstructed_candidates || []).some((run) => selectedRun(run, selected)) || runs.some((run) => run.n_proteoforms > 0);
    const hasRna = runs.some((run) => run.n_variant_rnas > 0);
    item.append(node("span", "sequence-source-exacto", !evaluated ? "Exacto: not evaluated" : hasProtein ? "Exacto: protein produced"
      : hasRna ? "Exacto: RNA called, no protein" : "Exacto: no variant RNA call or protein"));
    grid.append(item);
  }
  section.append(grid);
  const other = entries.filter((entry) => !entry.benchmark_sample &&
    (!filters.sample || entry.timepoint === filters.sample.split("-")[0]));
  if (other.length) {
    section.append(node("p", "sequence-methods", "Other RNA libraries in the Sid dataset (separate from the benchmark inputs):"));
    const table = node("div", "sequence-other-rna");
    for (const entry of other) {
      const row = node("div", "sequence-other-rna-row");
      row.append(node("span", "", entry.sample_label), node("span", "sequence-source-count", `${entry.alt_reads ?? "?"} / ${entry.total_reads ?? "?"}`),
        node("span", "", states[sourceRnaState([entry])]));
      table.append(row);
    }
    section.append(table);
  }
  return section;
}

function sequenceStrip(slice) {
  const strip = node("code", "sequence-strip");
  strip.append(node("span", "sequence-trim", slice.clippedLeft ? "…" : " "), document.createTextNode(slice.left));
  for (const piece of slice.pieces) strip.append(node(piece.matched ? "mark" : "span", "", piece.text));
  strip.append(document.createTextNode(slice.right), node("span", "sequence-trim", slice.clippedRight ? "…" : " "));
  return strip;
}

function referenceRow(entry, slice, number) {
  const { row } = entry;
  const item = node("div", "alignment-row alignment-reference");
  item.id = `${row.peptide.peptide_id}-part-${number}`;
  const label = node("div", "alignment-label");
  label.append(node("strong", "", `P${row.variant.published_vaccine_peptides.indexOf(row.peptide) + 1} · ${row.peptide.in_vaccines.join(", ")}`),
    node("span", "sequence-meta", `${row.peptide.is_mrna_minimal_epitope ? "Minimal epitope" : "Vaccine peptide"} · ${row.peptide.sequence.length} aa`),
    node("span", `sequence-status ${row.status}`, statusText[row.status]));
  item.append(label, sequenceStrip(slice));
  return item;
}

function reconstructionRow(report, window, slice, filters) {
  const { members, observations, number } = window;
  const fullSequences = new Set(members.map(({ protein }) => protein.sequence));
  const lengths = [...new Set([...fullSequences].map((seq) => seq.replace(/\*$/, "").length))].sort((a, b) => a - b);
  const outputs = new Set([...observations.values()].map((hit) => hit.reconstruction_id)).size;
  const item = node("div", "alignment-row alignment-reconstruction");
  const label = node("div", "alignment-label");
  label.append(node("strong", "", `Window sequence W${number}`), node("span", "sequence-meta", `${outputs} protein output${outputs === 1 ? "" : "s"}`));
  label.append(node("span", `sequence-status ${window.ranges.length ? "contained" : "sequence_disagreement"}`,
    window.ranges.length ? "Contains vaccine sequence" : "No contained vaccine sequence"));
  const content = node("div", "alignment-content");
  content.append(sequenceStrip(slice));
  const lengthText = lengths.length <= 6 ? lengths.join(", ") : `${lengths[0]}–${lengths[lengths.length - 1]} (${lengths.length} lengths)`;
  content.append(node("p", "sequence-coordinates",
    `${fullSequences.size} full protein sequence${fullSequences.size === 1 ? "" : "s"} · Full length${lengths.length === 1 ? "" : "s"}: ${lengthText} aa`));
  if (fullSequences.size > 1) content.append(node("p", "sequence-window-note",
    "Identical within this window; differences outside it are grouped together."));
  content.append(supportGrid(report, [...observations.values()], filters));
  item.append(label, content);
  return item;
}

function targetComparison(report, rows, filters, columns) {
  const fragment = document.createDocumentFragment();
  const groups = sequenceComparison(report, rows, filters);
  let nextWindow = 1;
  for (const group of groups) {
    group.windows = reconstructionWindows(group, group.start, group.end, nextWindow);
    nextWindow += group.windows.length;
  }
  const reconstructionIds = new Set(groups.flatMap((group) => group.reconstructions.flatMap((p) => [...p.observations.values()].map((hit) => hit.reconstruction_id))));
  const sequences = new Set(groups.flatMap((group) => group.reconstructions.map((p) => p.sequence)));
  fragment.append(node("p", "sequence-methods", `${rows.length} recorded peptide entr${rows.length === 1 ? "y" : "ies"} · ${reconstructionIds.size} target-linked protein outputs · ${sequences.size} distinct full sequence${sequences.size === 1 ? "" : "s"}, grouped into ${nextWindow - 1} window sequence${nextWindow === 2 ? "" : "s"}.`));
  const support = node("div", "sequence-peptide-support");
  support.append(node("h4", "", "Exacto sequence support for each vaccine peptide"));
  for (const row of rows) {
    const entry = node("div", "sequence-peptide-support-row");
    entry.append(node("strong", "", `P${row.variant.published_vaccine_peptides.indexOf(row.peptide) + 1} · ${row.peptide.in_vaccines.join(", ")} · ${row.peptide.sequence.length} aa`), supportGrid(report, row.matches, filters, row.variant));
    support.append(entry);
  }
  // Keep the reference rows immediately above the reconstructed sequences.
  groups.forEach((group, groupIndex) => {
    if (group.unaligned) fragment.append(node("p", "sequence-alignment-note", "Unaligned candidates: no shared vaccine-sequence anchor. Each 45-aa window is centered on the translated mutation; these rows are not aligned to the vaccine peptides."));
    else if (groups.filter((g) => !g.unaligned).length > 1) fragment.append(node("p", "sequence-alignment-note",
      `Sequence group ${groupIndex + 1}: no unambiguous shared stretch with the other reference group${groups.length > 2 ? "s" : ""}; shown separately.`));
    if (!group.unaligned) fragment.append(node("p", "sequence-alignment-note",
      `${group.end - group.start}-aa comparison window${group.end - group.start > 45 ? " (expanded to include all aligned vaccine peptides)" : " centered on the aligned vaccine peptides"}.`));
    const nParts = Math.ceil((group.end - group.start) / columns);
    for (let start = group.start, part = 1; start < group.end; start += columns, part++) {
      const end = Math.min(start + columns, group.end);
      const panel = node("div", "sequence-alignment");
      if (nParts > 1) panel.append(node("p", "sequence-alignment-note", `Region ${part} of ${nParts}${part < nParts ? " · continues below" : ""}`));
      for (const entry of group.references) {
        const slice = alignmentSlice(entry.row.peptide.sequence, entry.offset, start, end);
        if (slice) panel.append(referenceRow(entry, slice, part));
      }
      for (const window of group.windows) {
        const slice = windowSlice(window, start, end);
        if (slice) panel.append(reconstructionRow(report, window, slice, filters));
      }
      fragment.append(panel);
    }
  });
  fragment.append(support);
  if (!reconstructionIds.size) {
    const explanations = {
      peptide: "Translated candidates were produced, but none contains these vaccine sequences.",
      proteoform: "Translated candidates were produced, but none contains these vaccine sequences.",
      rna_only: "The target RNA variant was called, but no mutant protein carrying it was translated.",
      no_call: "No RNA call for the target allele was made in the available outputs.",
      no_reads: "No reads covered this target in the available outputs.",
    };
    const reason = rows.every((row) => row.status === "not_evaluated") ? "No completed output for this sample and method selection."
      : filters.sample || filters.method ? "No supporting reconstructions in this selection."
      : explanations[rows[0].variant.candidate_outcome] || "No supporting reconstructions in the available outputs.";
    fragment.append(node("p", "sequence-methods", reason));
  }
  return fragment;
}

function renderSequenceRows(report) {
  const filters = {
    query: queryNode("#sequence-search").value,
    status: queryNode("#sequence-status").value,
    sample: queryNode("#sequence-sample").value,
    method: queryNode("#sequence-method").value,
    vaccine: queryNode("#sequence-vaccine").value,
    sourceRna: queryNode("#sequence-source-rna").value,
  };
  const rows = sequenceRows(report, filters);
  const groups = new Map();
  rows.forEach((row) => {
    if (!groups.has(row.variant.variant_id)) groups.set(row.variant.variant_id, []);
    groups.get(row.variant.variant_id).push(row);
  });
  queryNode("#sequence-count").textContent = `${rows.length} peptide entries across ${groups.size} target${groups.size === 1 ? "" : "s"} shown. ` +
    `${rows.filter((row) => row.status === "contained").length} contained; ` +
    `${rows.filter((row) => row.status === "sequence_disagreement").length} sequence disagreement; ` +
    `${rows.filter((row) => row.status === "no_sequence").length} no reconstructed sequence; ` +
    `${rows.filter((row) => row.status === "not_evaluated").length} not evaluated in this selection.`;
  const container = queryNode("#sequence-targets");
  container.replaceChildren();
  const mobile = window.matchMedia("(max-width: 600px)").matches;
  const columns = Math.max(24, Math.min(72, Math.floor((container.clientWidth - (mobile ? 0 : 196)) / (mobile ? 7.3 : 7.9)) - 2));
  for (const entries of groups.values()) {
    const variant = entries[0].variant;
    const target = node("article", "sequence-target");
    const title = node("h3", "sequence-target-title");
    title.append(node("span", "", variant.gene), node("span", "mono", variant.protein_change || "Protein change not annotated"));
    target.append(title);
    target.append(sourceRnaEvidence(report, variant, filters));
    target.append(targetComparison(report, entries, filters, columns));
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
  if (report.source?.snapshot?.name) queryNode("#sequence-provenance").append(node("span", "",
    ` Sid RNA counts use the frozen osteosarc snapshot from ${report.source.snapshot.name}.`));
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
  for (const selector of ["#sequence-search", "#sequence-status", "#sequence-source-rna", ...Object.keys(options)]) {
    queryNode(selector).addEventListener("input", () => renderSequenceRows(report));
  }
  let resizeFrame;
  window.addEventListener("resize", () => {
    cancelAnimationFrame(resizeFrame);
    resizeFrame = requestAnimationFrame(() => renderSequenceRows(report));
  });
  renderSequenceRows(report);
}

async function main() {
  const response = await fetch("vaccine_peptide_analysis.json", { cache: "no-cache" });
  if (!response.ok) throw new Error(`HTTP ${response.status}`);
  populateSequenceReport(await response.json());
}

if (typeof module !== "undefined") module.exports = { sequenceRows, highlightedPieces, sharedOffset, referenceGroups, sequenceComparison, alignmentSlice, reconstructionWindows, windowSlice, rnaSupport, sequenceStatus, sourceRnaState };
if (typeof document !== "undefined") main().catch((error) => {
  queryNode("#sequence-summary").textContent = `The sequence report could not be loaded (${error.message}). Reload the page to try again.`;
});
