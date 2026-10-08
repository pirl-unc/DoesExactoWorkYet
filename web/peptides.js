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
      members: [], observations: new Map(), ranges: [], peptideIds: new Set(),
    });
    const window = identical.get(key);
    window.members.push({ protein, slice });
    for (const id of protein.peptideIds || []) window.peptideIds.add(id);
    for (const [key, hit] of protein.observations || []) window.observations.set(key, hit);
    for (const [a, b] of protein.ranges) {
      const from = Math.max(a, slice.from), to = Math.min(b, slice.to);
      if (to > from) window.ranges.push([from - slice.from, to - slice.from]);
    }
  }
  return [...identical.values()];
}

// Compare residues only where both windows cover the same displayed column.
// Missing ends are coverage differences, never inferred substitutions or gaps.
function windowDifferences(windows, start, end) {
  const comparisons = new Map(windows.map((window) => [window, {
    reference: windows[0].number, positions: new Map(), substitutions: [], missing: 0, additional: 0,
  }]));
  for (let column = start; column < end; column++) {
    const residues = windows.map((window) => window.sequence[column - window.offset] || null);
    const position = column - start + 1;
    windows.forEach((window, i) => {
      const comparison = comparisons.get(window), reference = residues[0], residue = residues[i];
      if (reference !== null && residue !== null && reference !== residue) {
        comparison.substitutions.push({ position, reference, residue });
        comparison.positions.set(column, `Window position ${position}: W${comparison.reference} ${reference} → W${window.number} ${residue}`);
      }
      if (reference !== null && residue === null) comparison.missing++;
      if (reference === null && residue !== null) comparison.additional++;
    });
  }
  return comparisons;
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

const methodCodes = { reads: "r", corrected: "c", assembly: "a", "assembly-ms1": "a1",
  "assembly-permissive": "ap", "assembly-unspliced": "au", isonform: "i", spades: "s" };
const statusSymbols = { contained: "✓", sequence_disagreement: "≠", no_sequence: "∅", not_evaluated: "—" };

function supportCells(report, matches, filters, variant = null) {
  const cells = [];
  for (const support of rnaSupport(report, matches, filters)) {
    const cell = node("td", "alignment-reads");
    const status = variant ? sequenceStatus(report, variant, matches, { ...filters, sample: support.sample }) : null;
    const description = `${support.sample}: ` + (!support.evaluated ? "Not evaluated" : support.methods.map(({ method, count }) => `${method}: ${count}`).join("; ")) +
      (status ? `; ${statusText[status]}` : "");
    cell.title = description;
    cell.setAttribute("aria-label", description);
    if (!support.evaluated) cell.append(node("span", "support-empty", "—"));
    else if (!support.methods.some(({ count }) => count)) {
      cell.append(node("span", "support-empty", "0"));
      if (status && status !== "contained") cell.append(node("span", `support-state ${status}`, statusSymbols[status]));
    } else {
      const displayed = support.methods.filter(({ count, method }) => count || method === "reads")
        .sort((a, b) => (a.method === "reads" ? -1 : b.method === "reads" ? 1 : a.method.localeCompare(b.method)));
      for (const { method, count } of displayed) {
        const value = node("span", count ? "support-count positive" : "support-count");
        value.append(node("strong", "", String(count)), node("span", "support-method", methodCodes[method] || method));
        cell.append(value);
      }
    }
    cells.push(cell);
  }
  return cells;
}

function sourceRnaState(entries) {
  if (entries.some((entry) => entry.alt_reads > 0)) return "supported";
  if (!entries.length || entries.some((entry) => entry.alt_reads === null || entry.alt_reads === undefined ||
    entry.total_reads === null || entry.total_reads === undefined)) return "not_reported";
  return entries.some((entry) => entry.total_reads > 0) ? "no_alt_reads" : "no_coverage";
}

function sourceRnaCount(entry) {
  const state = sourceRnaState(entry ? [entry] : []);
  const labels = { supported: "Mutant RNA observed", no_alt_reads: "No mutant reads observed",
    no_coverage: "No coverage", not_reported: "Not reported" };
  const value = entry && (entry.alt_reads != null || entry.total_reads != null)
    ? `${entry.alt_reads ?? "—"} / ${entry.total_reads ?? "—"}` : "—";
  const count = node("span", `sequence-source-count ${state}`, value);
  count.title = `${labels[state]} · ${value} mutant / total reads`;
  count.setAttribute("aria-label", count.title);
  return count;
}

function sourceRnaEvidence(report, variant, filters) {
  const section = node("div", "sequence-source-evidence");
  section.append(node("h4", "", "RNA at the locus"));
  const entries = variant.source_rna_support || [];
  const samples = rnaSupport(report, [], filters);
  const scroll = node("div", "sequence-table-scroll");
  const grid = node("table", "sequence-locus-table");
  const head = node("thead"), header = node("tr");
  header.append(node("th", "", "Sample"));
  samples.forEach(({ sample }) => header.append(node("th", "", sample)));
  head.append(header);
  const body = node("tbody"), counts = node("tr"), outcomes = node("tr");
  counts.append(node("th", "", "Sid mutant / total"));
  outcomes.append(node("th", "", "Exacto output"));
  for (const { sample, evaluated } of samples) {
    const rows = entries.filter((entry) => entry.benchmark_sample === sample);
    const item = node("td", "sequence-source-sample");
    for (const row of rows.length ? rows : [null]) item.append(sourceRnaCount(row));
    const selected = { ...filters, sample };
    const runs = (variant.exacto_runs || []).filter((run) => selectedRun(run, selected));
    const hasProtein = (variant.reconstructed_candidates || []).some((run) => selectedRun(run, selected)) || runs.some((run) => run.n_proteoforms > 0);
    const hasRna = runs.some((run) => run.n_variant_rnas > 0);
    const outcome = node("td", "sequence-source-exacto", !evaluated ? "—" : hasProtein ? "Protein" : hasRna ? "RNA only" : "No call");
    outcome.title = !evaluated ? "Not evaluated" : hasProtein ? "Target-linked protein produced" : hasRna ? "Variant RNA called, no translated protein" : "No variant RNA call or protein";
    counts.append(item);
    outcomes.append(outcome);
  }
  body.append(counts, outcomes);
  grid.append(head, body);
  scroll.append(grid);
  section.append(scroll);
  const other = entries.filter((entry) => !entry.benchmark_sample &&
    (!filters.sample || entry.timepoint === filters.sample.split("-")[0]));
  if (other.length) {
    const details = node("details", "sequence-details");
    details.append(node("summary", "", `Other Sid libraries (${other.length})`));
    const table = node("div", "sequence-other-rna");
    for (const entry of other) {
      const row = node("div", "sequence-other-rna-row");
      row.append(node("span", "", entry.sample_label), sourceRnaCount(entry));
      table.append(row);
    }
    details.append(table);
    section.append(details);
  }
  return section;
}

function sequenceStrip(slice, window = null, comparison = null) {
  const strip = node("code", "sequence-strip");
  strip.append(node("span", "sequence-trim", slice.clippedLeft ? "…" : " "), document.createTextNode(slice.left));
  let column = (window?.offset || 0) + slice.from;
  for (const piece of slice.pieces) {
    const segment = node(piece.matched ? "mark" : "span", "");
    for (const residue of piece.text) {
      const difference = comparison?.positions.get(column++);
      if (difference) {
        const aminoAcid = node("span", "sequence-difference", residue);
        aminoAcid.title = difference;
        segment.append(aminoAcid);
      } else segment.append(document.createTextNode(residue));
    }
    strip.append(segment);
  }
  strip.append(document.createTextNode(slice.right), node("span", "sequence-trim", slice.clippedRight ? "…" : " "));
  return strip;
}

function sequenceRuler(start, end, origin) {
  const columns = Array(end - start).fill(" ");
  const ticks = Array(end - start).fill(" ");
  for (let column = start; column < end; column++) {
    const position = column - origin + 1;
    if (column !== start && position % 10 !== 0) continue;
    const label = String(position);
    const index = column - start;
    const labelStart = Math.max(0, index - label.length + 1);
    if (labelStart + label.length <= columns.length) {
      [...label].forEach((digit, i) => { columns[labelStart + i] = digit; });
      ticks[index] = "│";
    }
  }
  return node("code", "sequence-ruler", ` ${columns.join("")} \n ${ticks.join("")} `);
}

function referenceRow(report, entry, slice, part, filters) {
  const { row } = entry;
  const item = node("tr", "alignment-reference");
  item.id = `${row.peptide.peptide_id}-part-${part}`;
  const label = node("th", "alignment-label", `P${row.variant.published_vaccine_peptides.indexOf(row.peptide) + 1}`);
  label.scope = "row";
  label.title = `${row.peptide.in_vaccines.join(", ")} · ${row.peptide.sequence.length} aa · ${row.peptide.is_mrna_minimal_epitope ? "Minimal epitope" : "Vaccine peptide"}`;
  const sequence = node("td", "alignment-sequence");
  sequence.append(sequenceStrip(slice));
  const match = node("td", `alignment-match ${row.status}`, statusSymbols[row.status]);
  match.title = statusText[row.status];
  match.setAttribute("aria-label", statusText[row.status]);
  item.append(label, sequence, match, node("td"), ...supportCells(report, row.matches, filters, row.variant));
  return item;
}

function windowMetadata(window) {
  const sequences = new Set(window.members.map(({ protein }) => protein.sequence));
  const lengths = [...new Set([...sequences].map((seq) => seq.replace(/\*$/, "").length))].sort((a, b) => a - b);
  const outputs = new Set([...window.observations.values()].map((hit) => hit.reconstruction_id)).size;
  return `${outputs} protein output${outputs === 1 ? "" : "s"}; ${sequences.size} full sequence${sequences.size === 1 ? "" : "s"}; lengths ${lengths.join(", ")} aa`;
}

function differenceDescription(comparison) {
  const changes = comparison.substitutions.map(({ position, reference, residue }) => `${reference}${position}${residue}`);
  if (comparison.missing) changes.push(`${comparison.missing} fewer covered positions`);
  if (comparison.additional) changes.push(`${comparison.additional} additional covered positions`);
  return changes.join(", ");
}

function reconstructionRow(report, window, slice, filters, rows, comparison) {
  const item = node("tr", "alignment-reconstruction");
  const label = node("th", "alignment-label", `W${window.number}`);
  label.scope = "row";
  label.title = windowMetadata(window);
  const sequence = node("td", "alignment-sequence");
  sequence.append(sequenceStrip(slice, window, comparison));
  const contained = rows.filter(({ peptide }) => window.peptideIds.has(peptide.peptide_id)).map(({ variant, peptide }) =>
    `P${variant.published_vaccine_peptides.indexOf(peptide) + 1}`);
  const match = node("td", `alignment-match ${contained.length ? "contained" : "sequence_disagreement"}`, contained.join(" ") || "≠");
  match.title = contained.length ? `Contains ${contained.join(", ")}` : "No contained vaccine sequence";
  const changes = node("td", "alignment-differences");
  if (comparison) {
    if (window.number === comparison.reference) changes.append(node("span", "comparison-reference", "ref"));
    else {
      const description = `vs W${comparison.reference}: ${differenceDescription(comparison)}`;
      changes.title = description;
      changes.setAttribute("aria-label", description);
      const { substitutions, missing, additional } = comparison;
      if (substitutions.length) changes.append(node("span", "difference-badge", substitutions.length === 1
        ? differenceDescription({ ...comparison, missing: 0, additional: 0 }) : `${substitutions.length} changes`));
      if (missing || additional) changes.append(node("span", "coverage-badge", "coverage"));
    }
  } else changes.append(node("span", "support-empty", "—"));
  item.append(label, sequence, match, changes, ...supportCells(report, [...window.observations.values()], filters));
  return item;
}

function targetComparison(report, rows, filters, columns) {
  const fragment = document.createDocumentFragment();
  const groups = sequenceComparison(report, rows, filters);
  const samples = rnaSupport(report, [], filters);
  let nextWindow = 1;
  for (const group of groups) {
    group.windows = reconstructionWindows(group, group.start, group.end, nextWindow);
    group.comparisons = group.unaligned ? null : windowDifferences(group.windows, group.start, group.end);
    nextWindow += group.windows.length;
  }
  const reconstructionIds = new Set(groups.flatMap((group) => group.reconstructions.flatMap((p) => [...p.observations.values()].map((hit) => hit.reconstruction_id))));
  // Reference and reconstruction support share the same sample columns.
  groups.forEach((group, groupIndex) => {
    if (group.unaligned) fragment.append(node("p", "sequence-alignment-note", "Unaligned candidates · mutation-centered windows; no shared vaccine anchor."));
    else if (groups.filter((g) => !g.unaligned).length > 1) fragment.append(node("p", "sequence-alignment-note", `Alignment group ${groupIndex + 1}`));
    const nParts = Math.ceil((group.end - group.start) / columns);
    for (let start = group.start, part = 1; start < group.end; start += columns, part++) {
      const end = Math.min(start + columns, group.end);
      const panel = node("div", "sequence-table-scroll");
      panel.tabIndex = 0;
      panel.setAttribute("role", "region");
      panel.setAttribute("aria-label", `${rows[0].variant.gene} sequence alignment${nParts > 1 ? `, region ${part}` : ""}`);
      const table = node("table", "sequence-alignment");
      const caption = node("caption", "sequence-alignment-note", `${group.end - group.start} aa${nParts > 1 ? ` · region ${part}/${nParts}` : ""}`);
      const head = node("thead"), header = node("tr");
      header.append(node("th", "", "Seq"), node("th", "", "Sequence"), node("th", "", "Match"),
        node("th", "", group.unaligned || !group.windows.length ? "Changes" : `Δ vs W${group.windows[0].number}`));
      samples.forEach(({ sample }) => header.append(node("th", "sample-heading", sample)));
      const ruler = node("tr", "alignment-ruler");
      const scale = node("td", "alignment-sequence");
      scale.append(sequenceRuler(start, end, group.start));
      ruler.append(node("td"), scale);
      const rest = node("td");
      rest.colSpan = 2 + samples.length;
      ruler.append(rest);
      head.append(header, ruler);
      const body = node("tbody");
      for (const entry of group.references) {
        const slice = alignmentSlice(entry.row.peptide.sequence, entry.offset, start, end);
        if (slice) body.append(referenceRow(report, entry, slice, part, filters));
      }
      for (const window of group.windows) {
        const slice = windowSlice(window, start, end);
        if (slice) body.append(reconstructionRow(report, window, slice, filters, rows, group.comparisons?.get(window)));
      }
      table.append(caption, head, body);
      panel.append(table);
      fragment.append(panel);
    }
  });
  const details = node("details", "sequence-details sequence-protein-details");
  details.append(node("summary", "", `Peptide & protein details (${rows.length} peptides, ${reconstructionIds.size} outputs)`));
  const metadata = node("dl", "sequence-metadata");
  for (const { variant, peptide } of rows) {
    metadata.append(node("dt", "", `P${variant.published_vaccine_peptides.indexOf(peptide) + 1}`),
      node("dd", "", `${peptide.in_vaccines.join(", ")} · ${peptide.sequence.length} aa · ${peptide.is_mrna_minimal_epitope ? "Minimal epitope" : "Vaccine peptide"}`));
  }
  for (const group of groups) for (const window of group.windows) {
    const comparison = group.comparisons?.get(window);
    metadata.append(node("dt", "", `W${window.number}`), node("dd", "", windowMetadata(window) +
      (comparison && comparison.reference !== window.number ? `. vs W${comparison.reference}: ${differenceDescription(comparison)}` : "")));
  }
  details.append(metadata);
  fragment.append(details);
  if (!reconstructionIds.size) {
    const explanations = {
      peptide: "Translated candidates were produced, but none contains these vaccine sequences.",
      proteoform: "Translated candidates were produced, but none contains these vaccine sequences.",
      rna_only: "Variant RNA called; no mutant protein translated.",
      no_call: "No target RNA variant called.",
      no_reads: "No reads covered this target.",
    };
    const reason = rows.every((row) => row.status === "not_evaluated") ? "No completed output for this selection."
      : filters.sample || filters.method ? "No reconstructions in this selection."
      : explanations[rows[0].variant.candidate_outcome] || "No reconstructions in the available outputs.";
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
  queryNode("#sequence-count").textContent = `${rows.length} peptides across ${groups.size} target${groups.size === 1 ? "" : "s"} shown`;
  const container = queryNode("#sequence-targets");
  container.replaceChildren();
  const columns = Math.max(45, Math.min(72, Math.floor((container.clientWidth - 650) / 9.6) - 2));
  for (const entries of groups.values()) {
    const variant = entries[0].variant;
    const target = node("article", "sequence-target");
    const title = node("h3", "sequence-target-title");
    title.append(node("span", "", variant.gene), node("span", "mono", variant.protein_change || "Protein change not annotated"));
    target.append(title);
    target.append(targetComparison(report, entries, filters, columns));
    target.append(sourceRnaEvidence(report, variant, filters));
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
  queryNode("#sequence-provenance").append(link, node("span", "", ` · ${analysis.n_completed_methods} completed methods.`));
  if (report.source?.snapshot?.name) queryNode("#sequence-provenance").append(node("span", "",
    ` · Sid snapshot ${report.source.snapshot.name}.`));
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

if (typeof module !== "undefined") module.exports = { sequenceRows, highlightedPieces, sharedOffset, referenceGroups, sequenceComparison, alignmentSlice, reconstructionWindows, windowDifferences, windowSlice, rnaSupport, sequenceStatus, sourceRnaState };
if (typeof document !== "undefined") main().catch((error) => {
  queryNode("#sequence-summary").textContent = `The sequence report could not be loaded (${error.message}). Reload the page to try again.`;
});
