const OUTCOMES = {
  peptide:    { label: "peptide",    short: "PEP", tone: "ok",   blurb: "mutant protein translated and novel peptides called" },
  proteoform: { label: "proteoform", short: "PRO", tone: "ok",   blurb: "mutant protein sequence translated" },
  rna_only:   { label: "RNA only",   short: "RNA", tone: "warn", blurb: "RNA variant called, no mutant protein" },
  no_call:    { label: "no call",    short: "—",   tone: "bad",  blurb: "reads present, nothing called" },
  no_reads:   { label: "no reads",   short: "·",   tone: "none", blurb: "locus not covered by the ONT data" },
};
const OUTCOME_ORDER = ["peptide", "proteoform", "rna_only", "no_call", "no_reads"];
const RECOVERED = new Set(["peptide", "proteoform"]);

const TONE_VAR = { ok: "--ok", warn: "--warn", bad: "--bad", none: "--none" };
const SEVERITY_TONE = {
  crash: "bad", "silent data loss": "bad",
  performance: "warn", limitation: "warn",
};

let DATA = null;
let sortKey = "gene";
let sortAsc = true;

const $ = (selector) => document.querySelector(selector);
const el = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};

/* ---------------------------------------------------------------- helpers */

function outcomeOf(variant, timepoint) {
  const recovery = variant.recovery;
  if (!recovery) return null;
  if (timepoint) return recovery.timepoints?.[timepoint]?.outcome ?? "no_reads";
  return recovery.outcome ?? "no_reads";
}

function ontReads(variant) {
  return Object.values(variant.ont_expectation || {})
    .reduce((total, entry) => total + (entry?.total_reads || 0), 0);
}

function latestVaf(variant) {
  const trend = variant.vaf_trend || [];
  return trend.length ? trend[trend.length - 1].value : -1;
}

function sparkline(trend) {
  if (!trend || trend.length < 2) {
    return el("span", "locus", trend?.length ? trend[0].value.toFixed(3) : "—");
  }
  const width = 62, height = 18, pad = 2;
  const values = trend.map((point) => point.value);
  const max = Math.max(...values, 0.001);
  const step = (width - pad * 2) / (values.length - 1);
  const points = values.map((value, index) => [
    pad + index * step,
    height - pad - (value / max) * (height - pad * 2),
  ]);
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("class", "spark");
  svg.setAttribute("width", width);
  svg.setAttribute("height", height);
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  path.setAttribute("d", points.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`).join(" "));
  svg.appendChild(path);
  const last = document.createElementNS("http://www.w3.org/2000/svg", "circle");
  last.setAttribute("cx", points.at(-1)[0].toFixed(1));
  last.setAttribute("cy", points.at(-1)[1].toFixed(1));
  last.setAttribute("r", "1.9");
  svg.appendChild(last);
  const wrapper = el("span");
  wrapper.title = trend.map((point) => `${point.label}: ${point.value}`).join("\n");
  wrapper.appendChild(svg);
  return wrapper;
}

function badge(outcome) {
  const spec = OUTCOMES[outcome] || OUTCOMES.no_reads;
  return el("span", `badge ${spec.tone}`, spec.label);
}

/* ------------------------------------------------------------- top matter */

function renderVerdict() {
  const node = $("#verdict");
  const summary = DATA.summary;

  if (!DATA.has_exacto_run) {
    node.style.setProperty("--verdict-color", `var(--none)`);
    node.append(
      el("span", "answer", "Not yet run"),
      el("span", "detail", `${summary.n_variants} vaccine variants catalogued; no Exacto run recorded.`),
    );
    return;
  }

  const recovered = summary.n_recovered ?? 0;
  const testable = summary.n_testable ?? 0;
  const answer = recovered === 0 ? "No"
    : recovered === testable ? "Yes"
    : "Partly";
  const tone = recovered === 0 ? "bad" : recovered === testable ? "ok" : "warn";
  node.style.setProperty("--verdict-color", `var(${TONE_VAR[tone]})`);

  const version = summary.exacto_version ? ` · Exacto ${summary.exacto_version}` : "";
  const residues = summary.n_residue_checkable
    ? `, ${summary.n_residue_confirmed} of ${summary.n_residue_checkable} with the right amino acid`
    : "";
  node.append(
    el("span", "answer", answer),
    el("span", "detail",
      `${recovered} of ${testable} covered vaccine mutations came back as a translated ` +
      `mutant protein sequence${residues}${version}`),
  );
}

function renderTiles() {
  const node = $("#tiles");
  node.innerHTML = "";
  const counts = DATA.summary.outcome_counts;

  const total = DATA.variants.length;
  const recovered = DATA.summary.n_recovered ?? 0;
  const testable = DATA.summary.n_testable ?? total;

  const tiles = [
    {
      label: "Mutant proteins recovered",
      value: DATA.has_exacto_run ? `${recovered}/${testable}` : "—",
      sub: DATA.has_exacto_run ? "of mutations the ONT data covers" : "no run yet",
      bar: DATA.has_exacto_run && counts ? OUTCOME_ORDER.map((outcome) => ({
        outcome, n: counts[outcome] || 0,
      })) : null,
    },
    {
      label: "Vaccine variants",
      value: String(total),
      sub: `${DATA.summary.n_peptide_entries} peptide entries across ${DATA.vaccine_names.length} vaccines`,
    },
  ];

  if (DATA.summary.n_residue_checkable) {
    tiles.splice(1, 0, {
      label: "…with the right amino acid",
      value: `${DATA.summary.n_residue_confirmed ?? 0}/${DATA.summary.n_residue_checkable}`,
      sub: "of the recovered proteoforms, vs. the annotated change",
    });
  }
  if (DATA.summary.n_with_vaccine_epitopes) {
    tiles.splice(2, 0, {
      label: "Exact vaccine peptides found",
      value: `${DATA.summary.n_epitope_confirmed ?? 0}/${DATA.summary.n_with_vaccine_epitopes}`,
      sub: "proteoform literally contains the manufactured epitope",
    });
  }

  for (const timepoint of DATA.timepoints) {
    const extraction = DATA.extraction?.[timepoint.name];
    const outcomes = DATA.variants.map((variant) => outcomeOf(variant, timepoint.name));
    const hits = outcomes.filter((outcome) => RECOVERED.has(outcome)).length;
    tiles.push({
      label: `${timepoint.name} · ${timepoint.biopsy_date}`,
      value: DATA.has_exacto_run ? String(hits) : "—",
      sub: extraction
        ? `${extraction.n_reads.toLocaleString()} reads · ${extraction.n_spanning_reads.toLocaleString()} spanning`
        : "not extracted",
    });
  }

  for (const spec of tiles) {
    const tile = el("div", "tile");
    tile.append(el("div", "label", spec.label), el("div", "value", spec.value));
    if (spec.sub) tile.append(el("div", "sub", spec.sub));
    if (spec.bar) {
      const bar = el("div", "bar");
      const sum = spec.bar.reduce((acc, item) => acc + item.n, 0) || 1;
      for (const item of spec.bar) {
        if (!item.n) continue;
        const segment = el("span");
        segment.style.width = `${(item.n / sum) * 100}%`;
        segment.style.background = `var(${TONE_VAR[OUTCOMES[item.outcome].tone]})`;
        segment.title = `${item.n} ${OUTCOMES[item.outcome].label}`;
        bar.appendChild(segment);
      }
      tile.append(bar);
    }
    node.appendChild(tile);
  }
}

function renderLegend() {
  const node = $("#legend");
  node.innerHTML = "";
  for (const outcome of OUTCOME_ORDER) {
    const spec = OUTCOMES[outcome];
    const item = el("span");
    const swatch = el("span", "swatch");
    swatch.style.background = `var(${TONE_VAR[spec.tone]})`;
    item.append(swatch, document.createTextNode(`${spec.label} — ${spec.blurb}`));
    node.appendChild(item);
  }
}

function renderHistory() {
  const history = DATA.history || [];
  if (history.length < 2) return;
  $("#history-block").hidden = false;
  const node = $("#history");
  node.innerHTML = "";
  const max = Math.max(...history.map((entry) => entry.n_testable || 1), 1);
  for (const entry of history.slice().reverse()) {
    const row = el("div", "history-row");
    const track = el("div", "track");
    const fill = el("span");
    fill.style.width = `${((entry.n_recovered || 0) / max) * 100}%`;
    track.appendChild(fill);
    row.append(
      el("span", "date", entry.date),
      el("span", "ver", entry.exacto_version || "—"),
      track,
      el("span", "locus", `${entry.n_recovered}/${entry.n_testable}`),
    );
    node.appendChild(row);
  }
}

/* ----------------------------------------------------------- variant table */

const BASE_COLUMNS = [
  { key: "gene", label: "Gene", sticky: true },
  { key: "protein", label: "Protein change", sticky: true },
  { key: "locus", label: "Locus (GRCh38)" },
  { key: "type", label: "Change" },
  { key: "vaccines", label: "Vaccines" },
  { key: "elispot", label: "ELISPOT" },
];

function renderHead() {
  const groups = $("#assay-group-row");
  const cols = $("#column-row");
  groups.innerHTML = "";
  cols.innerHTML = "";

  for (const col of BASE_COLUMNS) {
    const cell = el("th", col.sticky ? "sticky-col" : null);
    cell.rowSpan = 2;
    cell.dataset.sort = col.key;
    cell.textContent = col.label;
    groups.append(cell);
  }

  // One group header per assay, spanning its timepoints.
  let run = [];
  const flush = () => {
    if (!run.length) return;
    const spec = run[0];
    const cell = el("th", `assay-group ${spec.kind}`);
    cell.colSpan = run.length;
    cell.append(el("div", "assay-name", spec.assay_label));

    // Every column says its platform, its read length, and whether it is
    // single-cell — none of it left to be inferred from the name.
    const tags = el("div", "assay-tags");
    tags.append(el("span", `tag platform ${spec.long_read ? "longread" : ""}`, spec.platform));
    tags.append(el("span", `tag ${spec.long_read ? "longread" : ""}`, spec.read_length));
    if (spec.single_cell) tags.append(el("span", "tag sc", "single cell"));
    cell.append(tags);
    if (spec.tested) cell.append(el("span", "badge ok tested", "tested by Exacto"));
    groups.append(cell);
    for (const col of run) {
      const sub = el("th", `assay-col ${col.kind} num`);
      sub.dataset.sort = `assay:${col.key}`;
      sub.textContent = col.timepoint;
      cols.append(sub);
    }
    run = [];
  };
  for (const col of DATA.assay_columns || []) {
    if (run.length && run[0].assay !== col.assay) flush();
    run.push(col);
  }
  flush();

  const verdict = el("th", "sticky-right");
  verdict.rowSpan = 2;
  verdict.dataset.sort = "outcome";
  verdict.textContent = "Exacto";
  groups.append(verdict);
}

function assayCell(variant, column) {
  const cell = el("td", `num assay-cell ${column.kind}`);
  const data = variant.assay_matrix?.[column.key];
  if (!data || data.vaf === null || data.vaf === undefined) {
    cell.append(el("span", "locus", "—"));
    return cell;
  }
  const vaf = el("div", "vaf-value", data.vaf.toFixed(3));
  if (data.vaf === 0) vaf.classList.add("zero");
  else if (data.vaf >= 0.05) vaf.classList.add("present");
  cell.append(vaf);
  cell.append(el("div", "vaf-reads", `${data.alt.toLocaleString()}/${data.total.toLocaleString()}`));
  const germline = variant.germline_matrix?.[column.key];
  const descriptor = [column.platform, column.read_length]
    .concat(column.single_cell ? ["single cell"] : []).join(", ");
  cell.title =
    `${column.assay_label} ${column.timepoint} (${descriptor}): ` +
    `${data.alt} alt / ${data.total} total across ${data.samples} sample(s)` +
    (germline && germline.total
      ? `\nmatched normal: ${germline.alt}/${germline.total} (VAF ${(germline.vaf ?? 0).toFixed(3)})`
      : "");
  return cell;
}

function sortValue(variant, key) {
  if (key.startsWith("assay:")) {
    const cell = variant.assay_matrix?.[key.slice(6)];
    return cell && cell.vaf !== null && cell.vaf !== undefined ? -cell.vaf : 1;
  }
  switch (key) {
    case "gene": return variant.gene;
    case "protein": return variant.protein_change || "";
    case "locus": return [variant.chrom, String(variant.pos).padStart(12, "0")].join(":");
    case "type": return variant.variant_type;
    case "vaccines": return -variant.vaccines.length;
    case "elispot": return variant.elispot.status;
    case "vaf": return -latestVaf(variant);
    case "reads": return -ontReads(variant);
    case "outcome": {
      const outcome = outcomeOf(variant, null);
      return outcome ? OUTCOME_ORDER.indexOf(outcome) : 99;
    }
    default: return variant.gene;
  }
}

function visibleVariants() {
  const query = $("#filter").value.trim().toLowerCase();
  const onlyRecovered = $("#only-recovered").checked;
  const onlyElispot = $("#only-elispot").checked;

  return DATA.variants
    .filter((variant) => {
      if (onlyRecovered && !RECOVERED.has(outcomeOf(variant, null))) return false;
      if (onlyElispot && variant.elispot.status !== "positive") return false;
      if (!query) return true;
      const haystack = [
        variant.gene, variant.protein_change, variant.vaccine_label,
        `${variant.chrom}:${variant.pos}`, variant.consequence,
        variant.vaccines.join(" "),
      ].join(" ").toLowerCase();
      return haystack.includes(query);
    })
    .sort((a, b) => {
      const left = sortValue(a, sortKey);
      const right = sortValue(b, sortKey);
      const cmp = typeof left === "number" ? left - right : String(left).localeCompare(String(right));
      return sortAsc ? cmp : -cmp;
    });
}

function timepointCells(variant) {
  const wrapper = el("div", "tp-cells");
  for (const timepoint of DATA.timepoints) {
    const outcome = outcomeOf(variant, timepoint.name);
    const spec = OUTCOMES[outcome] || OUTCOMES.no_reads;
    const cell = el("span", "tp-cell", DATA.has_exacto_run ? spec.short : "·");
    cell.style.background = `var(${TONE_VAR[spec.tone]}-soft)`;
    cell.style.color = `var(${TONE_VAR[spec.tone]})`;
    cell.title = `${timepoint.name}: ${spec.label}`;
    wrapper.appendChild(cell);
  }
  return wrapper;
}

function highlightedSequence(form, expected) {
  const box = el("div", "seq");
  const context = form.context || "";
  const indices = form.mutant_residue_indices || [];
  const start = (form.context_start || 1) - 1;
  for (let offset = 0; offset < context.length; offset += 1) {
    const absolute = start + offset;
    const residue = context[offset];
    if (indices.includes(absolute)) {
      box.appendChild(el("span", "hit", residue));
    } else {
      box.appendChild(document.createTextNode(residue));
    }
  }
  if (expected?.kind === "missense" && expected.alt_aa) {
    box.appendChild(el("span", "locus", `  (expecting ${expected.alt_aa})`));
  }
  return box;
}

function detailCard(timepoint, variant) {
  const card = el("div", "detail-card");
  card.append(el("h4", null, `${timepoint.name} · ${timepoint.biopsy_date}`));

  const recovery = variant.recovery?.timepoints?.[timepoint.name];
  if (!recovery || !Object.keys(recovery.arms || {}).length) {
    card.append(el("div", "locus", "no run recorded"));
    return card;
  }

  for (const arm of DATA.arms) {
    const entry = recovery.arms[arm];
    if (!entry) continue;
    const list = el("dl", "kv");
    const add = (key, value) => {
      list.append(el("dt", null, key), el("dd", null, value));
    };
    add("arm", arm);
    add("outcome", (OUTCOMES[entry.outcome] || OUTCOMES.no_reads).label);
    const available = entry.spanning_reads_available ?? entry.spanning_reads;
    add(
      "spanning reads",
      available > entry.spanning_reads
        ? `${entry.spanning_reads.toLocaleString()} of ${available.toLocaleString()}`
        : entry.spanning_reads.toLocaleString(),
    );
    add(
      "RNA calls",
      `${entry.rna_variant_calls.length} exact` +
        (entry.integrated_rna_call_ids?.length
          ? ` · ${entry.integrated_rna_call_ids.length} integrated`
          : ""),
    );
    add("proteoforms", String(entry.n_proteoforms));
    if (entry.n_mutant_peptides) add("novel peptides", String(entry.n_mutant_peptides));
    if (entry.residue_confirmed !== null && entry.residue_confirmed !== undefined) {
      add("expected residue", entry.residue_confirmed ? "confirmed" : "not seen");
    }
    card.append(list);

    const form = entry.proteoforms?.[0];
    if (form) {
      card.append(highlightedSequence(form, entry.expected));
      const meta = el("div", "locus",
        `${form.reference_transcript_ids || "novel transcript"} · ${form.protein_length} aa` +
        (form.frameshift ? " · frameshift" : ""));
      card.append(meta);
    }
    if (entry.matched_epitopes?.length) {
      card.append(el("div", "locus", "vaccine epitopes found in this proteoform:"));
      const found = el("div", "peptide-list");
      for (const hit of entry.matched_epitopes.slice(0, 6)) {
        const chip = el("span", "epitope", hit.sequence);
        chip.title = `MHC class ${hit.mhc_class} · ${(hit.alleles || []).join(", ")}`;
        found.appendChild(chip);
      }
      card.append(found);
    }
    if (entry.mutant_peptides?.length) {
      const peptides = el("div", "peptide-list");
      for (const peptide of entry.mutant_peptides.slice(0, 6)) {
        peptides.appendChild(el("span", null, peptide.sequence));
      }
      card.append(peptides);
    }
    card.append(el("hr", null));
  }
  card.querySelector("hr:last-of-type")?.remove();
  return card;
}

function detailRow(variant, columns) {
  const row = el("tr", "detail");
  const cell = el("td");
  cell.colSpan = columns;
  const inner = el("div", "detail-inner");

  const facts = el("div", "detail-card");
  facts.append(el("h4", null, "Variant"));
  const list = el("dl", "kv");
  const add = (key, value) => list.append(el("dt", null, key), el("dd", null, value));
  add("locus", `${variant.chrom}:${variant.pos.toLocaleString()}`);
  add("change", `${variant.ref} → ${variant.alt}`);
  add("consequence", variant.consequence || "—");
  add("impact", variant.impact || "—");
  add("vaccines", variant.vaccines.join(", "));
  if (variant.peptide_classes?.length) add("peptide", variant.peptide_classes.join(", "));
  add("ELISPOT", variant.elispot.status.replace("_", " "));
  if (variant.vaccine_epitopes?.length) {
    const found = variant.recovery?.matched_epitopes?.length || 0;
    add("vaccine epitopes", `${variant.vaccine_epitopes.length} published, ${found} found`);
  }
  for (const [name, entry] of Object.entries(variant.ont_expectation || {})) {
    add(`${name} ONT depth`, entry ? `${entry.alt_reads}/${entry.total_reads} (VAF ${entry.vaf ?? "—"})` : "—");
  }
  facts.append(list);

  const grid = el("div", "detail-grid");
  grid.append(facts);
  for (const timepoint of DATA.timepoints) grid.append(detailCard(timepoint, variant));
  inner.append(grid);
  cell.append(inner);
  row.append(cell);
  return row;
}

function renderTable() {
  const body = $("#variant-table tbody");
  body.innerHTML = "";
  const rows = visibleVariants();
  const columns = BASE_COLUMNS.length + (DATA.assay_columns || []).length + 1;

  if (!rows.length) {
    const row = el("tr");
    const cell = el("td", "empty", "Nothing matches that filter.");
    cell.colSpan = columns;
    row.append(cell);
    body.append(row);
    return;
  }

  for (const variant of rows) {
    const row = el("tr", "row");

    row.append(el("td", "gene sticky-col", variant.gene));
    row.append(el("td", "change sticky-col", variant.protein_change || variant.vaccine_label || "—"));
    row.append(el("td", "locus", `${variant.chrom}:${variant.pos.toLocaleString()}`));

    const change = el("td", "change");
    const shortRef = variant.ref.length > 8 ? `${variant.ref.slice(0, 6)}…` : variant.ref;
    const shortAlt = variant.alt.length > 8 ? `${variant.alt.slice(0, 6)}…` : variant.alt;
    change.textContent = `${shortRef}>${shortAlt}`;
    change.title = `${variant.ref}>${variant.alt} (${variant.variant_type})`;
    row.append(change);

    const vaccines = el("td");
    const chips = el("div", "chips");
    for (const name of variant.vaccines) chips.append(el("span", "chip", name));
    vaccines.append(chips);
    row.append(vaccines);

    const elispot = el("td");
    if (variant.elispot.status === "positive") elispot.append(el("span", "badge ok", "responded"));
    else if (variant.elispot.status === "negative") elispot.append(el("span", "badge warn", "no response"));
    else elispot.append(el("span", "locus", "untested"));
    row.append(elispot);

    for (const column of DATA.assay_columns || []) row.append(assayCell(variant, column));

    const outcome = el("td", "sticky-right");
    outcome.append(timepointCells(variant));
    if (variant.recovery?.residue_confirmed === false) {
      const warn = el("span", "badge warn", "wrong residue");
      warn.style.marginLeft = ".35rem";
      warn.title = "A mutant proteoform came back, but not carrying the annotated change";
      outcome.append(warn);
    }
    row.append(outcome);

    body.append(row);

    let expanded = null;
    row.addEventListener("click", () => {
      if (expanded) {
        expanded.remove();
        expanded = null;
        row.classList.remove("open");
      } else {
        expanded = detailRow(variant, columns);
        row.after(expanded);
        row.classList.add("open");
      }
    });
  }
}

/* -------------------------------------------------------------- findings */

function renderFindings() {
  const node = $("#finding-list");
  node.innerHTML = "";
  const findings = DATA.findings || [];
  if (!findings.length) {
    node.append(el("div", "empty", "Nothing outstanding."));
    return;
  }
  for (const finding of findings) {
    const card = el("article", "finding");
    const head = el("div", "finding-head");
    head.append(
      el("h3", null, finding.title),
      el("span", `badge ${SEVERITY_TONE[finding.severity] || "warn"}`, finding.severity),
      el("span", "locus", `Exacto ${finding.observed_in}`),
    );
    card.append(head);
    if (finding.where) card.append(el("div", "finding-where mono", finding.where));
    card.append(el("p", null, finding.detail));
    if (finding.workaround) {
      const line = el("p", "finding-fix");
      line.append(el("strong", null, "Worked around: "), document.createTextNode(finding.workaround));
      card.append(line);
    }
    if (finding.suggested_fix) {
      const line = el("p", "finding-fix");
      line.append(el("strong", null, "Suggested fix: "), document.createTextNode(finding.suggested_fix));
      card.append(line);
    }
    node.append(card);
  }
}

/* ------------------------------------------------------------- run cards */

function renderRuns() {
  const node = $("#run-cards");
  node.innerHTML = "";
  if (!DATA.runs?.length) {
    node.append(el("div", "empty", "No Exacto run has been recorded yet."));
    return;
  }
  for (const run of DATA.runs) {
    const card = el("div", `run-card ${run.status}`);
    card.append(el("h3", null, `${run.timepoint} · ${run.arm}`));
    const counts = Object.entries(run.counts || {})
      .map(([key, value]) => `${key.replace(/_/g, " ")} ${value.toLocaleString()}`)
      .join(" · ");
    card.append(el("div", "meta",
      `${run.status} · ${Math.round((run.seconds || 0) / 60)} min${counts ? ` · ${counts}` : ""}`));

    const steps = el("ul", "steps");
    for (const step of run.steps || []) {
      const item = el("li", "step-row");
      const line = el("div", "step-line");
      line.append(el("span", `name${step.returncode ? " fail" : ""}`, step.name));
      line.append(el("span", "locus", `${step.seconds}s`));
      item.append(line);
      if (step.command?.length) {
        const detail = el("details", "step-command");
        detail.append(el("summary", "locus", "command"));
        detail.append(copyBlock(step.command.join(" \\\n  ")));
        item.append(detail);
      }
      steps.append(item);
    }
    card.append(steps);
    for (const note of run.workarounds || []) {
      card.append(el("div", "run-note", `worked around: ${note}`));
    }
    if (run.error) card.append(el("div", "run-error", run.error));
    node.append(card);
  }
}

/* ------------------------------------------------------------------ boot */

function copyBlock(text, label = "Copy") {
  const wrap = el("div", "copy-block");
  wrap.append(el("pre", "block", text));
  const button = el("button", "button ghost copy-button", label);
  button.addEventListener("click", async () => {
    await navigator.clipboard.writeText(text);
    button.textContent = "Copied";
    setTimeout(() => { button.textContent = label; }, 1600);
  });
  wrap.append(button);
  return wrap;
}

function renderInventory() {
  const node = $("#inventory-block");
  if (!node) return;
  node.innerHTML = "";
  const inv = DATA.inventory;
  if (!inv) {
    node.append(el("div", "empty", "No inventory recorded."));
    return;
  }

  const timepoints = [...new Set(
    inv.grid.flatMap((row) => Object.keys(row.timepoints)))]
    .sort((a, b) => (a === "unlabelled" ? 1 : b === "unlabelled" ? -1 : a.localeCompare(b)));

  const scroll = el("div", "table-scroll");
  const table = el("table", "inventory-table");
  const head = el("thead");
  const headRow = el("tr");
  headRow.append(el("th", null, "Platform"));
  for (const tp of timepoints) headRow.append(el("th", "num", tp));
  head.append(headRow);
  table.append(head);

  const body = el("tbody");
  for (const row of inv.grid) {
    const tr = el("tr");
    tr.append(el("td", "gene", row.platform));
    for (const tp of timepoints) {
      const data = row.timepoints[tp];
      const cell = el("td", "num");
      if (!data) {
        cell.append(el("span", "locus", "—"));
      } else {
        cell.append(el("div", "vaf-value present", data.files.toLocaleString()));
        cell.title = data.areas.join("\n");
      }
      body.append;
      tr.append(cell);
    }
    body.append(tr);
  }
  table.append(body);
  scroll.append(table);
  node.append(scroll);
  node.append(el("p", "section-note", inv.note));
  const link = el("a", "source-url mono", inv.manifest_url);
  link.href = inv.manifest_url;
  node.append(el("div", null, `${inv.n_files.toLocaleString()} files listed in `));
  node.append(link);
}

function renderReproduction() {
  const node = $("#reproduction");
  if (!node) return;
  node.innerHTML = "";
  for (const stage of DATA.reproduction || []) {
    const block = el("div", "source-group");
    block.append(el("h3", null, stage.stage));
    block.append(copyBlock(stage.commands.join("\n")));
    node.append(block);
  }
}

/* --------------------------------------------------------- live run status */

const RUN_TONE = {
  in_progress: "warn", queued: "warn", pending: "warn", waiting: "warn",
  requested: "warn", success: "ok", failure: "bad", cancelled: "none",
  skipped: "none", timed_out: "bad",
};

async function github(path) {
  const response = await fetch(`https://api.github.com/repos/${DATA.repo}/${path}`, {
    headers: { Accept: "application/vnd.github+json" },
  });
  if (!response.ok) throw new Error(`GitHub API ${response.status}`);
  return response.json();
}

function since(iso) {
  const seconds = Math.max(0, (Date.now() - new Date(iso).getTime()) / 1000);
  if (seconds < 90) return `${Math.round(seconds)}s`;
  if (seconds < 5400) return `${Math.round(seconds / 60)} min`;
  return `${(seconds / 3600).toFixed(1)} h`;
}

async function renderLiveRun() {
  const node = $("#live-run");
  if (!node || !DATA.repo) return;

  const WAITING = ["queued", "pending", "waiting", "requested"];

  let runs;
  try {
    const payload = await github(
      `actions/workflows/${DATA.workflow_file}/runs?per_page=5`);
    runs = payload.workflow_runs || [];
  } catch {
    return; // Rate-limited or offline: the page is still fine without this.
  }
  if (!runs.length) return;

  // The newest run can be sitting in the queue while an older one does the
  // work — the workflow serialises itself so two runs cannot both publish.
  const active = runs.find((r) => r.status === "in_progress");
  const waiting = runs.filter((r) => WAITING.includes(r.status));
  const run = active || runs[0];

  const live = run.status !== "completed";
  const state = live ? run.status : run.conclusion;
  const tone = RUN_TONE[state] || "none";

  node.hidden = false;
  node.className = `live-run ${tone}${live ? " pulsing" : ""}`;
  node.innerHTML = "";

  const head = el("div", "live-head");
  head.append(el("span", `badge ${tone}`,
    active ? "running" : live ? "queued" : state));
  head.append(el("span", "live-title",
    active ? "An Exacto run is in progress"
      : live ? "An Exacto run is queued"
      : `Last run ${state}`));
  const link = el("a", "live-link", `run #${run.run_number} →`);
  link.href = run.html_url;
  link.target = "_blank";
  link.rel = "noreferrer";
  head.append(link);
  head.append(el("span", "locus",
    `${live ? "started" : "finished"} ${since(run.updated_at)} ago · ${run.head_branch} · ${run.head_sha.slice(0, 7)}`));
  if (waiting.length && active) {
    head.append(el("span", "locus",
      `· ${waiting.length} more queued behind it`));
  }
  node.append(head);

  // Per-job progress, so "where in the process" is answerable at a glance.
  try {
    const { jobs } = await github(`actions/runs/${run.id}/jobs?per_page=30`);
    if (!jobs?.length) {
      node.append(el("div", "live-note",
        "Waiting for a runner — the workflow runs one at a time so results "
        + "cannot be published out of order."));
      if (live) setTimeout(renderLiveRun, 20000);
      return;
    }
    const grid = el("div", "live-jobs");
    for (const job of jobs) {
      const jobState = job.status !== "completed" ? job.status : job.conclusion;
      const jobTone = RUN_TONE[jobState] || "none";
      const card = el("div", `live-job ${jobTone}`);
      card.append(el("div", "live-job-name", job.name));
      const current = (job.steps || []).find((s) => s.status === "in_progress")
        || [...(job.steps || [])].reverse().find((s) => s.status === "completed");
      const done = (job.steps || []).filter((s) => s.status === "completed").length;
      const total = (job.steps || []).length;
      card.append(el("div", "live-job-step",
        current ? current.name : jobState));
      if (total) {
        const bar = el("div", "bar");
        const fill = el("span");
        fill.style.width = `${Math.round((done / total) * 100)}%`;
        fill.style.background = `var(${TONE_VAR[jobTone]})`;
        bar.append(fill);
        card.append(bar);
        card.append(el("div", "locus", `step ${Math.min(done + 1, total)} of ${total}`));
      }
      const jobLink = el("a", "live-link", "log →");
      jobLink.href = job.html_url;
      jobLink.target = "_blank";
      jobLink.rel = "noreferrer";
      card.append(jobLink);
      grid.append(card);
    }
    node.append(grid);
  } catch {
    // Jobs endpoint is a second request; skip it rather than lose the header.
  }

  if (live) setTimeout(renderLiveRun, 20000);
}

/* ------------------------------------------------------- data & method page */

function renderSources() {
  const node = $("#source-groups");
  node.innerHTML = "";
  for (const group of DATA.data_sources || []) {
    const block = el("div", "source-group");
    const head = el("div", "source-head");
    head.append(el("h3", null, group.group), el("span", "locus", group.origin));
    block.append(head);
    for (const entry of group.entries) {
      const row = el("div", "source-row");
      const title = el("div", "source-title");
      title.append(el("span", "source-label", entry.label));
      if (entry.size) title.append(el("span", "badge none", entry.size));
      row.append(title);
      const link = el("a", "source-url mono", entry.url);
      link.href = entry.url;
      link.rel = "noreferrer";
      row.append(link);
      if (entry.detail) row.append(el("p", null, entry.detail));
      if (entry.taken) row.append(el("p", "source-taken", entry.taken));
      block.append(row);
    }
    node.append(block);
  }
}

function renderParameters() {
  const node = $("#parameter-groups");
  node.innerHTML = "";
  for (const group of DATA.parameters || []) {
    const block = el("div", "source-group");
    block.append(el("h3", null, group.step));
    const list = el("dl", "kv wide");
    for (const [key, value] of group.settings) {
      list.append(el("dt", null, key), el("dd", null, value));
    }
    block.append(list);
    node.append(block);
  }
}

function renderEnvironment() {
  const node = $("#environment-block");
  node.innerHTML = "";
  const env = DATA.environment || {};
  const list = el("dl", "kv wide");
  const rows = [
    ["Exacto", env.exacto_version ? `${env.exacto_version} (${env.exacto_source || "?"} ${env.exacto_ref || ""})` : "not recorded — no run yet"],
    ["samtools", env.samtools],
    ["minimap2", env.minimap2],
    ["RNA-Bloom2", env.rnabloom],
    ["Python", env.python],
    ["Platform", env.platform],
    ["GENCODE", "v44, matching the 10x refdata-gex-GRCh38-2024-A the source BAMs used"],
    ["Harness commit", DATA.summary.commit],
  ];
  for (const [key, value] of rows) {
    if (!value) continue;
    list.append(el("dt", null, key), el("dd", null, String(value)));
  }
  node.append(list);
}

/* ------------------------------------------------------------- bugs page */

function issueMarkdown(finding) {
  const lines = [
    `**Observed in:** Exacto ${finding.observed_in}`,
    `**Location:** \`${finding.where}\``,
    `**Severity:** ${finding.severity}`,
    "",
    finding.detail,
  ];
  if (finding.suggested_fix) lines.push("", `**Suggested fix:** ${finding.suggested_fix}`);
  if (finding.workaround) lines.push("", `**Worked around by:** ${finding.workaround}`);
  lines.push(
    "",
    "---",
    "Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), " +
    "running Exacto over the ONT long-read RNA-seq at https://osteosarc.com.",
  );
  return `## ${finding.title}\n\n${lines.join("\n")}\n`;
}

function renderBugs() {
  const findings = DATA.findings || [];
  const summary = $("#bug-summary");
  const crashes = findings.filter((f) => f.severity === "crash").length;
  summary.style.setProperty("--verdict-color", `var(${crashes ? "--bad" : "--ok"})`);
  summary.append(
    el("span", "answer", String(findings.length)),
    el("span", "detail",
      `issues found so far — ${crashes} of them stop the pipeline outright` +
      (DATA.summary.exacto_version ? ` on Exacto ${DATA.summary.exacto_version}` : "")),
  );

  const node = $("#bug-list");
  node.innerHTML = "";
  for (const finding of findings) {
    const card = el("article", "finding bug");
    const head = el("div", "finding-head");
    head.append(
      el("h3", null, finding.title),
      el("span", `badge ${SEVERITY_TONE[finding.severity] || "warn"}`, finding.severity),
      el("span", "locus", `Exacto ${finding.observed_in}`),
    );
    card.append(head);
    card.append(el("div", "finding-where mono", finding.where));
    card.append(el("p", null, finding.detail));
    if (finding.workaround) {
      const line = el("p", "finding-fix");
      line.append(el("strong", null, "Worked around: "), document.createTextNode(finding.workaround));
      card.append(line);
    }
    if (finding.suggested_fix) {
      const line = el("p", "finding-fix");
      line.append(el("strong", null, "Suggested fix: "), document.createTextNode(finding.suggested_fix));
      card.append(line);
    }

    const actions = el("div", "bug-actions");
    const copy = el("button", "button", "Copy as issue");
    copy.addEventListener("click", async () => {
      await navigator.clipboard.writeText(issueMarkdown(finding));
      copy.textContent = "Copied";
      setTimeout(() => { copy.textContent = "Copy as issue"; }, 1600);
    });
    const file = el("a", "button ghost", "Open issue form →");
    file.href = "https://github.com/pirl-unc/exacto/issues/new";
    file.target = "_blank";
    file.rel = "noreferrer";
    actions.append(copy, file);
    card.append(actions);
    node.append(card);
  }
}

function renderObservedFailures() {
  const node = $("#observed-list");
  node.innerHTML = "";
  const failures = [];
  for (const run of DATA.runs || []) {
    for (const step of run.steps || []) {
      if (step.returncode) failures.push({ run, step });
    }
  }
  if (!failures.length) {
    node.append(el("div", "empty",
      DATA.has_exacto_run
        ? "No step exited non-zero in the latest run."
        : "No Exacto run has been recorded yet."));
    return;
  }
  for (const { run, step } of failures) {
    const card = el("details", "failure");
    const summary = el("summary");
    summary.append(
      el("span", "mono", step.name),
      el("span", "badge bad", `exit ${step.returncode}`),
      el("span", "locus", `${run.timepoint} · ${run.arm}`),
    );
    card.append(summary);
    if (step.command?.length) {
      card.append(el("h4", null, "Command"));
      card.append(el("pre", "block", step.command.join(" \\\n  ")));
    }
    if (step.log_tail) {
      card.append(el("h4", null, "Output (tail)"));
      card.append(el("pre", "block", step.log_tail.trim()));
    }
    node.append(card);
  }
}

function renderProvenance() {
  const parts = [`Generated ${DATA.generated_at.replace("T", " ").replace("+00:00", " UTC")}`];
  if (DATA.summary.commit) parts.push(`commit ${DATA.summary.commit}`);
  if (DATA.environment?.exacto_version) parts.push(`Exacto ${DATA.environment.exacto_version}`);
  if (DATA.environment?.gencode) parts.push(`GENCODE ${DATA.environment.gencode}`);
  $("#provenance").textContent = parts.join(" · ");
}

function wireControls() {
  $("#filter").addEventListener("input", renderTable);
  $("#only-recovered").addEventListener("change", renderTable);
  $("#only-elispot").addEventListener("change", renderTable);

  document.querySelectorAll("#variant-table thead th[data-sort]").forEach((header) => {
    header.addEventListener("click", () => {
      const key = header.dataset.sort;
      if (sortKey === key) sortAsc = !sortAsc;
      else { sortKey = key; sortAsc = true; }
      document.querySelectorAll("#variant-table thead th[data-sort]")
        .forEach((other) => other.classList.remove("sorted", "asc"));
      header.classList.add("sorted");
      if (sortAsc) header.classList.add("asc");
      renderTable();
    });
  });
}

async function main() {
  DATA = await (await fetch("data.json")).json();
  const page = document.body.dataset.page || "summary";

  if (page === "sources") {
    renderSources();
    renderInventory();
    renderReproduction();
    renderParameters();
    renderEnvironment();
  } else if (page === "bugs") {
    renderBugs();
    renderObservedFailures();
  } else {
    $("#variant-count").textContent = String(DATA.summary.n_variants);
    renderVerdict();
    renderTiles();
    renderLegend();
    renderHistory();
    renderHead();
    wireControls();
    renderTable();
    renderFindings();
    renderRuns();
  }
  renderProvenance();
  renderLiveRun();
}

main().catch((error) => {
  document.querySelector("main").prepend(
    el("div", "empty", `Could not load data.json — ${error}`));
});
