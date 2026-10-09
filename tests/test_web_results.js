// Exercise the real rendering functions without a browser or DOM dependency.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const test = require("node:test");
const vm = require("node:vm");

class Element {
  constructor(tag) {
    this.tag = tag;
    this.children = [];
    this.text = "";
    this.style = { setProperty() {} };
  }
  set textContent(value) { this.text = String(value); this.children = []; }
  get textContent() { return this.text + this.children.map((child) => child.textContent).join(""); }
  set innerHTML(value) { this.text = value; this.children = []; }
  append(...children) { this.children.push(...children); }
  appendChild(child) { this.append(child); }
}

function renderer(data, globals = {}) {
  const nodes = {};
  const context = vm.createContext({
    document: {
      querySelector(selector) { return nodes[selector] ||= new Element("div"); },
      createElement(tag) { return new Element(tag); },
    },
    data,
    ...globals,
  });
  const source = fs.readFileSync(path.join(__dirname, "../web/app.js"), "utf8");
  // Omit only the async page bootstrap; run the actual helpers and renderers.
  vm.runInContext(source.slice(0, source.lastIndexOf("\nmain().catch")), context);
  vm.runInContext("DATA = data;", context);
  return { nodes, context };
}

function render(data, functions) {
  const { nodes, context } = renderer(data);
  vm.runInContext(functions.map((name) => `${name}();`).join(" "), context);
  return nodes;
}

function result(testable, recovered, arms = {}) {
  return {
    has_exacto_run: true,
    summary: { n_variants: 37, n_peptide_entries: 38, n_testable: testable,
      n_recovered: recovered, exacto_version: "0.4.6a1", outcome_counts: { no_reads: 37 } },
    variants: [{ recovery: { samples: { "T1-ONT": { arms } } } }],
    samples: [{ name: "T1-ONT", label: "T1 · ONT" }],
    vaccine_names: [], runs: [],
    history: [{ n_testable: 37, n_recovered: 22, date: "2026-09-28", exacto_version: "0.4.6a1" }],
  };
}

test("old failed-only 0/0 shows failure and the last evaluated result", () => {
  const data = result(0, 0);
  data.runs.push({ sample: "T1-ONT", status: "failed" });
  data.summary.n_with_vaccine_epitopes = 10;
  data.history.push({ n_testable: 0, n_recovered: 0, date: "2026-10-05" });
  const nodes = render(data, ["renderVerdict", "renderTiles", "renderHistory"]);
  assert.match(nodes["#verdict"].textContent, /^Run failed/);
  assert.match(nodes["#verdict"].textContent, /Last evaluated result: 22\/37/);
  assert.doesNotMatch(nodes["#verdict"].textContent, /0 of 0/);
  assert.match(nodes["#tiles"].textContent, /Mutant proteins recovered—no evaluable results/);
  assert.match(nodes["#tiles"].textContent, /T1 · ONT—/);
  assert.doesNotMatch(nodes["#tiles"].textContent, /Exact vaccine peptides found/);
  assert.doesNotMatch(nodes["#history"].textContent, /0\/0/);
  assert.match(nodes["#history"].textContent, /No evaluation/);
});

test("completed runs without coverage get a separate status", () => {
  const data = result(0, 0, { reads: { outcome: "no_reads" } });
  data.summary.evaluation_status = "no_coverage";
  const nodes = render(data, ["renderVerdict", "renderTiles"]);
  assert.match(nodes["#verdict"].textContent, /^No covered mutations/);
  assert.match(nodes["#tiles"].textContent, /no mutations covered/);
  assert.doesNotMatch(nodes["#verdict"].textContent, /0 of 0/);
});

test("a measured zero recovery remains a negative finding", () => {
  const nodes = render(result(37, 0), ["renderVerdict"]);
  assert.match(nodes["#verdict"].textContent, /^No0 of 37 covered/);
});

test("partial recovery keeps its measured fraction", () => {
  const nodes = render(result(37, 22), ["renderVerdict"]);
  assert.match(nodes["#verdict"].textContent, /^Partly22 of 37 covered/);
});

test("recorded vaccine sequence report is distinct from predicted epitopes", () => {
  const data = result(51, 32);
  data.vaccine_sequence_report = { primary: { recovered: 12, targets: 36 }, summary: {
    n_variants_any_peptide_matched: 15, n_variants: 36,
    n_peptide_entries_matched: 37, n_peptide_entries: 78,
  } };
  const nodes = render(data, ["renderVaccineSequenceSummary"]);
  assert.match(nodes["#vaccine-sequence-summary"].textContent, /12\/36 targets recovered by the top RNA-supported sequence/);
  assert.doesNotMatch(nodes["#vaccine-sequence-summary"].textContent, /37\/78|15\/36/);
  assert.equal(nodes["#vaccine-sequence-summary"].hidden, false);
});

test("a pipeline that has never run stays pending", () => {
  const data = result(0, 0);
  data.has_exacto_run = false;
  const nodes = render(data, ["renderVerdict"]);
  assert.match(nodes["#verdict"].textContent, /^Not yet run/);
});

test("input errors are visible even when no external step ran", () => {
  const data = result(0, 0);
  data.runs.push({ sample: "T1-ONT", arm: "reads", status: "failed",
    error: "cannot encode MAP2", steps: [] });
  const nodes = render(data, ["renderObservedFailures"]);
  assert.match(nodes["#observed-list"].textContent, /pipelinefailed/);
  assert.match(nodes["#observed-list"].textContent, /cannot encode MAP2/);
});

async function liveRun(jobs, { failPage, status = "completed", total = jobs.length } = {}) {
  const requests = [];
  const timers = [];
  const run = {
    id: 49, run_number: 49, status, conclusion: "success",
    html_url: "https://github.com/example/benchmark/actions/runs/49",
    run_started_at: new Date(Date.now() - 7 * 60000).toISOString(),
    updated_at: new Date(Date.now() - 60000).toISOString(),
    head_branch: "main", head_sha: "811778362708e65e01716f954b26124b0ff3b0e5",
  };
  const { nodes, context } = renderer({ repo: "example/benchmark", workflow_file: "exacto-test.yml" }, {
    fetch: async (url) => {
      requests.push(url);
      const { pathname, searchParams } = new URL(url);
      if (pathname.endsWith("/runs")) {
        return { ok: true, json: async () => ({ workflow_runs: [run] }) };
      }
      assert.match(pathname, /\/runs\/49\/jobs$/);
      const page = Number(searchParams.get("page") || 1);
      const size = Number(searchParams.get("per_page"));
      if (page === failPage) return { ok: false, status: 503 };
      return { ok: true, json: async () => ({
        total_count: total, jobs: jobs.slice((page - 1) * size, page * size),
      }) };
    },
    setTimeout: (callback, delay) => timers.push({ callback, delay }),
  });
  await vm.runInContext("renderLiveRun()", context);
  return { node: nodes["#live-run"], requests, timers };
}

function workflowJobs(count) {
  return Array.from({ length: count }, (_, index) => ({
    name: `job-${index + 1}`, status: "completed", conclusion: "success", steps: [],
    html_url: `https://github.com/example/benchmark/actions/runs/49/job/${index + 1}`,
  }));
}

test("live status includes all 36 preparation, method and publishing jobs", async () => {
  const jobs = workflowJobs(36);
  const { node, requests } = await liveRun(jobs);
  const cards = node.children.find((child) => child.className === "live-jobs").children;
  assert.equal(cards.length, 36);
  assert.equal(cards.at(-1).children[0].textContent, "job-36");
  assert.equal(requests.length, 2); // One workflow request and one jobs request.
});

test("live status follows every jobs page for a larger matrix", async () => {
  const { node, requests } = await liveRun(workflowJobs(205));
  const cards = node.children.find((child) => child.className === "live-jobs").children;
  assert.equal(cards.length, 205);
  assert.equal(new Set(cards.map((card) => card.children[0].textContent)).size, 205);
  assert.equal(cards.at(-1).children[0].textContent, "job-205");
  assert.equal(requests.length, 4);
});

test("preparation explains the waiting methods and uses the run start time", async () => {
  const jobs = workflowJobs(5).map((job, index) => ({ ...job,
    name: `prepare (sample-${index})`, status: "in_progress", conclusion: null,
  }));
  const { node, timers } = await liveRun(jobs, { status: "in_progress" });
  assert.match(node.textContent, /Method jobs start after the preparation phase finishes/);
  assert.match(node.textContent, /started 7 min ago/);
  assert.equal(timers.length, 1);
});

test("a failed later page shows an error instead of an incomplete job list", async () => {
  const { node, timers } = await liveRun(workflowJobs(101), { failPage: 2, status: "in_progress" });
  assert.match(node.textContent, /Could not load job progress/);
  assert.ok(!node.children.some((child) => child.className === "live-jobs"));
  assert.equal(timers.length, 1);
});

test("an empty later page stops polling that jobs snapshot", async () => {
  const { node, requests } = await liveRun(workflowJobs(100), { total: 101 });
  assert.equal(node.children.find((child) => child.className === "live-jobs").children.length, 100);
  assert.equal(requests.length, 3);
});

test('amino-acid failures are named beside the score and show the actual failed residue', () => {
  const data = result(3, 2, { reads: { outcome: 'proteoform' } });
  data.variants = [
    { gene: 'ABI3BP', protein_change: 'p.Arg472Met', rna_support: { category: 'multiple' }, recovery: {
      outcome: 'proteoform', residue_confirmed: false, expected: { alt_aa: 'M' }, samples: {
        T2: { arms: { reads: { proteoforms: [{ context: 'AAWBB', context_start: 10,
          mutant_residue_indices: [11], mutant_residues: 'W' }] } } },
      },
    } },
    { gene: 'CONTROL', rna_support: { category: 'single' }, recovery: { outcome: 'peptide', residue_confirmed: true } },
    { gene: 'NO_PROTEIN', rna_support: { category: 'zero' }, recovery: { outcome: 'rna_only', residue_confirmed: false } },
  ];
  const nodes = render(data, ['renderVerdict', 'renderEvidenceBreakdown']);
  assert.match(nodes['#verdict'].textContent, /1 amino-acid failures: ABI3BP/);
  assert.doesNotMatch(nodes['#residue-failures'].textContent, /CONTROL|NO_PROTEIN/);
  assert.match(nodes['#residue-failures'].textContent, /Expected M; observed W/);
  const sequence = nodes['#residue-failures'].children[1].children.find(x => x.tag === 'code');
  assert.equal(sequence.children[2].className, 'sequence-vaccine-difference');
  assert.equal(sequence.children[2].textContent, 'W');
  assert.match(nodes['#rna-breakdown'].textContent, /2\+ mutant RNA reads11\/10\/1/);
  assert.match(nodes['#rna-breakdown'].textContent, /Single-read evidence11\/11\/1/);
});

test('RNA support groups do not score failed or missing outputs as failed recovery', () => {
  const data = result(0, 0);
  data.variants[0].rna_support = { category: 'multiple' };
  data.variants[0].recovery.outcome = 'no_reads';
  const nodes = render(data, ['renderEvidenceBreakdown']);
  assert.match(nodes['#rna-breakdown'].textContent, /2\+ mutant RNA reads1——/);
  assert.doesNotMatch(nodes['#rna-breakdown'].textContent, /0\/1/);
});
