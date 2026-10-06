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

function render(data, functions) {
  const nodes = {};
  const context = vm.createContext({
    document: {
      querySelector(selector) { return nodes[selector] ||= new Element("div"); },
      createElement(tag) { return new Element(tag); },
    },
    data,
  });
  const source = fs.readFileSync(path.join(__dirname, "../web/app.js"), "utf8");
  // Omit only the async page bootstrap; run the actual helpers and renderers.
  vm.runInContext(source.slice(0, source.lastIndexOf("\nmain().catch")), context);
  vm.runInContext(`DATA = data; ${functions.map((name) => `${name}();`).join(" ")}`, context);
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
