// Use exactly the same window grouping/ranking as the report, without a DOM.
const fs = require("node:fs");
const { primaryWindowReport } = require("../web/peptides.js");
const report = JSON.parse(fs.readFileSync(0, "utf8"));
process.stdout.write(JSON.stringify(primaryWindowReport(report)));
