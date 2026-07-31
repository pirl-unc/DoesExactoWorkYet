"""Generate ready-to-file issue bodies from results/findings.json.

Drafts only -- nothing here files anything. Regenerate after editing findings:

    python scripts/build_issues.py

Kept as a script rather than hand-written files so a finding cannot be corrected
in one place and stay wrong in the other.
"""

from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
FINDINGS = ROOT / "results" / "findings.json"
OUT = ROOT / "docs" / "issues"

FOOTER = (
    "Found by [DoesExactoWorkYet](https://github.com/pirl-unc/DoesExactoWorkYet), "
    "an automated end-to-end test of Exacto on open data from "
    "[osteosarc.com](https://osteosarc.com). Full detail, including the failing "
    "command and its stderr: https://pirl-unc.github.io/DoesExactoWorkYet/bugs.html"
)


def body(finding: dict) -> str:
    lines = [
        f"# {finding['title']}",
        "",
        f"**Severity:** {finding['severity']}",
        f"**Observed in:** Exacto {finding['observed_in']}",
        f"**Where:** `{finding['where']}`",
        "",
        "## What happens",
        "",
        finding["detail"],
        "",
    ]
    if finding.get("workaround"):
        lines += ["## Workaround currently in use", "", finding["workaround"], ""]
    if finding.get("suggested_fix"):
        lines += ["## Suggested fix", "", finding["suggested_fix"], ""]
    lines += [
        "## Reproduction and a suggested test",
        "",
        "See [docs/reproductions.md](https://github.com/pirl-unc/DoesExactoWorkYet/"
        "blob/main/docs/reproductions.md) for a minimal reproduction and a "
        "suggested unit test in the crate where this lives.",
        "",
        "---",
        "",
        FOOTER,
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    findings = json.loads(FINDINGS.read_text())["findings"]
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for finding in findings:
        (OUT / f"{finding['id']}.md").write_text(body(finding))
        rows.append((finding["id"], finding["severity"], finding["title"]))

    readme = [
        "# Ready-to-file issue bodies",
        "",
        "Generated from `results/findings.json` by `scripts/build_issues.py`.",
        "Nothing here has been filed \u2014 these are drafts.",
        "",
        "| file | severity | title |",
        "|---|---|---|",
    ]
    readme += [f"| [`{i}.md`]({i}.md) | {s} | {t} |" for i, s, t in rows]
    (OUT / "README.md").write_text("\n".join(readme) + "\n")
    print(f"{len(rows)} issue bodies -> {OUT}")


if __name__ == "__main__":
    main()
