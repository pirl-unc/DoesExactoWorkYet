"""Assemble the GitHub Pages site from whatever results exist.

The variant table only needs ``results/vaccine_variants.json``, so the site
builds and deploys even when the Exacto run has not happened yet (or failed) —
it just says so instead of showing recovery columns.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from . import assays, inventory
from .config import ARMS, REPO_ROOT, RESULTS_DIR, SITE_DIR, TIMEPOINTS
from .extract_reads import stats_path
from .sources import configuration, data_sources, reproduction

WEB_DIR = REPO_ROOT / "web"
HISTORY_PATH = RESULTS_DIR / "history.json"


def git(*args: str) -> str | None:
    try:
        return subprocess.run(
            ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def load(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def extraction_summary(exacto_payload: dict | None) -> dict:
    """How many reads went in, per timepoint.

    Prefer the copy carried through the scored results: the site is often built
    in a job that never had ``work/``.
    """
    summary = dict((exacto_payload or {}).get("extraction") or {})
    for timepoint in TIMEPOINTS:
        if timepoint.name in summary:
            continue
        stats = load(stats_path(timepoint))
        if stats:
            summary[timepoint.name] = {
                key: stats[key]
                for key in (
                    "n_reads",
                    "n_spanning_reads",
                    "n_context_reads",
                    "mean_read_length",
                )
                if key in stats
            }
    return summary


def update_history(summary: dict) -> list[dict]:
    """Record whether the answer changed, not that the site was rebuilt.

    The site rebuilds on any push to ``web/``, so appending unconditionally would
    fill the track record with rows that all say the same thing. A run whose
    Exacto version and score match the last entry refreshes that entry in place;
    only a different answer earns a new row.
    """
    history = load(HISTORY_PATH) or []
    entry = {
        "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "exacto_version": summary.get("exacto_version"),
        "commit": summary.get("commit"),
        "n_variants": summary.get("n_variants"),
        "n_testable": summary.get("n_testable"),
        "n_recovered": summary.get("n_recovered"),
    }
    if entry["n_recovered"] is None:
        return history

    def score(item: dict) -> tuple:
        return (
            item.get("exacto_version"),
            item.get("n_testable"),
            item.get("n_recovered"),
        )

    if history and score(history[-1]) == score(entry):
        # Same answer as last time: keep the row, move it to now.
        entry["date"] = history[-1]["date"]
        history[-1] = entry
    else:
        history.append(entry)

    HISTORY_PATH.write_text(json.dumps(history, indent=2) + "\n")
    return history


def build_payload() -> dict:
    variants_payload = load(RESULTS_DIR / "vaccine_variants.json")
    if variants_payload is None:
        raise SystemExit(
            "results/vaccine_variants.json missing — run pipeline.fetch_osteosarc"
        )
    exacto_payload = load(RESULTS_DIR / "exacto_results.json")
    environment = load(RESULTS_DIR / "environment.json") or {}

    by_variant_id = {}
    if exacto_payload:
        by_variant_id = {
            item["variant_id"]: item for item in exacto_payload["variants"]
        }

    assay_columns = assays.columns(variants_payload["variants"])

    variants = []
    for variant in variants_payload["variants"]:
        recovery = by_variant_id.get(variant["variant_id"])
        # assay_support is 200 rows per variant, so it is collapsed into a VAF
        # grid for the page; the per-sample rows stay in
        # results/vaccine_variants.json for anyone who wants them.
        trimmed = {
            key: value for key, value in variant.items() if key != "assay_support"
        }
        variants.append(
            {
                **trimmed,
                "assay_matrix": assays.matrix(variant),
                "germline_matrix": assays.germline(variant),
                "recovery": recovery,
            }
        )

    summary = {
        "n_variants": variants_payload["n_variants"],
        "n_peptide_entries": variants_payload["n_peptide_entries"],
        "exacto_version": environment.get("exacto_version"),
        "commit": git("rev-parse", "--short", "HEAD"),
    }
    if exacto_payload:
        summary.update(
            {
                "n_testable": exacto_payload["n_testable"],
                "n_recovered": exacto_payload["n_recovered"],
                "outcome_counts": exacto_payload["outcome_counts"],
                "recovered_genes": exacto_payload["recovered_genes"],
                "n_residue_checkable": exacto_payload.get("n_residue_checkable"),
                "n_residue_confirmed": exacto_payload.get("n_residue_confirmed"),
                "n_with_vaccine_epitopes": exacto_payload.get("n_with_vaccine_epitopes"),
                "n_epitope_confirmed": exacto_payload.get("n_epitope_confirmed"),
            }
        )

    history = update_history(summary)
    extraction = extraction_summary(exacto_payload)

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        # Lets the page ask the GitHub API what is running right now. The site
        # is static, so live status has to come from the browser.
        "repo": "pirl-unc/DoesExactoWorkYet",
        "workflow_file": "exacto-test.yml",
        "has_exacto_run": exacto_payload is not None,
        "summary": summary,
        "environment": environment,
        "sources": variants_payload["source"],
        "vaccine_names": variants_payload["vaccine_names"],
        "vaccine_set_sizes": variants_payload["vaccine_set_sizes"],
        "timepoints": [
            {
                "name": timepoint.name,
                "label": timepoint.label,
                "biopsy_date": timepoint.biopsy_date,
                "bam_url": timepoint.bam_url,
            }
            for timepoint in TIMEPOINTS
        ],
        "arms": list(ARMS),
        "extraction": extraction,
        "data_sources": data_sources(extraction),
        "assay_columns": assay_columns,
        "inventory": inventory.load(),
        "configuration": configuration(),
        "reproduction": reproduction(),
        "findings": (load(RESULTS_DIR / "findings.json") or {}).get("findings", []),
        "runs": (exacto_payload or {}).get("runs", []),
        "variants": variants,
        "history": history,
    }


def main() -> None:
    payload = build_payload()

    if SITE_DIR.exists():
        shutil.rmtree(SITE_DIR)
    shutil.copytree(WEB_DIR, SITE_DIR)
    (SITE_DIR / "data.json").write_text(json.dumps(payload, indent=2) + "\n")
    # Pages would otherwise run the files through Jekyll.
    (SITE_DIR / ".nojekyll").touch()

    print(f"site -> {SITE_DIR}")
    print(f"  variants: {payload['summary']['n_variants']}")
    if payload["has_exacto_run"]:
        print(
            f"  recovered: {payload['summary']['n_recovered']}"
            f"/{payload['summary']['n_testable']} testable"
        )
    else:
        print("  no Exacto run yet — site shows the variant table only")


if __name__ == "__main__":
    main()
