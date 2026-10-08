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

from . import assays, inventory, paths
from .config import (
    ARMS,
    REPO_ROOT,
    RESULTS_DIR,
    SAMPLES,
    SITE_DIR,
    arms_for,
    minimap2_preset,
)
from .extract_reads import stats_path
from .sources import configuration, data_sources, reproduction
from .vaccine_peptides import results_fingerprint

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


def current_peptide_report(catalogue: dict, results: dict | None) -> dict | None:
    report = load(RESULTS_DIR / "vaccine_peptide_analysis.json")
    if (not report or not results
            or report.get("source", {}).get("input_id") != catalogue.get("source", {}).get("input_id")
            or report.get("analysis", {}).get("results_sha256") != results_fingerprint(results)):
        return None
    return report


def source_rna_evidence(variant: dict) -> list[dict]:
    samples = {sample.bam_url.rsplit("/", 1)[-1]: sample.name for sample in SAMPLES}
    return [
        {**entry, "benchmark_sample": samples.get((entry.get("bam_file") or "").rsplit("/", 1)[-1])}
        for entry in variant.get("assay_support", [])
        if entry.get("tissue") == "tumor"
        and assays.ASSAY_META.get(entry.get("assay_type"), {}).get("kind") == "rna"
    ]


def rna_support_tier(evidence: list[dict]) -> dict:
    """Maximum observed mutant count, never pooled across libraries or methods."""
    def classify(entries):
        counts = [entry["alt_reads"] for entry in entries
                  if isinstance(entry.get("alt_reads"), int) and entry["alt_reads"] >= 0]
        maximum = max(counts, default=None)
        return {"max_alt_reads": maximum, "category": "unknown" if maximum is None else
                "multiple" if maximum >= 2 else "single" if maximum == 1 else "zero"}

    by_sample = {sample.name: classify([entry for entry in evidence
                                      if entry.get("benchmark_sample") == sample.name])
                 for sample in SAMPLES}
    return {
        **classify([entry for entry in evidence if entry.get("benchmark_sample")]),
        "by_sample": by_sample,
        "missing_samples": [name for name, counts in by_sample.items() if counts["category"] == "unknown"],
        "basis": "Sid mutant reads in benchmark-matched tumor RNA libraries; maximum, not sum",
    }


def peptide_report_with_rna_evidence(report: dict, catalogue: dict, results: dict | None = None) -> dict:
    """Attach source counts without pooling libraries or inferring missing counts.

    These are allele counts at the locus, not vaccine peptide observations.
    Match benchmark samples by the actual BAM basename, not timepoint alone.
    """
    by_id = {variant["variant_id"]: variant for variant in catalogue["variants"]}
    report = json.loads(json.dumps(report))
    scored = {variant["variant_id"]: variant for variant in (results or {}).get("variants", [])}
    for variant in report["variants"]:
        if "exacto_runs" not in variant:
            # Older analyses contain matching proteins only. Retain their run
            # outcomes without mislabelling an unarchived mismatch as no protein.
            variant["exacto_runs"] = [
                {"sample": sample, "arm": arm, "n_proteoforms": entry.get("n_proteoforms", 0),
                 "n_variant_rnas": len({call["transcript_model_id"] for call in entry.get("rna_variant_calls", [])})}
                for sample, values in scored.get(variant["variant_id"], {}).get("samples", {}).items()
                for arm, entry in values.get("arms", {}).items()
            ]
        variant["source_rna_support"] = source_rna_evidence(by_id[variant["variant_id"]])
        variant["rna_support"] = rna_support_tier(variant["source_rna_support"])
    return report


# Runs are keyed on the sequencing sample now, not the biopsy — but a run that
# was already in flight when that changed wrote "T1" where it would now write
# "T1-ONT". Those results are still perfectly good; only their labels are stale.
# ONT was the only platform a timepoint-keyed run could have used, so the
# mapping is exact rather than a guess.
LEGACY_SAMPLE_NAMES = {
    sample.timepoint: sample.name for sample in SAMPLES if sample.platform == "ONT"
}


def migrate_name(name: str) -> str:
    return LEGACY_SAMPLE_NAMES.get(name, name)


def migrate_keys(mapping: dict) -> dict:
    return {migrate_name(key): value for key, value in mapping.items()}


def migrate_payload(payload: dict | None) -> dict | None:
    """Rewrite a timepoint-keyed Exacto result in place as a sample-keyed one."""
    if not payload:
        return payload
    for run in payload.get("runs", []):
        if "sample" not in run and "timepoint" in run:
            run["sample"] = migrate_name(run["timepoint"])
            run.setdefault("platform", "ONT")
            run.setdefault("label", run["sample"])
    for variant in payload.get("variants", []):
        if "samples" not in variant and "timepoints" in variant:
            variant["samples"] = migrate_keys(variant.pop("timepoints"))
    if "extraction" in payload:
        payload["extraction"] = migrate_keys(payload["extraction"])
    return payload


def extraction_summary(exacto_payload: dict | None) -> dict:
    """How many reads went in, per sample.

    Prefer the copy carried through the scored results: the site is often built
    in a job that never had ``work/``.
    """
    summary = migrate_keys(dict((exacto_payload or {}).get("extraction") or {}))
    for sample in SAMPLES:
        if sample.name in summary:
            continue
        stats = load(stats_path(sample))
        if stats:
            summary[sample.name] = {
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
        "input_id": summary.get("input_id"),
        "n_variants": summary.get("n_variants"),
        "n_testable": summary.get("n_testable"),
        "n_recovered": summary.get("n_recovered"),
    }
    # A failed evaluation or a completed run with no covered loci has no
    # recovery fraction. Keep the last measured answer in the track record.
    if entry["n_recovered"] is None or not entry["n_testable"]:
        return history

    def score(item: dict) -> tuple:
        return (
            item.get("exacto_version"),
            item.get("input_id"),
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


# Enough to cover a few weeks of work without turning data.json into an archive.
WORKLOG_LIMIT = 30

# git log's own record/field separators, so a commit body containing newlines,
# tabs or pipes cannot break the parse.
_RECORD, _FIELD = "\x1e", "\x1f"


def worklog(limit: int = WORKLOG_LIMIT) -> list[dict]:
    """What has been changed in this harness, and why, from the commit log.

    The commit messages already carry the diagnosis — what failed, what the
    evidence was, what was done about it. Publishing them means the site says
    why it reports what it reports, rather than presenting a number as though
    it fell out of the sky. Read from git rather than hand-maintained so it
    cannot drift from what actually shipped.

    Distinct from findings.json, which is about bugs in *Exacto*; this is about
    changes to the harness testing it.
    """
    raw = git(
        "log",
        f"-{limit}",
        "--date=short",
        f"--pretty=format:%h{_FIELD}%ad{_FIELD}%an{_FIELD}%s{_FIELD}%b{_RECORD}",
    )
    if not raw:
        return []

    entries = []
    for record in raw.split(_RECORD):
        record = record.strip("\n")
        if not record.strip():
            continue
        parts = record.split(_FIELD)
        if len(parts) < 4:
            continue
        sha, date, author, subject = parts[:4]
        body = parts[4] if len(parts) > 4 else ""
        entries.append(
            {
                "sha": sha,
                "date": date,
                "author": author,
                "subject": subject,
                # Blank-line-separated paragraphs, rewrapped in the browser.
                "body": [
                    " ".join(block.split())
                    for block in body.strip().split("\n\n")
                    if block.strip() and not block.startswith("Claude-Session:")
                ],
                # Results commits are written by CI, not by a person changing
                # how the test works; the site separates the two.
                "kind": "results" if subject.startswith("results:") else "change",
            }
        )
    return entries


def build_payload() -> dict:
    variants_payload = load(RESULTS_DIR / "vaccine_variants.json")
    if variants_payload is None:
        raise SystemExit(
            "results/vaccine_variants.json missing — run pipeline.fetch_osteosarc"
        )
    raw_results = load(RESULTS_DIR / "exacto_results.json")
    peptide_report = current_peptide_report(variants_payload, raw_results)
    exacto_payload = migrate_payload(raw_results)
    source = variants_payload.get("source", {})
    mismatched_result = bool(source.get("input_id") and exacto_payload and (
        exacto_payload.get("input_id") != source["input_id"]
    ))
    if mismatched_result:
        # A previous 37-locus evaluation is history, not a score for the new
        # 51-locus corrected panel. Never attach it by a coincidentally equal ID.
        exacto_payload = None
    environment = {} if mismatched_result else (load(RESULTS_DIR / "environment.json") or {})

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
                "rna_support": rna_support_tier(source_rna_evidence(variant)),
            }
        )

    summary = {
        "input_id": source.get("input_id"),
        "osteosarc_snapshot": source.get("snapshot"),
        "n_variants": variants_payload["n_variants"],
        "n_peptide_entries": variants_payload["n_peptide_entries"],
        "exacto_version": environment.get("exacto_version"),
        "commit": git("rev-parse", "--short", "HEAD"),
    }
    if exacto_payload:
        # Older payloads reported a failed-only merge as 0/0 with every locus
        # labelled no_reads. Infer evaluation status from actual mutation
        # scores so rebuilding the current site corrects that display too.
        n_evaluated = exacto_payload.get("n_evaluated")
        if n_evaluated is None:
            n_evaluated = sum(
                any(entry.get("arms") for entry in item.get("samples", {}).values())
                for item in exacto_payload["variants"]
            )
        evaluation_status = exacto_payload.get("evaluation_status") or (
            "failed" if not n_evaluated
            else "no_coverage" if not exacto_payload["n_testable"] else "ok"
        )
        outcome_counts = exacto_payload["outcome_counts"]
        if not n_evaluated:
            outcome_counts = {
                **dict.fromkeys(outcome_counts, 0),
                "not_run": exacto_payload["n_variants"],
            }
        summary.update(
            {
                "n_evaluated": n_evaluated,
                "evaluation_status": evaluation_status,
                "n_testable": exacto_payload["n_testable"],
                "n_recovered": exacto_payload["n_recovered"],
                "outcome_counts": outcome_counts,
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
        "samples": [
            {
                "name": sample.name,
                "timepoint": sample.timepoint,
                "platform": sample.platform,
                "read_type": sample.read_type,
                "assay": sample.assay,
                "label": sample.label,
                "biopsy_date": sample.biopsy_date,
                "library": sample.library,
                "biosample": sample.biosample,
                "bam_url": sample.bam_url,
                "portal_genotyped": sample.portal_genotyped,
                "provenance": sample.provenance,
                # Only the arms this read type can actually run: asking for an
                # Illumina "reads" preset is a question with no answer, and
                # minimap2_preset is deliberately strict about that.
                "arms": arms_for(sample),
                "minimap2_preset": {
                    arm: minimap2_preset(sample.platform, arm)
                    for arm in arms_for(sample)
                },
            }
            for sample in SAMPLES
        ],
        "arms": list(ARMS),
        "extraction": extraction,
        "data_sources": data_sources(extraction),
        "assay_columns": assay_columns,
        "absent_platforms": assays.absent_platforms(
            variants_payload["variants"]
        ),
        "inventory": inventory.load(),
        "configuration": configuration(),
        "reproduction": reproduction(),
        "findings": (load(RESULTS_DIR / "findings.json") or {}).get("findings", []),
        "worklog": worklog(),
        "paths": paths.analyse(exacto_payload, variants_payload["variants"]),
        "runs": (exacto_payload or {}).get("runs", []),
        "variants": variants,
        "history": history,
        "vaccine_sequence_report": {
            "summary": peptide_report["summary"],
            "analysis": {
                key: value for key, value in peptide_report["analysis"].items()
                if key != "artifacts"
            },
        } if peptide_report else None,
    }


def main() -> None:
    payload = build_payload()

    if SITE_DIR.exists():
        shutil.rmtree(SITE_DIR)
    shutil.copytree(WEB_DIR, SITE_DIR)
    (SITE_DIR / "data.json").write_text(json.dumps(payload, indent=2) + "\n")
    if payload["vaccine_sequence_report"]:
        report = peptide_report_with_rna_evidence(
            load(RESULTS_DIR / "vaccine_peptide_analysis.json"),
            load(RESULTS_DIR / "vaccine_variants.json"),
            load(RESULTS_DIR / "exacto_results.json"),
        )
        (SITE_DIR / "vaccine_peptide_analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    else:
        (SITE_DIR / "vaccine_peptide_analysis.json").write_text(json.dumps({
            "status": "unavailable",
            "message": "A vaccine-sequence analysis is not available for the current benchmark run.",
        }) + "\n")
    for name in ("vaccine_peptide_subset.json", "vaccine_peptide_subset.fasta"):
        if (RESULTS_DIR / name).exists():
            shutil.copyfile(RESULTS_DIR / name, SITE_DIR / name)
    # Pages would otherwise run the files through Jekyll.
    (SITE_DIR / ".nojekyll").touch()

    print(f"site -> {SITE_DIR}")
    print(f"  variants: {payload['summary']['n_variants']}")
    status = payload["summary"].get("evaluation_status")
    if status == "failed":
        print("  evaluation failed: no vaccine mutations scored")
    elif status == "no_coverage":
        print("  no vaccine mutations covered in completed runs")
    elif payload["has_exacto_run"]:
        print(
            f"  recovered: {payload['summary']['n_recovered']}"
            f"/{payload['summary']['n_testable']} testable"
        )
    else:
        print("  no Exacto run yet — site shows the variant table only")


if __name__ == "__main__":
    main()
