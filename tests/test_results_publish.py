"""Run the real publishing shell step against local Git repositories."""

import json
import os
import subprocess
import textwrap
from pathlib import Path

import pytest

WORKFLOW = Path(__file__).parents[1] / ".github/workflows/exacto-test.yml"


def git(cwd, *args):
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


@pytest.mark.parametrize("new_results", [True, False])
def test_results_publication_preserves_code_merged_during_run(tmp_path, new_results):
    remote = tmp_path / "remote.git"
    git(tmp_path, "init", "--bare", str(remote))
    run = tmp_path / "run"
    git(tmp_path, "init", "-b", "main", str(run))
    git(run, "config", "user.name", "Test")
    git(run, "config", "user.email", "test@example.invalid")
    (run / "web").mkdir()
    (run / "results").mkdir()
    (run / "web/report.html").write_text("old report")
    (run / "results/count.json").write_text("1")
    (run / "results/environment.json").write_text(
        json.dumps({"exacto_version": "test"})
    )
    git(run, "add", ".")
    git(run, "commit", "-m", "run starts here")
    git(run, "remote", "add", "origin", str(remote))
    git(run, "push", "origin", "main")

    newer = tmp_path / "newer"
    git(tmp_path, "clone", "--branch", "main", str(remote), str(newer))
    git(newer, "config", "user.name", "Test")
    git(newer, "config", "user.email", "test@example.invalid")
    (newer / "web/report.html").write_text("clarified report")
    (newer / "web/new.js").write_text("new report feature")
    git(newer, "add", ".")
    git(newer, "commit", "-m", "report update during benchmark")
    git(newer, "push", "origin", "main")
    newer_head = git(newer, "rev-parse", "HEAD")

    if new_results:
        (run / "results/count.json").write_text("2")
    script = WORKFLOW.read_text().split(
        "      - name: Commit results\n        run: |\n", 1
    )[1].split("\n      - name:", 1)[0]
    subprocess.run(
        ["bash", "-eu", "-c", textwrap.dedent(script)],
        cwd=run,
        env={**os.environ, "GITHUB_REF_NAME": "main"},
        check=True,
        capture_output=True,
        text=True,
    )

    # The published commit and the checkout used for the next site build must
    # both retain newer website code, including files absent from the old run.
    assert git(remote, "show", "main:web/report.html") == "clarified report"
    assert git(remote, "show", "main:web/new.js") == "new report feature"
    assert (run / "web/report.html").read_text() == "clarified report"
    assert (run / "web/new.js").read_text() == "new report feature"
    assert git(remote, "show", "main:results/count.json") == (
        "2" if new_results else "1"
    )
    assert git(run, "status", "--porcelain") == ""
    if new_results:
        assert git(run, "diff", "--name-only", newer_head, "HEAD") == (
            "results/count.json"
        )
    else:
        assert git(run, "rev-parse", "HEAD") == newer_head
