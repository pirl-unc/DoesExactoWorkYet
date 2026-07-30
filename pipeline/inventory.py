"""What sequencing data osteosarc.com actually holds, per platform and timepoint.

Derived from the portal's own file manifest rather than asserted, because the
obvious question — "isn't there PacBio too, and bulk as well as single-cell?" —
deserves an answer anyone can re-derive rather than take on trust.

The short version, at the time of writing: PacBio exists for T1 only, ONT for
all three biopsies, and every long-read RNA dataset is single-cell. Both
long-read platforms are run; the bulk RNA-seq is Illumina short-read, which
Exacto is not built for.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.request
from collections import defaultdict
from typing import Any

from .config import B2, RESULTS_DIR

MANIFEST_URL = f"{B2}/manifest.txt"
INVENTORY_PATH = RESULTS_DIR / "data_inventory.json"

# Matched on path segments so "cont" cannot look like ONT, nor "long" like a
# platform name.
PLATFORMS = [
    ("PacBio", r"pacbio|pbmm2|isoseq|iso-seq|flnc|hifi|revio|sequel"),
    ("ONT", r"ont|nanopore|minion|promethion|pod5"),
    ("Illumina", r"cellranger|illumina|spaceranger|bostongene|tempus|personalis"),
]
PLATFORM_RE = [
    (name, re.compile(rf"(?:^|[/_.-])(?:{pattern})(?:$|[/_.-])", re.IGNORECASE))
    for name, pattern in PLATFORMS
]
PBI_RE = re.compile(r"\.pbi$", re.IGNORECASE)
TIMEPOINT_RE = re.compile(r"(?:^|[/_.-])(T[0-3])(?:$|[/_.-])")

TIMEPOINTS = ["T0", "T1", "T2", "T3"]


def fetch_manifest() -> list[str] | None:
    """The portal's full file listing. ~30 MB, so failure is not fatal."""
    request = urllib.request.Request(
        MANIFEST_URL,
        headers={"User-Agent": "DoesExactoWorkYet/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
            body = response.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError) as error:
        print(f"  manifest unavailable ({error}); keeping the committed inventory")
        return None
    prefix = f"{B2}/"
    return [
        line.strip().removeprefix(prefix) for line in body.splitlines() if line.strip()
    ]


def build(paths: list[str]) -> dict[str, Any]:
    """Platform × timepoint file counts, plus where each lives."""
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    areas: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))

    for path in paths:
        timepoints = set(TIMEPOINT_RE.findall(path))
        matched = [name for name, regex in PLATFORM_RE if regex.search(path)]
        if PBI_RE.search(path) and "PacBio" not in matched:
            matched.append("PacBio")
        if not matched:
            continue
        area = "/".join(path.split("/")[:2])
        for name in matched:
            for timepoint in timepoints or {"unlabelled"}:
                counts[name][timepoint] += 1
                areas[name][timepoint].add(area)

    grid = []
    for name, _ in PLATFORMS:
        row = {"platform": name, "timepoints": {}}
        for timepoint in [*TIMEPOINTS, "unlabelled"]:
            n = counts[name].get(timepoint, 0)
            if not n:
                continue
            row["timepoints"][timepoint] = {
                "files": n,
                "areas": sorted(areas[name][timepoint])[:4],
            }
        grid.append(row)

    return {
        "manifest_url": MANIFEST_URL,
        "n_files": len(paths),
        "grid": grid,
        "note": (
            "Counts are files whose path names a platform, so they include "
            "intermediates and index files, not just primary data. What matters "
            "for this test: PacBio RNA exists for T1 only, ONT for all three "
            "biopsies, and every long-read RNA dataset on the portal is "
            "single-cell. Both long-read platforms are run — the bulk RNA-seq is "
            "Illumina short-read, which Exacto is not built for."
        ),
    }


def load() -> dict[str, Any] | None:
    return json.loads(INVENTORY_PATH.read_text()) if INVENTORY_PATH.exists() else None


def refresh() -> dict[str, Any] | None:
    """Rebuild from the live manifest, falling back to what is committed."""
    paths = fetch_manifest()
    if paths is None:
        return load()
    inventory = build(paths)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    INVENTORY_PATH.write_text(json.dumps(inventory, indent=2) + "\n")
    return inventory


def main() -> None:
    inventory = refresh()
    if not inventory:
        raise SystemExit("no inventory: the manifest could not be read")
    print(f"{inventory['n_files']:,} files on the portal")
    for row in inventory["grid"]:
        spread = ", ".join(
            f"{tp} ({data['files']:,})" for tp, data in row["timepoints"].items()
        )
        print(f"  {row['platform']:10} {spread}")


if __name__ == "__main__":
    main()
