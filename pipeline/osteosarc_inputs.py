"""Frozen osteosarc catalogue and source identities used by every benchmark job."""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
from email.utils import parsedate_to_datetime
from importlib.metadata import version
from pathlib import Path

INPUTS_DIR = Path(__file__).resolve().parents[1] / "data" / "osteosarc"


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def manifest() -> dict:
    return json.loads((INPUTS_DIR / "manifest.json").read_text())


def input_id() -> str:
    return digest(INPUTS_DIR / "manifest.json")


def frozen_inputs() -> dict:
    metadata = manifest()
    path = INPUTS_DIR / "inputs.json.gz"
    if digest(path) != metadata["inputs_sha256"]:
        raise ValueError("Frozen osteosarc catalogue checksum changed")
    with gzip.open(path, "rt") as handle:
        inputs = json.load(handle)
    if inputs["snapshot"] != metadata["snapshot"]:
        raise ValueError("Frozen osteosarc snapshot identity differs from its manifest")
    return inputs


def source_file(sample: str):
    from osteosarc import File
    from osteosarc.models import SampleClaim

    values = dict(manifest()["sources"][sample])
    values["claims"] = tuple(SampleClaim(**claim) for claim in values["claims"])
    values["index_urls"] = tuple(values["index_urls"])
    return File(**values)


def require_version() -> None:
    expected = manifest()["osteosarc_version"]
    if version("osteosarc") != expected:
        raise ValueError(f"This input set requires osteosarc=={expected}")


def regional_reads(sample: str, regions: list[dict], cache_path: Path):
    """Verified, cached whole-gene extracts, retaining original BAM records."""
    from osteosarc import Cache, ReadFilter, Region, extract_reads, inspect_alignment

    require_version()
    source = source_file(sample)
    cache = Cache(cache_path)
    snapshot_id = manifest()["snapshot"]["id"]
    info = inspect_alignment(source, cache=cache, snapshot_id=snapshot_id)
    check_source_identity(source, info.receipt["remote_identity"])
    lengths = {row["SN"]: row["LN"] for row in info.header["SQ"]}
    targets = [
        Region(r["chrom"], r["start"] - 1, min(r["end"], lengths[r["chrom"]]), "GRCh38")
        for r in regions
    ]
    try:
        subset = extract_reads(
            source,
            targets,
            cache=cache,
            snapshot_id=snapshot_id,
            filters=ReadFilter(exclude_flags=0x900),
            fetch_pairs=False,
            timeout=1800,
        )
    except subprocess.CalledProcessError as error:
        # osteosarc captures stderr, but CalledProcessError.__str__ omits it.
        # Preserve the diagnostic in CI rather than guessing why I/O failed.
        stderr = error.stderr or "No stderr was returned by the extraction command."
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        raise RuntimeError(
            f"{sample}: read extraction failed (exit {error.returncode}).\n{stderr}"
        ) from error
    check_source_identity(source, subset.receipt["remote_identity"])
    return subset


def check_source_identity(source, identity: dict) -> None:
    """Match the recorded inventory, as well as checking stability during I/O."""
    if (
        source.size is not None
        and int(identity.get("content-length") or -1) != source.size
    ):
        raise ValueError(
            f"{source.key}: remote size differs from the frozen osteosarc inventory"
        )
    if isinstance(source.modified, (int, float)):
        modified = identity.get("last-modified")
        if (
            not modified
            or abs(parsedate_to_datetime(modified).timestamp() - source.modified) > 1
        ):
            raise ValueError(
                f"{source.key}: remote modification date differs from the frozen inventory"
            )
