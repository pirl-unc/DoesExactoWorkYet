"""Prepare each sample once; verify identical inputs before every method runs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import build_reference, extract_reads, fetch_osteosarc
from .config import REPO_ROOT, WORK_DIR, samples_named
from .osteosarc_inputs import digest, input_id


def preparation_id() -> str:
    # Changing a region rule, sampling cap or adapter invalidates the cache.
    hasher = hashlib.sha256(input_id().encode())
    for name in (
        "config",
        "fetch_osteosarc",
        "build_reference",
        "extract_reads",
        "osteosarc_inputs",
        "prepare_inputs",
    ):
        hasher.update((REPO_ROOT / "pipeline" / f"{name}.py").read_bytes())
    hasher.update((REPO_ROOT / "environment.yml").read_bytes())
    return hasher.hexdigest()


def receipt_path(sample: str) -> Path:
    return WORK_DIR / "prepared" / f"{sample}.json"


def prepared_files(sample) -> list[Path]:
    return [
        build_reference.MASKED_FASTA,
        Path(str(build_reference.MASKED_FASTA) + ".fai"),
        build_reference.SUBSET_GTF,
        build_reference.GENE_PROTEINS,
        build_reference.REGIONS_JSON,
        *extract_reads.extraction_outputs(sample),
        extract_reads.stats_path(sample),
    ]


def verified(sample: str) -> bool:
    path = receipt_path(sample)
    if not path.exists():
        return False
    receipt = json.loads(path.read_text())
    if receipt.get("preparation_id") != preparation_id():
        return False
    selected = samples_named([sample])[0]
    expected = {str(path.relative_to(WORK_DIR)) for path in prepared_files(selected)}
    if set(receipt.get("files", {})) != expected:
        return False
    for name, checksum in receipt.get("files", {}).items():
        file = WORK_DIR / name
        if not file.is_file() or digest(file) != checksum:
            return False
    return bool(receipt.get("files"))


def prepare(samples: list[str] | None, *, verify_only: bool = False) -> None:
    selected = samples_named(samples)
    if all(verified(sample.name) for sample in selected):
        print("Prepared osteosarc input checksums verified")
        return
    if verify_only:
        raise ValueError(
            "Missing or changed prepared inputs; run pipeline.prepare_inputs"
        )
    fetch_osteosarc.main()
    build_reference.main()
    regions = extract_reads.load_regions()
    variants = build_reference.load_variants()
    for sample in selected:
        extract_reads.extract(sample, regions, variants)
        files = prepared_files(sample)
        receipt = {
            "input_id": input_id(),
            "preparation_id": preparation_id(),
            "files": {str(path.relative_to(WORK_DIR)): digest(path) for path in files},
        }
        path = receipt_path(sample.name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(receipt, indent=2) + "\n")
    print("Prepared and pinned osteosarc sample inputs")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="*")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    prepare(args.samples, verify_only=args.verify_only)


if __name__ == "__main__":
    main()
