"""Reproducible catalogue, read acquisition and preparation boundaries."""

from __future__ import annotations

import copy
import gzip
import json
from types import SimpleNamespace

import pytest

from pipeline import (
    build_site,
    evaluate,
    extract_reads,
    fetch_osteosarc,
    osteosarc_inputs,
    prepare_inputs,
)
from pipeline.run_exacto import as_graph_operation


def test_all_51_memberships_and_corrected_alleles_are_retained():
    payload = fetch_osteosarc.build_variant_records()
    assert payload["n_variants"] == payload["n_ready"] == 51
    by_gene = {v["gene"]: v for v in payload["variants"]}
    assert by_gene["MAP2"]["ref"] == "CCTGGGCTACTGTGTGTTCAATAAGTACACAGT"
    assert by_gene["MAP2"]["alt"] == "CAGGG"
    assert "allele-MAP2-chr2-209694768" in by_gene["MAP2"]["corrections"]
    assert by_gene["TRMO"]["gene"] == "TRMO"
    assert (
        len([v for v in payload["variants"] if v["vaccine_membership"]["overlap_rows"]])
        == 37
    )
    assert (
        "ATRX" in by_gene and not by_gene["ATRX"]["vaccine_membership"]["overlap_rows"]
    )
    for v in payload["variants"]:
        as_graph_operation(v)


def test_unresolved_membership_stays_visible():
    inputs = copy.deepcopy(osteosarc_inputs.frozen_inputs())
    extra = next(v for v in inputs["variants"] if v["gene"] == "ATRX")
    extra["status"] = "unresolved"
    extra["alleles"] = []
    payload = fetch_osteosarc.build_variant_records(inputs)
    assert payload["n_variants"] == 51 and payload["n_ready"] == 50
    extra = next(v for v in payload["variants"] if v["gene"] == "ATRX")
    assert extra["allele_status"] == "unresolved" and extra["ref"] is None


def test_modified_frozen_catalogue_is_rejected(tmp_path, monkeypatch):
    manifest = osteosarc_inputs.manifest()
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    (tmp_path / "inputs.json.gz").write_bytes(gzip.compress(b"{}"))
    monkeypatch.setattr(osteosarc_inputs, "INPUTS_DIR", tmp_path)
    with pytest.raises(ValueError, match="checksum"):
        osteosarc_inputs.frozen_inputs()


def test_source_inventory_mismatch_is_rejected():
    source = SimpleNamespace(key="example.bam", size=123, modified=0)
    with pytest.raises(ValueError, match="size"):
        osteosarc_inputs.check_source_identity(source, {"content-length": "124"})
    with pytest.raises(ValueError, match="modification"):
        osteosarc_inputs.check_source_identity(source, {"content-length": "123"})
    osteosarc_inputs.check_source_identity(
        source,
        {
            "content-length": "123",
            "last-modified": "Thu, 01 Jan 1970 00:00:00 GMT",
        },
    )


def test_region_padding_is_clipped_and_the_snapshot_is_passed_to_osteosarc(
    tmp_path, monkeypatch
):
    import osteosarc

    source = SimpleNamespace(key="example.bam", size=123, modified=None)
    identity = {"content-length": "123"}
    monkeypatch.setattr(osteosarc_inputs, "source_file", lambda name: source)
    monkeypatch.setattr(
        osteosarc,
        "inspect_alignment",
        lambda *a, **k: SimpleNamespace(
            header={"SQ": [{"SN": "chrM", "LN": 16569}]},
            receipt={"remote_identity": identity},
        ),
    )
    calls = []

    def acquire(file, regions, **kwargs):
        calls.append((file, regions, kwargs))
        return SimpleNamespace(receipt={"remote_identity": identity})

    monkeypatch.setattr(osteosarc, "extract_reads", acquire)
    osteosarc_inputs.regional_reads(
        "T1-ONT", [{"chrom": "chrM", "start": 1671, "end": 24148}], tmp_path
    )
    _, regions, kwargs = calls[0]
    assert (regions[0].start, regions[0].end) == (1670, 16569)
    assert kwargs["snapshot_id"] == osteosarc_inputs.manifest()["snapshot"]["id"]
    assert kwargs["filters"].exclude_flags == 0x900
    assert kwargs["fetch_pairs"] is False


def test_gzip_fastq_rebuilds_have_identical_bytes(tmp_path):
    a, b = tmp_path / "a.fastq.gz", tmp_path / "b.fastq.gz"
    for p in (a, b):
        with extract_reads.fastq_writer(p) as out:
            out.write("@read\nACGT\n+\nIIII\n")
    assert a.read_bytes() == b.read_bytes()


def test_prepared_inputs_reject_changed_files_and_recipes(tmp_path, monkeypatch):
    monkeypatch.setattr(prepare_inputs, "WORK_DIR", tmp_path)
    f = tmp_path / "reads.fastq.gz"
    f.write_bytes(b"original")
    monkeypatch.setattr(prepare_inputs, "prepared_files", lambda sample: [f])
    p = prepare_inputs.receipt_path("T1-ONT")
    p.parent.mkdir()
    p.write_text(
        json.dumps(
            {
                "preparation_id": prepare_inputs.preparation_id(),
                "files": {"reads.fastq.gz": osteosarc_inputs.digest(f)},
            }
        )
    )
    assert prepare_inputs.verified("T1-ONT")
    f.write_bytes(b"changed")
    assert not prepare_inputs.verified("T1-ONT")
    with pytest.raises(ValueError, match="prepared inputs"):
        prepare_inputs.prepare(["T1-ONT"], verify_only=True)


def test_mixed_catalogues_cannot_be_merged(tmp_path, monkeypatch):
    monkeypatch.setattr(evaluate, "SCORED_DIR", tmp_path)
    p = tmp_path / "score.json"
    p.write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "sample": "T1-ONT",
                        "arm": "reads",
                        "status": "ok",
                        "input_id": "other",
                        "variants": {},
                    }
                ]
            }
        )
    )
    with pytest.raises(ValueError, match="different osteosarc"):
        evaluate.merge()


def test_old_37_locus_result_is_not_attached_to_the_new_panel(tmp_path, monkeypatch):
    # Existing fixture exercises the rest of the site builder. Its published
    # result is intentionally from an earlier, incompatible panel.
    source = fetch_osteosarc.build_variant_records()
    (tmp_path / "vaccine_variants.json").write_text(json.dumps(source))
    (tmp_path / "exacto_results.json").write_text(json.dumps({"input_id": "old"}))
    monkeypatch.setattr(build_site, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(build_site, "HISTORY_PATH", tmp_path / "history.json")
    payload = build_site.build_payload()
    assert payload["summary"]["n_variants"] == 51
    assert payload["has_exacto_run"] is False
    assert not any(v["recovery"] for v in payload["variants"])
