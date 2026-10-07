"""Reproducible catalogue, read acquisition and preparation boundaries."""

from __future__ import annotations

import copy
import gzip
import json
import subprocess
from types import SimpleNamespace

import pytest

from pipeline import (
    build_reference,
    build_site,
    evaluate,
    extract_reads,
    fetch_osteosarc,
    osteosarc_inputs,
    prepare_inputs,
    run_exacto,
)
from pipeline.config import SAMPLES_BY_NAME
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


@pytest.mark.parametrize("gene", ["MAP2", "TECPR1"])
def test_unresolved_overlap_members_keep_assertions_and_skip_exacto(
    gene, tmp_path, monkeypatch
):
    inputs = copy.deepcopy(osteosarc_inputs.frozen_inputs())
    baseline = fetch_osteosarc.build_variant_records(inputs)
    original = next(v for v in baseline["variants"] if v["gene"] == gene)
    variant = next(v for v in inputs["variants"] if v["gene"] == gene)
    variant["status"] = "unresolved"
    variant["alleles"] = []
    payload = fetch_osteosarc.build_variant_records(inputs)
    assert payload["n_variants"] == 51 and payload["n_ready"] == 50
    assert payload["n_peptide_entries"] == baseline["n_peptide_entries"] == 38
    unresolved = next(v for v in payload["variants"] if v["gene"] == gene)
    assert unresolved["ref"] is None and unresolved["alt"] is None
    assert unresolved["vaccines"] == original["vaccines"]
    assert unresolved["elispot"] == original["elispot"]
    assert unresolved["peptide_classes"] == original["peptide_classes"]
    assert unresolved["vaccine_membership"] == {
        **original["vaccine_membership"],
        "overlap_join_status": "unresolved_allele",
    }
    # Removing rows from the strict join must not shift other members' indices.
    assert [v for v in payload["variants"] if v["gene"] != gene] == [
        v for v in baseline["variants"] if v["gene"] != gene
    ]
    (tmp_path / "vaccine_variants.json").write_text(json.dumps(payload))
    monkeypatch.setattr(build_reference, "RESULTS_DIR", tmp_path)
    assert len(build_reference.load_variants()) == 50
    assert all(v["gene"] != gene for v in build_reference.load_variants())
    assert len(build_reference.load_variants(include_unresolved=True)) == 51


def test_overlap_assertions_alone_include_an_unresolved_member():
    inputs = copy.deepcopy(osteosarc_inputs.frozen_inputs())
    variant = next(v for v in inputs["variants"] if v["gene"] == "MAP2")
    variant.update(status="unresolved", alleles=[], vaccines=[], vaccine_count=0)
    variant["annotations"]["source_vaccines"] = []
    payload = fetch_osteosarc.build_variant_records(inputs)
    member = next(v for v in payload["variants"] if v["gene"] == "MAP2")
    assert payload["n_variants"] == 51 and payload["n_ready"] == 50
    assert member["vaccine_membership"]["included"]
    assert member["vaccines"] == ["JLF V1", "JLF V2", "JLF V3", "mRNA"]


@pytest.mark.parametrize(
    "problem", ["unknown", "ambiguous_ready", "ambiguous_unresolved"]
)
def test_overlap_join_still_rejects_unknown_or_ambiguous_loci(problem):
    from osteosarc.errors import IntegrityError

    inputs = copy.deepcopy(osteosarc_inputs.frozen_inputs())
    variant = next(v for v in inputs["variants"] if v["gene"] == "MAP2")
    if problem == "unknown":
        row = next(r for r in inputs["vaccine_rows"] if r["gene"] == "MAP2")
        row["pos"] = "1"
    else:
        if problem == "ambiguous_unresolved":
            variant.update(status="unresolved", alleles=[])
        inputs["variants"].append({**copy.deepcopy(variant), "id": "other-MAP2"})
    with pytest.raises(IntegrityError, match="joins"):
        fetch_osteosarc.build_variant_records(inputs)


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


def test_regional_extraction_surfaces_captured_samtools_stderr(tmp_path, monkeypatch):
    import osteosarc

    source = SimpleNamespace(key="example.bam", size=123, modified=None)
    monkeypatch.setattr(osteosarc_inputs, "source_file", lambda name: source)
    monkeypatch.setattr(osteosarc, "inspect_alignment", lambda *a, **k: SimpleNamespace(
        header={"SQ": [{"SN": "chr1", "LN": 1000}]},
        receipt={"remote_identity": {"content-length": "123"}},
    ))

    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, ["samtools", "view"], stderr=b"HTTP connection reset")

    monkeypatch.setattr(osteosarc, "extract_reads", fail)
    with pytest.raises(RuntimeError, match="T2-ONT: read extraction failed.*\\nHTTP connection reset"):
        osteosarc_inputs.regional_reads("T2-ONT", [{"chrom": "chr1", "start": 10, "end": 20}], tmp_path)


@pytest.mark.parametrize(
    "changed_input",
    [
        None,
        "spanning_fastq",
        "reads_arm_fastq",
        "context_fastq",
        "context_reads_per_region_cap",
        "spanning_reads_per_variant_cap",
        "reads_arm_reads_per_variant_cap",
        "synthetic_base_quality",
        "regions",
    ],
)
def test_rebuilt_reads_merge_only_when_inputs_and_settings_match(
    tmp_path, monkeypatch, changed_input
):
    sample = SAMPLES_BY_NAME["T1-ONT"]
    monkeypatch.setattr(run_exacto, "EXACTO_DIR", tmp_path / "exacto")
    monkeypatch.setattr(evaluate, "SCORED_DIR", tmp_path / "scored")
    monkeypatch.setattr(evaluate, "RESULTS_DIR", tmp_path / "results")
    monkeypatch.setattr(evaluate, "load_variants", lambda **kwargs: [])

    def stop_before_tools(*args, **kwargs):
        raise run_exacto.StepFailed("stopped after input validation")

    monkeypatch.setattr(run_exacto, "align", stop_before_tools)
    monkeypatch.setattr(run_exacto.Runner, "run", stop_before_tools)
    stats = {
        "sample": sample.name,
        "input_id": osteosarc_inputs.input_id(),
        "n_spanning_reads": 1,
        "context_reads_per_region_cap": 5000,
        "spanning_reads_per_variant_cap": 3000,
        "reads_arm_reads_per_variant_cap": 600,
        "synthetic_base_quality": 30,
        "regions": [{"chrom": "chr1", "start": 1, "end": 100, "genes": ["TEST"]}],
    }
    runs = []
    for index, arm in enumerate(("reads", "assembly")):
        # A fresh cache/workspace changes receipts and absolute paths even when
        # the source reads and sampling settings produce identical FASTQ bytes.
        reads_dir = tmp_path / f"workspace-{index}" / "reads"
        reads_dir.mkdir(parents=True)
        monkeypatch.setattr(extract_reads, "READS_DIR", reads_dir)
        stats["acquisition"] = {
            "index_receipt": {"retrieved_at": f"2026-10-0{index + 1}T00:00:00Z"},
            "path": str(reads_dir / "regional.bam"),
        }
        for name in ("spanning_fastq", "reads_arm_fastq", "context_fastq"):
            path = getattr(extract_reads, name)(sample)
            stats[name] = str(path)
            sequence = "TGCA" if index and changed_input == name else "ACGT"
            with extract_reads.fastq_writer(path) as out:
                out.write(f"@read\n{sequence}\n+\nIIII\n")
        if index and changed_input == "regions":
            stats["regions"][0]["end"] += 1
        elif index and changed_input and not changed_input.endswith("_fastq"):
            stats[changed_input] += 1
        # Dictionary order is serialization detail, not an input difference.
        outputs = extract_reads.extraction_outputs(sample)
        stats["files"] = {
            path.name: osteosarc_inputs.digest(path)
            for path in (reversed(outputs) if index else outputs)
        }
        extract_reads.stats_path(sample).write_text(json.dumps(stats, indent=index + 1))
        run = run_exacto.run_arm(sample, arm, [], 1)
        assert run["error"] == "stopped after input validation"
        assert run["reads_input_id"]
        runs.append(run)

    # Exercise the identity that is actually saved and passed through scoring,
    # with external bioinformatics tools stopped just after input validation.
    scored = evaluate.score_samples([sample.name])
    assert {run["reads_input_id"] for run in scored} == {
        run["reads_input_id"] for run in runs
    }
    if changed_input is None:
        assert runs[0]["reads_input_id"] == runs[1]["reads_input_id"]
        assert len(evaluate.merge()["runs"]) == 2
    else:
        with pytest.raises(ValueError, match="Cannot merge different prepared reads"):
            evaluate.merge()


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


def test_incremental_preparation_preserves_earlier_sample_receipts(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(prepare_inputs, "WORK_DIR", tmp_path)
    monkeypatch.setattr(fetch_osteosarc, "RESULTS_DIR", tmp_path / "results")
    monkeypatch.setattr(build_reference, "RESULTS_DIR", tmp_path / "results")
    monkeypatch.setattr(extract_reads, "READS_DIR", tmp_path / "reads")
    reference = tmp_path / "reference"
    reference.mkdir()
    for name, filename in (
        ("MASKED_FASTA", "genome.fa"),
        ("SUBSET_GTF", "genes.gtf.gz"),
        ("GENE_PROTEINS", "proteins.fa"),
        ("REGIONS_JSON", "regions.json"),
    ):
        monkeypatch.setattr(build_reference, name, reference / filename)
    monkeypatch.setattr(extract_reads, "REGIONS_JSON", build_reference.REGIONS_JSON)
    gtf = tmp_path / "source.gtf.gz"
    content = '# header\nchr1\ttest\tgene\t10\t20\t.\t+\t.\tgene_id "TEST";\n'
    gtf.write_bytes(gzip.compress(content.encode()))
    clock = {"time": 1000}
    monkeypatch.setattr(gzip.time, "time", lambda: clock["time"])

    def build():
        build_reference.MASKED_FASTA.write_text(">chr1\nACGT\n")
        build_reference.MASKED_FASTA.with_suffix(".fa.fai").write_text(
            "chr1\t4\t6\t4\t5\n"
        )
        build_reference.GENE_PROTEINS.write_text(">TEST\nM\n")
        build_reference.REGIONS_JSON.write_text("[]")
        build_reference.write_subset_gtf(
            [build_reference.Region("chr1", 1, 100, ("TEST",))],
            gtf,
            build_reference.SUBSET_GTF,
        )

    def extract(sample, regions, variants):
        for path in extract_reads.extraction_outputs(sample):
            path.parent.mkdir(parents=True, exist_ok=True)
            with extract_reads.fastq_writer(path) as out:
                out.write("@read\nACGT\n+\nIIII\n")
        extract_reads.stats_path(sample).write_text("{}")

    monkeypatch.setattr(build_reference, "main", build)
    monkeypatch.setattr(extract_reads, "extract", extract)
    prepare_inputs.prepare(["T1-ONT"])
    earlier_receipt = prepare_inputs.receipt_path("T1-ONT").read_bytes()
    reference_bytes = build_reference.SUBSET_GTF.read_bytes()
    clock["time"] = 2000
    prepare_inputs.prepare(["T2-ONT"])
    assert build_reference.SUBSET_GTF.read_bytes() == reference_bytes
    assert gzip.decompress(reference_bytes).decode() == content
    assert prepare_inputs.receipt_path("T1-ONT").read_bytes() == earlier_receipt
    prepare_inputs.prepare(["T1-ONT", "T2-ONT"], verify_only=True)


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
