"""Compare the routes a read can take to become something Exacto calls on.

Exacto does not read FASTQs. It reads a BAM of *transcripts*, and how those
transcripts are produced is a choice made before Exacto ever runs. This repo
takes two routes and the difference between them turns out to dominate the
result, so the site states it rather than leaving it to be inferred from two
numbers in a table:

``reads``
    No assembly. Every variant-spanning read is handed over as its own
    transcript. Nothing is filtered on splicing.

``assembly``
    The documented route, and Andy Lee's Nexus subworkflow: RNA-Bloom2 stitches
    reads into full-length contigs, ``nexus_filter_rnabloom2_transcripts``
    drops the poorly supported ones, and ``remove-unspliced-rnas`` drops those
    that do not look spliced.

Two things are measured for each route. The **sequence funnel** — how much of
what went in survives to variant calling — and the **variant ladder**, how many
of the 37 vaccine mutations each route gets to each rung of the verdict. The
second is the question worth asking: a route that keeps more sequence is not
automatically a route that recovers more mutations.

Only one assembler is wired up today. The structure is per-arm rather than
per-assembler so that a second one is a config entry and a row here, not a
rewrite — but the page says "one assembly method" rather than implying a
survey that has not been run.
"""

from __future__ import annotations

from typing import Any

# Rungs of the verdict ladder, in order. Mirrors evaluate.OUTCOME_ORDER.
_RANK = {"no_reads": 0, "no_call": 1, "rna_only": 2, "proteoform": 3, "peptide": 4}

PATH_META = {
    "reads": {
        "label": "Reads",
        "subtitle": "no assembly",
        "assembler": None,
        "description": (
            "Each variant-spanning read is handed to Exacto as its own "
            "transcript. Nothing is assembled, nothing is filtered on splicing. "
            "This began as a control — a way to tell an Exacto miss from an "
            "assembler miss — and is not what the Exacto docs prescribe."
        ),
    },
    "assembly": {
        "label": "Assembly",
        "subtitle": "RNA-Bloom2 + Nexus filter",
        "assembler": "RNA-Bloom2 2.0.1",
        "description": (
            "The documented route: RNA-Bloom2 assembles the reads into "
            "full-length contigs, nexus_filter_rnabloom2_transcripts drops "
            "those without enough read support, and remove-unspliced-rnas drops "
            "what does not look spliced. This is the canonical pipeline."
        ),
    },
}

PATH_ORDER = ["reads", "assembly"]


def _stages(arm: str, runs: list[dict], extraction: dict) -> list[dict]:
    """The sequence funnel for one route, summed over samples.

    Summed rather than averaged: these are absolute quantities of sequence, and
    a mean across samples of very different depth would be a number that
    describes none of them.
    """
    def total(key: str) -> int:
        return sum(run.get("counts", {}).get(key, 0) for run in runs)

    samples = {run["sample"] for run in runs}
    reads_in = sum(
        extraction.get(name, {}).get("n_reads", 0) for name in samples
    )

    if arm == "assembly":
        rows = [
            ("Reads in", reads_in,
             "spanning + context, everything RNA-Bloom2 is given"),
            ("Contigs assembled", total("assembled_transcripts"),
             "RNA-Bloom2 output"),
            ("Kept by the Nexus filter", total("query_sequences"),
             "min read support 3, min mapping quality 30"),
            ("Aligned", total("aligned"),
             "minimap2 -ax splice:hq against the masked reference"),
            ("Spliced", total("after_unspliced_filter"),
             "what remove-unspliced-rnas leaves for call-rna-vars"),
        ]
    else:
        rows = [
            ("Spanning reads available", total("spanning_reads_available"),
             "reads whose alignment covers a vaccine variant"),
            ("Handed to Exacto", total("query_sequences"),
             "capped per variant to fit call-rna-vars in memory"),
            ("Aligned", total("aligned"),
             "minimap2 against the masked reference"),
            ("Spliced", total("aligned"),
             "not applicable — this route runs no unspliced filter"),
        ]

    first = rows[0][1] or 1
    stages = []
    previous = None
    for label, n, note in rows:
        stages.append(
            {
                "label": label,
                "n": n,
                "note": note,
                "of_input": round(n / first, 4),
                # Retention against the previous stage is what shows *where* the
                # sequence goes, which a cumulative fraction alone hides.
                "of_previous": round(n / previous, 4) if previous else None,
            }
        )
        previous = n or None
    return stages


def _ladder(runs: list[dict]) -> dict[str, set[str]]:
    """Which variants this route reaches each rung for, in any sample."""
    reached: dict[str, set[str]] = {
        key: set() for key in ("covered", "rna_call", "proteoform", "residue", "epitope")
    }
    for run in runs:
        for variant_id, entry in run.get("variants", {}).items():
            rank = _RANK.get(entry.get("outcome"), 0)
            if entry.get("spanning_reads", 0) > 0:
                reached["covered"].add(variant_id)
            if rank >= _RANK["rna_only"]:
                reached["rna_call"].add(variant_id)
            if rank >= _RANK["proteoform"]:
                reached["proteoform"].add(variant_id)
            if entry.get("residue_confirmed"):
                reached["residue"].add(variant_id)
            if entry.get("n_matched_epitopes"):
                reached["epitope"].add(variant_id)
    return reached


RUNGS = [
    ("covered", "Reads cover the locus",
     "at least one read spans the mutation — nobody could call it otherwise"),
    ("rna_call", "RNA variant called",
     "Exacto called this exact locus and allele in the RNA"),
    ("proteoform", "Mutant protein translated",
     "a translated primary structure carries the mutation"),
    ("residue", "Correct residue",
     "the amino acid matches what the portal's annotation predicts"),
    ("epitope", "Vaccine epitope verbatim",
     "the manufactured pVACtools peptide appears inside the proteoform"),
]


def _quality(runs: list[dict]) -> dict[str, Any]:
    """How good the proteoforms are, not just how many.

    A frameshifted proteoform is the signature the reads route is expected to
    produce and the assembly route is expected to avoid: an indel miscalled in a
    homopolymer shifts the frame and everything downstream of the mutation is
    wrong. Exacto labels it, so it can be counted rather than assumed.
    """
    lengths: list[int] = []
    frameshifted = 0
    residue_ok = residue_checked = 0
    epitope_hit = epitope_possible = 0
    # How many candidate proteins a route hands you for one mutation. The reads
    # route gives one per read that carries the allele, so a downstream user has
    # to choose among them and Exacto does not rank them; the assembly route
    # gives about one. Recall and per-answer precision trade against each other
    # here, and the count is the visible form of that trade.
    per_variant: list[int] = []
    for run in runs:
        for entry in run.get("variants", {}).values():
            if entry.get("n_proteoforms"):
                per_variant.append(entry["n_proteoforms"])
            for form in entry.get("proteoforms", []):
                lengths.append(form.get("protein_length", 0))
                if form.get("frameshift"):
                    frameshifted += 1
            if entry.get("residue_confirmed") is not None:
                residue_checked += 1
                residue_ok += bool(entry["residue_confirmed"])
            if entry.get("n_vaccine_epitopes") and entry.get("n_proteoforms"):
                epitope_possible += 1
                epitope_hit += bool(entry.get("n_matched_epitopes"))
    lengths.sort()
    per_variant.sort()
    return {
        "n_proteoforms": len(lengths),
        "median_per_variant": per_variant[len(per_variant) // 2] if per_variant else None,
        "max_per_variant": per_variant[-1] if per_variant else None,
        "n_recovered_variant_runs": len(per_variant),
        "median_length": lengths[len(lengths) // 2] if lengths else None,
        "frameshift_fraction": round(frameshifted / len(lengths), 4) if lengths else None,
        "residue_ok": residue_ok,
        "residue_checked": residue_checked,
        "epitope_hit": epitope_hit,
        "epitope_possible": epitope_possible,
    }


def _vaf_profile(
    reached: dict[str, dict[str, set[str]]], variants: list[dict], n_variants: int
) -> list[dict] | None:
    """Where on the VAF scale each route stops working.

    If assembly is a consensus step, a subclonal mutation is a minority signal
    inside its own contig and should drop out first. That predicts a gradient,
    which is checkable against the portal's own ONT genotyping.
    """
    if len(reached) < 2:
        return None
    vaf_by_id = {}
    for variant in variants:
        values = [
            entry["vaf"]
            for entry in (variant.get("ont_expectation") or {}).values()
            if entry and entry.get("vaf") is not None
        ]
        if values:
            vaf_by_id[variant["variant_id"]] = max(values)

    reads = reached.get("reads", {}).get("proteoform", set())
    assembly = reached.get("assembly", {}).get("proteoform", set())
    groups = [
        ("Both routes", assembly & reads),
        ("Reads only", reads - assembly),
        ("Neither", set(vaf_by_id) - (reads | assembly)),
    ]
    rows = []
    for label, members in groups:
        values = sorted(vaf_by_id[v] for v in members if v in vaf_by_id)
        rows.append(
            {
                "label": label,
                "n": len(members),
                "median_vaf": round(values[len(values) // 2], 4) if values else None,
            }
        )
    return rows


def _supported(variant_id: str, entry: dict, sample_timepoint: str | None,
               platform: str | None, portal: dict) -> bool:
    """Was the allele in this sample's RNA at all?

    Either source suffices: Exacto called the variant de novo, or the portal's
    own genotyping of the same BAM counted alt reads. The portal only genotyped
    ONT, so for other platforms the first is the only source available — which
    can understate support but cannot invent it.
    """
    if entry.get("rna_variant_calls"):
        return True
    if platform == "ONT" and sample_timepoint:
        seen = (portal.get(variant_id) or {}).get(sample_timepoint)
        if seen and (seen.get("alt_reads") or 0) > 0:
            return True
    return False


def benchmark(payload: dict | None, variants: list[dict] | None) -> dict | None:
    """Sensitivity against specificity, per method and per sample.

    The two pull against each other and no single number captures the trade, so
    four are reported rather than one score:

    * **sensitivity** — recovered / mutations whose allele is actually in this
      sample's RNA. Scoring against all 37 would punish a method for mutations
      that are simply not expressed.
    * **residue precision** — of the mutations recovered, how many carry the
      amino acid the annotation predicts. A protein at the right codon in the
      wrong frame is not a recovered neoantigen.
    * **in-frame fraction** — of the proteoforms emitted, how many are not
      frameshifted. The direct measure of sequence integrity.
    * **candidates per mutation** — how many distinct proteins the method hands
      back for one mutation. Exacto does not rank them, so this is the work a
      downstream user inherits, and lower is better.
    """
    if not payload or not payload.get("runs"):
        return None
    portal = {
        v["variant_id"]: (v.get("ont_expectation") or {}) for v in (variants or [])
    }
    ok_runs = [r for r in payload["runs"] if r.get("status") == "ok"]
    if not ok_runs:
        return None

    rows = []
    for run in ok_runs:
        supported = recovered = residue_ok = residue_seen = 0
        consensus_ok = consensus_seen = 0
        candidates_total = candidates_correct = 0
        forms = frameshifted = 0
        epi_any = epi_consensus = epi_possible = 0
        per_variant: list[int] = []
        for variant_id, entry in run.get("variants", {}).items():
            if _supported(variant_id, entry, run.get("timepoint"),
                          run.get("platform"), portal):
                supported += 1
            rank = _RANK.get(entry.get("outcome"), 0)
            if rank >= _RANK["proteoform"]:
                recovered += 1
                per_variant.append(entry.get("n_proteoforms", 0))
                if entry.get("residue_confirmed") is not None:
                    residue_seen += 1
                    residue_ok += bool(entry["residue_confirmed"])
                if entry.get("consensus_residue_confirmed") is not None:
                    consensus_seen += 1
                    consensus_ok += bool(entry["consensus_residue_confirmed"])
                if entry.get("n_vaccine_epitopes"):
                    epi_possible += 1
                    epi_any += bool(entry.get("n_matched_epitopes"))
                    epi_consensus += bool(entry.get("n_consensus_matched_epitopes"))
                if entry.get("n_residue_correct") is not None:
                    candidates_total += entry.get("n_proteoforms", 0)
                    candidates_correct += entry["n_residue_correct"]
            for form in entry.get("proteoforms", []):
                forms += 1
                frameshifted += bool(form.get("frameshift"))
        per_variant.sort()
        method = run.get("method") or {}
        rows.append({
            "sample": run["sample"],
            "platform": run.get("platform"),
            "method": run.get("arm"),
            "method_label": method.get("label", run.get("arm")),
            "family": method.get("family"),
            "tool": method.get("tool"),
            "params": method.get("params") or {},
            "supported": supported,
            "recovered": recovered,
            "sensitivity": round(recovered / supported, 4) if supported else None,
            "residue_precision": round(residue_ok / residue_seen, 4)
            if residue_seen else None,
            # What a caller gets picking the modal translation, with no
            # knowledge of the answer. The honest single-answer number.
            "consensus_precision": round(consensus_ok / consensus_seen, 4)
            if consensus_seen else None,
            # And what one gets picking blindly: the fraction of all candidates
            # that are right.
            "candidate_precision": round(candidates_correct / candidates_total, 4)
            if candidates_total else None,
            # The whole manufactured peptide, not one residue. In any candidate,
            # and in the one a caller would pick.
            "epitope_any": round(epi_any / epi_possible, 4) if epi_possible else None,
            "epitope_consensus": round(epi_consensus / epi_possible, 4)
            if epi_possible else None,
            "inframe_fraction": round(1 - frameshifted / forms, 4) if forms else None,
            "candidates_per_variant": per_variant[len(per_variant) // 2]
            if per_variant else None,
            "seconds": run.get("seconds"),
        })

    # Rolled up per method across samples, which is what the headline compares.
    by_method: dict[str, dict] = {}
    for row in rows:
        agg = by_method.setdefault(row["method"], {
            "method": row["method"], "method_label": row["method_label"],
            "family": row["family"], "tool": row["tool"], "params": row["params"],
            "supported": 0, "recovered": 0, "samples": 0, "seconds": 0,
            "_inframe": [], "_precision": [], "_candidates": [],
            "_consensus": [], "_candprec": [], "_epiany": [], "_epicons": [],
        })
        agg["samples"] += 1
        agg["supported"] += row["supported"]
        agg["recovered"] += row["recovered"]
        agg["seconds"] += row["seconds"] or 0
        for key, source in (("_inframe", "inframe_fraction"),
                            ("_precision", "residue_precision"),
                            ("_consensus", "consensus_precision"),
                            ("_candprec", "candidate_precision"),
                            ("_epiany", "epitope_any"),
                            ("_epicons", "epitope_consensus"),
                            ("_candidates", "candidates_per_variant")):
            if row[source] is not None:
                agg[key].append(row[source])

    summary = []
    for agg in by_method.values():
        def mean(values):
            return round(sum(values) / len(values), 4) if values else None
        summary.append({
            k: v for k, v in agg.items() if not k.startswith("_")
        } | {
            "sensitivity": round(agg["recovered"] / agg["supported"], 4)
            if agg["supported"] else None,
            "residue_precision": mean(agg["_precision"]),
            "consensus_precision": mean(agg["_consensus"]),
            "candidate_precision": mean(agg["_candprec"]),
            "epitope_any": mean(agg["_epiany"]),
            "epitope_consensus": mean(agg["_epicons"]),
            "inframe_fraction": mean(agg["_inframe"]),
            "candidates_per_variant": mean(agg["_candidates"]),
        })
    summary.sort(key=lambda r: -(r["sensitivity"] or 0))
    return {"by_method": summary, "by_sample": rows}


def analyse(
    payload: dict | None, variants: list[dict] | None = None
) -> dict[str, Any] | None:
    """Per-route funnel, variant ladder, and what each route uniquely finds."""
    if not payload or not payload.get("runs"):
        return None

    n_variants = payload.get("n_variants") or len(payload.get("variants", []))
    extraction = payload.get("extraction") or {}
    ok_runs = [run for run in payload["runs"] if run.get("status") == "ok"]
    if not ok_runs:
        return None

    paths = []
    reached_by_arm: dict[str, dict[str, set[str]]] = {}
    for arm in PATH_ORDER:
        runs = [run for run in ok_runs if run.get("arm") == arm]
        if not runs:
            continue
        reached = _ladder(runs)
        reached_by_arm[arm] = reached
        meta = PATH_META.get(arm, {"label": arm, "subtitle": "", "assembler": None,
                                   "description": ""})
        paths.append(
            {
                "arm": arm,
                **meta,
                "quality": _quality(runs),
                "by_platform": {
                    platform: _quality(
                        [r for r in runs if r.get("platform") == platform]
                    )
                    for platform in sorted({r.get("platform") for r in runs if r.get("platform")})
                },
                "n_samples": len({run["sample"] for run in runs}),
                "seconds": sum(run.get("seconds") or 0 for run in runs),
                "stages": _stages(arm, runs, extraction),
                "ladder": [
                    {
                        "key": key,
                        "label": label,
                        "note": note,
                        "n": len(reached[key]),
                        "fraction": round(len(reached[key]) / n_variants, 4)
                        if n_variants
                        else 0,
                    }
                    for key, label, note in RUNGS
                ],
            }
        )

    comparison = None
    if len(reached_by_arm) == 2:
        first, second = PATH_ORDER
        a = reached_by_arm.get(first, {}).get("proteoform", set())
        b = reached_by_arm.get(second, {}).get("proteoform", set())
        comparison = {
            "a": first,
            "b": second,
            "a_label": PATH_META[first]["label"],
            "b_label": PATH_META[second]["label"],
            "both": sorted(_gene(v) for v in a & b),
            "a_only": sorted(_gene(v) for v in a - b),
            "b_only": sorted(_gene(v) for v in b - a),
            "neither": n_variants - len(a | b),
        }

    return {
        "n_variants": n_variants,
        "paths": paths,
        "comparison": comparison,
        "vaf_profile": _vaf_profile(reached_by_arm, variants or [], n_variants),
        "benchmark": benchmark(payload, variants),
        "n_assembly_methods": sum(
            1 for path in paths if path.get("assembler")
        ),
    }


def _gene(variant_id: str) -> str:
    """Variant ids are GENE-chrom-pos; the gene is what a reader recognises."""
    return variant_id.split("-")[0]
