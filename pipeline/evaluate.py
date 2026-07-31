"""Decide, per vaccine mutation and sample, whether Exacto recovered it.

The question this repo exists to answer is narrow: *does a translated mutant
protein sequence come out the other end?* So the verdict ladder is graded rather
than binary, and each rung is read straight out of an Exacto output file:

===============  =========================================================
``no_reads``     nothing in this sample's reads covers the locus — nobody
                 could have called it, so it is not counted against Exacto
``no_call``      reads cover it, Exacto called no RNA variant there
``rna_only``     Exacto called the variant at that exact locus and allele in
                 the RNA, but translated no mutant protein carrying it
``proteoform``   a translated primary structure carries the mutation
``peptide``      ...and ``call-peptide-vars`` emitted novel mutant peptides
===============  =========================================================

Every rung is keyed on the RNA variant call Exacto made at the mutation's exact
locus with its exact allele, never on ``integrate-vars`` output — that links a
DNA variant to anything within 10 kb of a transcript edge or 100 kb
intergenically, and scoring off it turns 5 recovered mutations into 28.

For missense mutations the residue Exacto produced is compared against the one
the portal's HGVS annotation predicts — a proteoform that carries *some* change
at the right codon but the wrong amino acid is reported as an unconfirmed
residue rather than quietly counted as a win. Where the portal publishes the
pVACtools epitope sequences the vaccine designs came from, the proteoform is
also searched for those peptides verbatim.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

from .build_reference import load_variants
from .config import ARMS, RESULTS_DIR, SAMPLES, SAMPLES_BY_NAME, samples_named
from .extract_reads import stats_path
from .run_exacto import EXACTO_DIR, as_graph_operation, collect_runs

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

AMINO_ACIDS = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

MISSENSE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

# How much of the translated protein to show either side of the mutation.
CONTEXT_RESIDUES = 12

OUTCOME_ORDER = ["no_reads", "no_call", "rna_only", "proteoform", "peptide"]


def expected_change(variant: dict) -> dict:
    """What the portal's protein annotation says this mutation should do."""
    protein_change = variant.get("protein_change") or ""
    match = MISSENSE.match(protein_change)
    if match:
        reference, position, alternate = match.groups()
        return {
            "kind": "missense",
            "protein_change": protein_change,
            "position": int(position),
            "ref_aa": AMINO_ACIDS.get(reference),
            "alt_aa": AMINO_ACIDS.get(alternate),
        }
    if "fs" in protein_change or variant.get("consequence") == "frameshift_variant":
        return {"kind": "frameshift", "protein_change": protein_change or None}
    if "del" in protein_change or variant.get("consequence") == "inframe_deletion":
        return {"kind": "inframe_deletion", "protein_change": protein_change or None}
    return {"kind": "other", "protein_change": protein_change or None}


def read_tsv(path: Path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def split_ids(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item for item in value.strip('"').split(",") if item}


def rna_calls_by_variant(path: Path, variants: list[dict]) -> dict[str, list[dict]]:
    """Match Exacto's de-novo RNA calls back to the vaccine mutations."""
    wanted: dict[tuple[str, int, int, str], dict] = {}
    for variant in variants:
        position_1, position_2, _, variant_type, sequence = as_graph_operation(variant)
        key = (variant["chrom"], position_1, position_2, variant_type)
        wanted[key] = {"variant": variant, "sequence": sequence.upper()}

    hits: dict[str, list[dict]] = {}
    if not path.exists():
        return hits

    for row in read_tsv(path):
        key = (
            row["chromosome_1"],
            int(row["position_1"]),
            int(row["position_2"]),
            row["variant_type"],
        )
        target = wanted.get(key)
        if target is None:
            continue
        called_sequence = (row.get("sequence") or "").strip('"').upper()
        if called_sequence != target["sequence"]:
            continue
        hits.setdefault(target["variant"]["variant_id"], []).append(
            {
                "rna_variant_call_id": row["variant_call_id"],
                "transcript_model_id": row["transcript_model_id"],
                "reference_transcript_ids": (
                    row.get("reference_transcript_ids") or ""
                ).strip('"'),
                # Exacto names the reads behind each call. Counting them gives
                # allele support derived from this run rather than read off the
                # portal — the only way to get a number for PacBio, which the
                # portal's variant table never genotyped.
                "n_supporting_reads": len(
                    split_ids(row.get("consensus_read_names"))
                ),
            }
        )
    return hits


def integrated_pairs(path: Path) -> dict[str, set[str]]:
    """dna_variant_call_id -> rna_variant_call_ids Exacto tied it to."""
    pairs: dict[str, set[str]] = {}
    if not path.exists():
        return pairs
    for row in read_tsv(path):
        pairs.setdefault(row["dna_variant_call_id"], set()).add(
            row["rna_variant_call_id"]
        )
    return pairs


def proteoforms_by_rna_call(
    path: Path, call_ids: set[str], references: dict[str, dict] | None = None
) -> dict[str, list[dict]]:
    """Pull the translated protein around each mutation.

    Keyed on *RNA* variant call ids, not DNA ones. Exacto's ``integrate-vars``
    links a DNA variant to any RNA variant within ``--max-exon-offset`` of an
    exon, ``--max-transcript-boundary-offset`` (10 kb by default) of a transcript
    edge, or ``--max-intergenic-distance`` (100 kb) otherwise — in one T1 arm
    only 19 of 3,359 integrations were exact. Scoring off those links would count
    a mutation as recovered whenever Exacto found *some* other variant in the
    neighbourhood. So the question asked here is the narrow one: did the
    RNA variant Exacto called at this exact locus, with this exact allele, change
    an amino acid in a protein it translated?

    The primary-structures TSV is one row per *nucleotide* per peptide, so it can
    run to tens of millions of rows. It is grouped by ``peptide_id``, so a single
    streaming pass that buffers one peptide at a time keeps memory flat regardless
    of how big the file gets.
    """
    if not path.exists():
        return {}

    by_call: dict[str, list[dict]] = {}

    def flush(peptide_id: int | None, residues: dict[int, str], hits: dict[str, dict]) -> None:
        if peptide_id is None or not hits:
            return
        ordered = [residues[index] for index in sorted(residues)]
        for call_id, entry in hits.items():
            indices = sorted(entry["amino_acid_indices"])
            start = max(0, indices[0] - CONTEXT_RESIDUES)
            end = min(len(ordered), indices[-1] + CONTEXT_RESIDUES + 1)
            by_call.setdefault(call_id, []).append(
                {
                    "peptide_id": peptide_id,
                    "transcript_model_id": entry["transcript_model_id"],
                    "reference_transcript_ids": entry["reference_transcript_ids"],
                    "protein_length": len(ordered),
                    "mutant_residue_indices": indices,
                    "mutant_residues": "".join(
                        ordered[index] for index in indices if index < len(ordered)
                    ),
                    "context": "".join(ordered[start:end]),
                    "context_start": start + 1,
                    "frameshift": sorted(entry["frameshift_states"]) != ["inframe"],
                    # Kept for epitope matching; stripped before serialising.
                    "protein": "".join(ordered),
                }
            )

    current: int | None = None
    residues: dict[int, str] = {}
    hits: dict[str, dict] = {}
    completed: set[int] = set()

    for row in read_tsv(path):
        if row["type"] != "base":
            continue
        peptide_id = int(row["peptide_id"])
        if peptide_id != current:
            flush(current, residues, hits)
            if current is not None:
                completed.add(current)
            if peptide_id in completed:
                raise SystemExit(
                    f"{path} is not grouped by peptide_id (saw {peptide_id} again) — "
                    "the streaming reader would drop rows"
                )
            current, residues, hits = peptide_id, {}, {}

        if row["codon_index"] == "0":
            residues[int(row["amino_acid_index"])] = row["amino_acid"]

        linked = split_ids(row["codon_rna_variant_call_ids"]) | split_ids(
            row["rna_variant_call_ids"]
        )
        if references is not None:
            for call_id in linked & call_ids:
                tally = references.setdefault(call_id, {})
                key = row["amino_acid_change"] or "?"
                tally[key] = tally.get(key, 0) + 1

        if row["amino_acid_change"] != "mutant":
            continue
        for call_id in linked & call_ids:
            entry = hits.setdefault(
                call_id,
                {
                    "amino_acid_indices": set(),
                    "reference_transcript_ids": (
                        row.get("reference_transcript_ids") or ""
                    ).strip('"'),
                    "transcript_model_id": row.get("transcript_model_id"),
                    "frameshift_states": set(),
                },
            )
            entry["amino_acid_indices"].add(int(row["amino_acid_index"]))
            entry["frameshift_states"].add(row["frameshift_state"])

    flush(current, residues, hits)
    return by_call


def peptides_by_rna_call(path: Path, call_ids: set[str]) -> dict[str, list[dict]]:
    """Novel mutant peptides, keyed on the RNA call that produced them."""
    if not path.exists():
        return {}
    found: dict[str, list[dict]] = {}
    for row in read_tsv(path):
        matched = split_ids(row.get("rna_variant_call_ids")) & call_ids
        for call_id in matched:
            found.setdefault(call_id, []).append(
                {
                    "sequence": row["mutant_peptide_sequence"],
                    "k": int(row["k"]),
                }
            )
    for call_id, peptides in found.items():
        unique = {peptide["sequence"]: peptide for peptide in peptides}
        found[call_id] = sorted(unique.values(), key=lambda item: item["sequence"])
    return found


def evaluate_arm(
    run: dict,
    variants: list[dict],
    spanning: dict[str, int],
    spanning_seen: dict[str, int] | None = None,
) -> dict[str, dict]:
    """Grade one sample/arm run against every vaccine mutation."""
    call_ids = {str(index): variant for index, variant in enumerate(variants, start=1)}
    outputs = run.get("outputs", {})

    rna_hits = rna_calls_by_variant(
        Path(outputs.get("rna_variant_calls", "/nonexistent")), variants
    )
    integrated = integrated_pairs(
        Path(outputs.get("integrated_variants", "/nonexistent"))
    )
    # The RNA call ids Exacto assigned to each vaccine variant's exact locus and
    # allele. Everything downstream is keyed off these, not off integrate-vars.
    rna_ids_by_variant = {
        variant_id: {hit["rna_variant_call_id"] for hit in hits}
        for variant_id, hits in rna_hits.items()
    }
    all_rna_ids = set().union(*rna_ids_by_variant.values()) if rna_ids_by_variant else set()

    # Why did a call that Exacto made produce no protein? Answering that needs
    # to know whether translate-structs referenced the call at all.
    references: dict[str, dict] = {}
    proteoforms_by_rna = proteoforms_by_rna_call(
        Path(outputs.get("primary_structures_tsv", "/nonexistent")),
        all_rna_ids,
        references,
    )
    peptides_by_rna = peptides_by_rna_call(
        Path(outputs.get("peptide_variants", "/nonexistent")), all_rna_ids
    )

    graded: dict[str, dict] = {}
    for call_id, variant in call_ids.items():
        variant_id = variant["variant_id"]
        expectation = expected_change(variant)
        spanning_reads = spanning.get(variant_id, 0)

        rna_ids = rna_ids_by_variant.get(variant_id, set())
        variant_proteoforms = _dedupe_proteoforms(
            [form for rna_id in rna_ids for form in proteoforms_by_rna.get(rna_id, [])]
        )
        variant_peptides = _dedupe_peptides(
            [item for rna_id in rna_ids for item in peptides_by_rna.get(rna_id, [])]
        )
        variant_rna = rna_hits.get(variant_id, [])
        # Alt-read support for this allele, as Exacto itself saw it. Deduped
        # across calls would be better, but Exacto does not tell us whether two
        # calls share reads, so this is a ceiling and is labelled as one.
        alt_reads = sum(hit.get("n_supporting_reads", 0) for hit in variant_rna)
        # Recorded for reference, not used to decide anything: integrate-vars is
        # far too permissive to mean "Exacto found this variant" (see
        # proteoforms_by_rna_call).
        variant_integrated = sorted(integrated.get(call_id, set()))

        if variant_peptides and variant_proteoforms:
            outcome = "peptide"
        elif variant_proteoforms:
            outcome = "proteoform"
        elif variant_rna:
            outcome = "rna_only"
        elif spanning_reads == 0:
            outcome = "no_reads"
        else:
            outcome = "no_call"

        # Three different questions, deliberately kept apart.
        #
        # residue_confirmed  -- is the right answer anywhere in the candidate
        #                       set? A ceiling. Generous, and it flatters a
        #                       method that emits many candidates: one correct
        #                       proteoform among a hundred wrong ones scores the
        #                       same as one correct proteoform on its own.
        # residue_precision  -- what fraction of the candidates are right? What
        #                       a user picking blindly would get.
        # consensus_*        -- what a user gets picking deterministically, with
        #                       no knowledge of the answer. See _consensus_form.
        residue_confirmed = None
        n_residue_correct = None
        consensus_residue_confirmed = None
        # Picked for every variant, not only missense ones: the epitope check
        # below applies whatever the consequence class.
        consensus = _consensus_form(variant_proteoforms)
        if expectation["kind"] == "missense" and variant_proteoforms:
            expected_aa = expectation["alt_aa"]
            correct = [
                form for form in variant_proteoforms
                if expected_aa and expected_aa in form["mutant_residues"]
            ]
            residue_confirmed = bool(correct)
            n_residue_correct = len(correct)
            if consensus is not None:
                consensus_residue_confirmed = bool(
                    expected_aa and expected_aa in consensus["mutant_residues"]
                )

        # Two strengths of evidence, both reported.
        #
        #   residue    one amino acid: does the codon under test change to what
        #              the annotation predicts.
        #   epitope    the entire manufactured peptide -- an 8-11mer for class I
        #              or up to a 17mer for class II -- found verbatim as a
        #              substring of the whole translated protein, not of the
        #              trimmed display context. Far stronger: it requires the
        #              sequence either side of the mutation to be right too,
        #              which is exactly what a frameshift destroys.
        epitope_hits = matched_epitopes(variant, variant_proteoforms)
        # And the same question asked of the single candidate a caller would
        # actually pick, rather than of the whole set.
        consensus_epitopes = (
            matched_epitopes(variant, [consensus]) if consensus else []
        )
        # One proteoform per variant, chosen without knowing the answer, with
        # the mutant interval and the vaccine peptide located inside it.
        consensus_peptide = vaccine_peptide(variant, consensus)

        graded[variant_id] = {
            "dna_variant_call_id": call_id,
            "outcome": outcome,
            "spanning_reads": spanning_reads,
            "spanning_reads_available": (spanning_seen or {}).get(
                variant_id, spanning_reads
            ),
            "expected": expectation,
            "residue_confirmed": residue_confirmed,
            "n_residue_correct": n_residue_correct,
            "consensus_residue_confirmed": consensus_residue_confirmed,
            "consensus_peptide_id": (consensus or {}).get("peptide_id"),
            "consensus_support": (consensus or {}).get("consensus_support"),
            "n_vaccine_epitopes": len(variant.get("vaccine_epitopes") or []),
            # Not truncated: the per-variant roll-up counts from this list, and
            # a variant can legitimately match dozens of overlapping epitopes.
            "matched_epitopes": epitope_hits,
            "n_matched_epitopes": len(epitope_hits),
            "n_consensus_matched_epitopes": len(consensus_epitopes),
            "consensus_vaccine_peptide": consensus_peptide,
            "consensus_proteoform": (
                {
                    key: value
                    for key, value in consensus.items()
                    if key != "protein"
                }
                if consensus
                else None
            ),
            "rna_variant_calls": variant_rna,
            "alt_reads_in_calls": alt_reads,
            # Rows in the primary-structures table that name this variant's RNA
            # calls, by what translate-structs said the amino acid did. Empty
            # means the call was never referenced; {"reference": n} means it was
            # seen and judged not to change the protein.
            "primary_structure_rows": {
                rna_id: references[rna_id] for rna_id in sorted(rna_ids)
                if rna_id in references
            },
            "integrated_rna_call_ids": variant_integrated,
            "proteoforms": [
                {key: value for key, value in form.items() if key != "protein"}
                for form in variant_proteoforms[:3]
            ],
            "n_proteoforms": len(variant_proteoforms),
            "mutant_peptides": variant_peptides[:10],
            "n_mutant_peptides": len(variant_peptides),
        }
    return graded


def _consensus_form(forms: list[dict]) -> dict | None:
    """The candidate a caller would pick with no knowledge of the answer.

    Take the modal protein sequence. Basecalling errors are independent between
    reads, so a spurious indel appears in the read that carried it and nowhere
    else, while the true sequence recurs across every read that covers the
    locus. Counting identical translations is therefore a consensus taken at the
    protein level rather than the nucleotide level — and unlike assembly, it
    happens *after* the allele has already been separated onto its own reads, so
    it cannot average a minority allele away.

    Ties break on the lowest peptide id, which is arbitrary but deterministic;
    a tie means the reads genuinely disagree and no rule saves you.

    This is a rule the harness applies, not something Exacto does. Exacto emits
    the candidates unranked.
    """
    if not forms:
        return None
    counts: dict[str, int] = {}
    for form in forms:
        protein = form.get("protein") or form.get("context") or ""
        counts[protein] = counts.get(protein, 0) + 1
    best = max(
        forms,
        key=lambda form: (
            counts.get(form.get("protein") or form.get("context") or "", 0),
            -form["peptide_id"],
        ),
    )
    return {
        **best,
        "consensus_support": counts.get(
            best.get("protein") or best.get("context") or "", 0
        ),
    }


def _dedupe_proteoforms(forms: list[dict]) -> list[dict]:
    """One entry per translated peptide, however many RNA calls point at it."""
    unique: dict[int, dict] = {}
    for form in forms:
        unique.setdefault(form["peptide_id"], form)
    return [unique[key] for key in sorted(unique)]


def _dedupe_peptides(peptides: list[dict]) -> list[dict]:
    unique = {peptide["sequence"]: peptide for peptide in peptides}
    return sorted(unique.values(), key=lambda item: item["sequence"])


def vaccine_peptide(variant: dict, form: dict | None) -> dict | None:
    """The whole vaccine peptide, as observed in this proteoform.

    The portal publishes pVACtools epitopes: 9- to 17-mers tiling the mutation,
    59 of them for some variants. Individually they are the predicted binders;
    together they span the long peptide a construct would carry. So rather than
    reconstructing that span from the epitope list, it is read straight out of
    the translated protein — the region from the first matched epitope's start
    to the last one's end. That is observed sequence, not a stitch, and it is
    the strongest statement this test can make: not "the right residue appeared"
    but "the entire manufactured peptide is here, verbatim, in a protein Exacto
    translated from the patient's RNA".

    The mutant interval is reported as an offset within that peptide, so a
    reader can see where the change sits inside what was injected.
    """
    epitopes = variant.get("vaccine_epitopes") or []
    protein = (form or {}).get("protein") or ""
    if not epitopes or not protein:
        return None

    start, end, found = len(protein), 0, 0
    for epitope in epitopes:
        sequence = epitope.get("sequence") or ""
        if not sequence:
            continue
        index = protein.find(sequence)
        if index < 0:
            continue
        found += 1
        start = min(start, index)
        end = max(end, index + len(sequence))
    if not found:
        return None

    indices = [i for i in form.get("mutant_residue_indices", []) if start <= i < end]
    return {
        "sequence": protein[start:end],
        "start": start + 1,
        "n_epitopes_matched": found,
        "n_epitopes_total": len(epitopes),
        "complete": found == len(epitopes),
        # Where the mutation sits inside the peptide, 0-based within it.
        "mutant_offsets": [i - start for i in indices],
        "peptide_id": form.get("peptide_id"),
    }


def matched_epitopes(variant: dict, proteoforms: list[dict]) -> list[dict]:
    """Vaccine epitopes that appear verbatim in a recovered proteoform.

    The strongest form of the question this repo asks: not "did Exacto produce
    a mutant protein at the right codon" but "is the peptide that was actually
    manufactured a substring of what Exacto translated". Only answerable for the
    loci where the portal published pVACtools epitope sequences.

    Matching is done against the whole translated protein, not the trimmed
    display context, since a 17-mer class II epitope easily runs past it.
    """
    epitopes = variant.get("vaccine_epitopes") or []
    if not epitopes or not proteoforms:
        return []

    hits: dict[str, dict] = {}
    for form in proteoforms:
        protein = form.get("protein") or form.get("context") or ""
        if not protein:
            continue
        for epitope in epitopes:
            sequence = epitope["sequence"]
            if sequence and sequence in protein:
                hits.setdefault(
                    sequence,
                    {
                        "sequence": sequence,
                        "mhc_class": epitope["mhc_class"],
                        "alleles": epitope["alleles"][:4],
                        "peptide_ids": [],
                    },
                )["peptide_ids"].append(form["peptide_id"])
    return sorted(hits.values(), key=lambda item: (len(item["sequence"]), item["sequence"]))


def best_outcome(outcomes: list[str]) -> str:
    return max(outcomes, key=OUTCOME_ORDER.index) if outcomes else "no_reads"


SCORED_DIR = RESULTS_DIR / "scored"


def score_samples(samples: list[str]) -> list[dict]:
    """Grade the runs on disk for these samples and cache them per sample.

    Scoring happens where the run happened, because the primary-structures TSVs
    it reads are far too big to hand between CI jobs.
    """
    variants = load_variants()
    runs = collect_runs(samples)
    if not runs:
        raise SystemExit(
            f"no run.json under {EXACTO_DIR} for {', '.join(samples)} — "
            "run pipeline.run_exacto first"
        )

    graded_runs = []
    for run in runs:
        stats_file = stats_path(SAMPLES_BY_NAME[run["sample"]])
        stats = json.loads(stats_file.read_text()) if stats_file.exists() else {}
        graded = (
            evaluate_arm(
                run,
                variants,
                stats.get("spanning_reads_by_variant", {}),
                stats.get("spanning_reads_seen_by_variant", {}),
            )
            if run["status"] == "ok"
            else {}
        )
        graded_runs.append(
            {
                "sample": run["sample"],
                "timepoint": run.get("timepoint"),
                "platform": run.get("platform"),
                "label": run.get("label", run["sample"]),
                "arm": run["arm"],
                "status": run["status"],
                "error": run.get("error"),
                "seconds": run.get("seconds"),
                "counts": run.get("counts", {}),
                "workarounds": run.get("workarounds", []),
                # Every step carries its full command so the site can offer it
                # verbatim; failures additionally carry the stderr tail, so a
                # bug report can be filed without digging through CI logs.
                "steps": [
                    {
                        key: step[key]
                        for key in (
                            "name",
                            "returncode",
                            "seconds",
                            "command",
                            *(("log_tail",) if step.get("returncode") else ()),
                        )
                        if key in step
                    }
                    for step in run.get("steps", [])
                ],
                "variants": graded,
            }
        )

    # One file per sample *and* arm. CI runs them as separate jobs, so two legs
    # of the same sample would otherwise write the same filename and the second
    # artifact downloaded would quietly replace the first.
    SCORED_DIR.mkdir(parents=True, exist_ok=True)
    for run in graded_runs:
        name = run["sample"]
        # Carry the read-extraction summary along: work/ does not survive between
        # CI jobs, and the site wants to show what went in.
        stats_file = stats_path(SAMPLES_BY_NAME[name])
        stats = json.loads(stats_file.read_text()) if stats_file.exists() else {}
        out_path = SCORED_DIR / f"{name}.{run['arm']}.json"
        out_path.write_text(
            json.dumps(
                {
                    "sample": name,
                    "arm": run["arm"],
                    "extraction": {
                        key: stats[key]
                        for key in (
                            "n_reads",
                            "n_spanning_reads",
                            "n_context_reads",
                            "mean_read_length",
                            "spanning_reads_per_variant_cap",
                            "context_reads_per_region_cap",
                        )
                        if key in stats
                    },
                    "runs": [run],
                },
                indent=2,
            )
            + "\n"
        )
        print(f"scored {name}/{run['arm']} -> {out_path}")
    return graded_runs


def merge() -> dict:
    """Combine every cached per-sample score into the final verdict."""
    variants = load_variants()
    graded_runs: list[dict] = []
    extraction: dict[str, dict] = {}
    for path in sorted(SCORED_DIR.glob("*.json")):
        scored = json.loads(path.read_text())
        graded_runs.extend(scored["runs"])
        if scored.get("extraction"):
            extraction[scored["sample"]] = scored["extraction"]
    if not graded_runs:
        raise SystemExit(f"no scored samples in {SCORED_DIR}")

    # Roll up to one verdict per variant across every sample and arm.
    summary = []
    for variant in variants:
        variant_id = variant["variant_id"]
        per_sample = {}
        for sample in SAMPLES:
            per_arm = {}
            for arm in ARMS:
                run = next(
                    (
                        item
                        for item in graded_runs
                        if item["sample"] == sample.name and item["arm"] == arm
                    ),
                    None,
                )
                if run and variant_id in run["variants"]:
                    per_arm[arm] = run["variants"][variant_id]
            per_sample[sample.name] = {
                "arms": per_arm,
                "outcome": best_outcome(
                    [entry["outcome"] for entry in per_arm.values()]
                ),
            }
        overall = best_outcome(
            [entry["outcome"] for entry in per_sample.values() if entry["arms"]]
        )
        residue_confirmed = None
        checks = [
            arm_entry.get("residue_confirmed")
            for entry in per_sample.values()
            for arm_entry in entry["arms"].values()
            if arm_entry.get("residue_confirmed") is not None
        ]
        if checks:
            residue_confirmed = any(checks)

        epitopes = sorted(
            {
                hit["sequence"]
                for entry in per_sample.values()
                for arm_entry in entry["arms"].values()
                for hit in arm_entry.get("matched_epitopes", [])
            },
            key=lambda sequence: (len(sequence), sequence),
        )
        summary.append(
            {
                "variant_id": variant_id,
                "gene": variant["gene"],
                "protein_change": variant["protein_change"],
                "expected": expected_change(variant),
                "outcome": overall,
                "residue_confirmed": residue_confirmed,
                "n_vaccine_epitopes": len(variant.get("vaccine_epitopes") or []),
                "matched_epitopes": epitopes,
                "samples": per_sample,
            }
        )

    recovered = [item for item in summary if item["outcome"] in ("proteoform", "peptide")]
    with_epitopes = [item for item in summary if item["n_vaccine_epitopes"]]
    epitope_confirmed = [item for item in with_epitopes if item["matched_epitopes"]]
    # Of the proteoforms that came back, how many carry the residue the vaccine
    # was designed around? A translated protein at the right codon in the wrong
    # frame is not a recovered neoantigen.
    residue_checkable = [
        item for item in recovered if item["residue_confirmed"] is not None
    ]
    residue_confirmed = [item for item in residue_checkable if item["residue_confirmed"]]
    testable = [
        item
        for item in summary
        if any(entry["outcome"] != "no_reads" for entry in item["samples"].values())
    ]

    payload = {
        "n_variants": len(summary),
        "n_testable": len(testable),
        "n_recovered": len(recovered),
        "recovered_genes": sorted(item["gene"] for item in recovered),
        "n_residue_checkable": len(residue_checkable),
        "n_residue_confirmed": len(residue_confirmed),
        "n_with_vaccine_epitopes": len(with_epitopes),
        "n_epitope_confirmed": len(epitope_confirmed),
        "epitope_confirmed_genes": sorted(item["gene"] for item in epitope_confirmed),
        "outcome_counts": {
            outcome: sum(1 for item in summary if item["outcome"] == outcome)
            for outcome in OUTCOME_ORDER
        },
        "extraction": extraction,
        "runs": graded_runs,
        "variants": summary,
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / "exacto_results.json"
    out_path.write_text(json.dumps(payload, indent=2) + "\n")

    print(f"{len(recovered)}/{len(testable)} testable vaccine mutations recovered "
          f"as mutant protein sequences ({len(summary)} total)")
    print(f"{len(residue_confirmed)}/{len(residue_checkable)} of the recovered "
          "proteoforms carry the amino acid the annotation predicts")
    print(f"{len(epitope_confirmed)}/{len(with_epitopes)} of the mutations with a "
          "published vaccine epitope had that exact peptide inside the proteoform")
    for outcome in OUTCOME_ORDER:
        print(f"  {outcome:12} {payload['outcome_counts'][outcome]}")
    print(f"-> {out_path}")
    return payload


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        nargs="*",
        help=(
            "Score these sequencing samples from the Exacto outputs on disk, "
            "e.g. T1-ONT T1-PacBio (default: all)."
        ),
    )
    parser.add_argument(
        "--merge-only",
        action="store_true",
        help="Skip scoring and just combine the cached per-sample results.",
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        help="Score only; leave the final verdict to a later --merge-only pass.",
    )
    args = parser.parse_args()

    if not args.merge_only:
        score_samples([s.name for s in samples_named(args.samples)])
    if not args.no_merge:
        merge()


if __name__ == "__main__":
    main()
