"""Checks on the streaming reader for Exacto's primary-structures TSV.

That file is one row per nucleotide per translated peptide and can run to tens
of millions of rows, so it is read in a single grouped pass that buffers one
peptide at a time. These pin down that the pass finds the right residues, keeps
the surrounding context, and refuses to guess if the file ever stops being
grouped by peptide.

It keys on RNA variant call ids: the DNA ids come from ``integrate-vars``, which
is far too permissive to base a verdict on.
"""

from __future__ import annotations

import pytest

from pipeline.evaluate import CONTEXT_RESIDUES, proteoforms_by_rna_call

COLUMNS = [
    "peptide_id", "primary_structure_index", "type", "amino_acid",
    "amino_acid_index", "codon_index", "nucleotide", "transcript_model_id",
    "reference_transcript_ids", "transcript_structure_index", "read_start",
    "read_end", "net_variant_nucleotides_count", "frameshift_state",
    "rna_variant_call_ids", "dna_variant_call_ids", "codon_rna_variant_call_ids",
    "codon_dna_variant_call_ids", "frameshift_rna_variant_call_ids",
    "frameshift_dna_variant_call_ids", "amino_acid_change",
]


def write_peptides(path, peptides):
    """Render peptides as Exacto does: three nucleotide rows per residue.

    ``peptides`` maps peptide_id -> (protein, {residue_index: rna_call_id}).
    """
    lines = ["\t".join(COLUMNS)]
    for peptide_id, (protein, mutations) in peptides.items():
        index = 0
        for residue_index, residue in enumerate(protein):
            call_id = mutations.get(residue_index)
            for codon_index in range(3):
                row = {name: "" for name in COLUMNS}
                row.update(
                    peptide_id=str(peptide_id),
                    primary_structure_index=str(index),
                    type="base",
                    amino_acid=residue,
                    amino_acid_index=str(residue_index),
                    codon_index=str(codon_index),
                    nucleotide="a",
                    transcript_model_id="7",
                    reference_transcript_ids="ENST00000000001.1",
                    frameshift_state="inframe",
                    codon_rna_variant_call_ids=call_id or "",
                    amino_acid_change="mutant" if call_id else "reference",
                )
                lines.append("\t".join(row[name] for name in COLUMNS))
                index += 1
    path.write_text("\n".join(lines) + "\n")
    return path


def test_finds_the_mutant_residue(tmp_path):
    protein = "MAAKSDGRLKMKKSSDVAFTPLQNSDNSGSVQG"
    path = write_peptides(tmp_path / "ps.tsv", {1: (protein, {17: "3"})})

    found = proteoforms_by_rna_call(path, {"3"})

    assert list(found) == ["3"]
    (form,) = found["3"]
    assert form["mutant_residue_indices"] == [17]
    assert form["mutant_residues"] == protein[17]
    assert form["protein_length"] == len(protein)
    assert not form["frameshift"]


def test_context_window_is_centred_and_bounded(tmp_path):
    protein = "M" + "ACDEFGHIKL" * 6
    path = write_peptides(tmp_path / "ps.tsv", {1: (protein, {40: "1"})})

    (form,) = proteoforms_by_rna_call(path, {"1"})["1"]

    assert form["context"] == protein[40 - CONTEXT_RESIDUES : 40 + CONTEXT_RESIDUES + 1]
    assert form["context_start"] == 40 - CONTEXT_RESIDUES + 1


def test_context_clamps_at_the_protein_start(tmp_path):
    protein = "MACDEFGHIKL"
    path = write_peptides(tmp_path / "ps.tsv", {1: (protein, {1: "1"})})

    (form,) = proteoforms_by_rna_call(path, {"1"})["1"]

    assert form["context"] == protein
    assert form["context_start"] == 1


def test_ignores_peptides_without_a_matching_call(tmp_path):
    path = write_peptides(
        tmp_path / "ps.tsv",
        {1: ("MACDEFG", {2: "9"}), 2: ("MKKLLL", {})},
    )

    assert proteoforms_by_rna_call(path, {"1", "2"}) == {}


def test_collects_one_entry_per_peptide_carrying_the_call(tmp_path):
    path = write_peptides(
        tmp_path / "ps.tsv",
        {1: ("MACDEFG", {2: "4"}), 2: ("MACDEFGH", {3: "4"})},
    )

    found = proteoforms_by_rna_call(path, {"4"})

    assert [form["peptide_id"] for form in found["4"]] == [1, 2]


def test_missing_file_is_not_an_error(tmp_path):
    assert proteoforms_by_rna_call(tmp_path / "absent.tsv", {"1"}) == {}


def test_ungrouped_file_is_rejected_rather_than_silently_truncated(tmp_path):
    path = write_peptides(tmp_path / "ps.tsv", {1: ("MACD", {1: "1"})})
    # Splice a second block for peptide 1 after peptide 2 has been flushed.
    lines = path.read_text().splitlines()
    header, body = lines[0], lines[1:]
    fields = body[0].split("\t")
    fields[0] = "2"
    path.write_text("\n".join([header, *body, "\t".join(fields), *body]) + "\n")

    with pytest.raises(SystemExit, match="not grouped by peptide_id"):
        proteoforms_by_rna_call(path, {"1"})


# --------------------------------------------------------------------------
# Resolving many candidates to one, without knowing the answer
# --------------------------------------------------------------------------

from pipeline.evaluate import _consensus_form


def _form(peptide_id, protein, residues):
    return {"peptide_id": peptide_id, "protein": protein, "mutant_residues": residues}


def test_consensus_picks_the_recurring_sequence_not_the_erroneous_one():
    """Errors are independent between reads; the true sequence recurs.

    Three reads agree on the same translation and two disagree in different
    ways, which is what a basecalling indel looks like: present in the read that
    carried it and nowhere else.
    """
    forms = [
        _form(1, "MKTAYIAKQR", "T"),
        _form(2, "MKTAYIAKQR", "T"),
        _form(3, "MKTGYIAKQR", "G"),   # one read's error
        _form(4, "MKTAYIAKQR", "T"),
        _form(5, "MKTWYIAKQR", "W"),   # a different read's error
    ]

    consensus = _consensus_form(forms)

    assert consensus["protein"] == "MKTAYIAKQR"
    assert consensus["consensus_support"] == 3


def test_consensus_is_deterministic_under_a_tie():
    forms = [_form(7, "AAA", "A"), _form(3, "BBB", "B")]

    assert _consensus_form(forms)["peptide_id"] == _consensus_form(forms)["peptide_id"]


def test_consensus_of_nothing_is_none():
    assert _consensus_form([]) is None
