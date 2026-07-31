"""Checks on the strictest rung of the test: does the proteoform contain the
peptide that was actually put in the vaccine?

The portal publishes the curated pVACtools run the vaccine designs were picked
from, so for the loci it covers this is an exact substring question rather than
an argument about annotation.
"""

from __future__ import annotations

from pipeline.evaluate import matched_epitopes

PROTEIN = "MSTAENECGKSFGRSCHLIQHQTIHTGEKPYKCNE"


def variant(*sequences, **overrides):
    record = {
        "gene": "ZNF436",
        "vaccine_epitopes": [
            {
                "sequence": sequence,
                "mhc_class": "I",
                "alleles": ["HLA-A*01:01", "HLA-B*07:02"],
                "wild_type": None,
            }
            for sequence in sequences
        ],
    }
    return {**record, **overrides}


def proteoform(protein=PROTEIN, peptide_id=1):
    return {"peptide_id": peptide_id, "protein": protein, "context": protein[:20]}


def test_finds_an_epitope_inside_the_proteoform():
    hits = matched_epitopes(variant("CGKSFGRSC"), [proteoform()])

    assert [hit["sequence"] for hit in hits] == ["CGKSFGRSC"]
    assert hits[0]["mhc_class"] == "I"
    assert hits[0]["peptide_ids"] == [1]


def test_reports_nothing_when_the_epitope_is_absent():
    # One residue off — the wild-type version of the same epitope.
    assert matched_epitopes(variant("CGKSFGRSS"), [proteoform()]) == []


def test_matches_against_the_whole_protein_not_the_display_context():
    # A class II 17-mer runs past the trimmed context window.
    long_epitope = "NECGKSFGRSCHLIQHQ"
    assert long_epitope in PROTEIN
    hits = matched_epitopes(variant(long_epitope), [proteoform()])
    assert [hit["sequence"] for hit in hits] == [long_epitope]


def test_orders_hits_shortest_first_and_dedups_across_proteoforms():
    hits = matched_epitopes(
        variant("CHLIQHQTI", "CGKSFGRSCHLIQHQTI", "CGKSFGRSC"),
        [proteoform(peptide_id=1), proteoform(peptide_id=2)],
    )

    assert [len(hit["sequence"]) for hit in hits] == [9, 9, 17]
    assert all(hit["peptide_ids"] == [1, 2] for hit in hits)


def test_no_epitopes_published_means_no_claim():
    assert matched_epitopes(variant(), [proteoform()]) == []


def test_no_proteoform_means_no_claim():
    assert matched_epitopes(variant("CGKSFGRSC"), []) == []


def test_falls_back_to_context_when_the_protein_is_absent():
    form = {"peptide_id": 3, "context": "NECGKSFGRSCHLIQ"}
    hits = matched_epitopes(variant("CGKSFGRSC"), [form])
    assert [hit["sequence"] for hit in hits] == ["CGKSFGRSC"]


# --------------------------------------------------------------------------
# The whole vaccine peptide, located inside the chosen proteoform
# --------------------------------------------------------------------------

from pipeline.evaluate import vaccine_peptide


def test_vaccine_peptide_spans_from_first_epitope_to_last():
    """The span is read out of observed sequence, not reconstructed."""
    protein = "MAAAKQWERTYCHLIQHQTIPPPGGG"
    variant = {"vaccine_epitopes": [
        {"sequence": "QWERTYC"},        # protein[5:12]
        {"sequence": "CHLIQHQTI"},      # protein[11:20]
    ]}
    form = {"protein": protein, "mutant_residue_indices": [12], "peptide_id": 1}

    peptide = vaccine_peptide(variant, form)

    assert peptide["sequence"] == protein[5:20] == "QWERTYCHLIQHQTI"
    assert peptide["n_epitopes_matched"] == 2
    assert peptide["complete"] is True
    # The mutation is at protein index 12, so 7 residues into the peptide.
    assert peptide["mutant_offsets"] == [7]
    assert peptide["start"] == 6  # 1-based position in the protein


def test_vaccine_peptide_reports_partial_matches_as_incomplete():
    variant = {"vaccine_epitopes": [{"sequence": "AAAA"}, {"sequence": "ZZZZ"}]}
    form = {"protein": "GGAAAAGG", "mutant_residue_indices": [3], "peptide_id": 2}

    peptide = vaccine_peptide(variant, form)

    assert peptide["n_epitopes_matched"] == 1
    assert peptide["complete"] is False


def test_no_peptide_when_nothing_matches():
    variant = {"vaccine_epitopes": [{"sequence": "ZZZZ"}]}
    assert vaccine_peptide(variant, {"protein": "AAAA"}) is None
    assert vaccine_peptide({"vaccine_epitopes": []}, {"protein": "AAAA"}) is None


# --------------------------------------------------------------------------
# The expectation check has to be defined for every consequence class
# --------------------------------------------------------------------------

def test_frameshift_expectation_is_checkable():
    """A frameshift has no single expected residue; being frameshifted is the check.

    Before this, residue_confirmed was computed only for missense, so a
    frameshift that translated perfectly scored as unverifiable and dropped out
    of the funnel's residue rung — indistinguishable from a failure.
    """
    from pipeline.evaluate import expected_change

    assert expected_change({"protein_change": "p.Ser775fs",
                            "consequence": "frameshift_variant"})["kind"] == "frameshift"
    assert expected_change({"protein_change": "p.Ala197_Lys201del",
                            "consequence": "inframe_deletion"})["kind"] == "inframe_deletion"
    assert expected_change({"protein_change": "p.Ser207Cys",
                            "consequence": "missense_variant"})["kind"] == "missense"
