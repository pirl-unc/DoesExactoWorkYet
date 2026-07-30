"""Checks on the two translations this repo does by hand.

Exacto brackets a variant with the untouched bases either side rather than using
VCF's anchor-base convention, and getting that off by one would quietly produce
a callset that never integrates with anything. These pin it down.
"""

from __future__ import annotations

import pytest

from pipeline.evaluate import expected_change
from pipeline.run_exacto import as_graph_operation


def variant(**overrides):
    base = {
        "gene": "TEST",
        "chrom": "chr1",
        "pos": 1000,
        "ref": "C",
        "alt": "T",
        "protein_change": None,
        "consequence": None,
    }
    return {**base, **overrides}


class TestGraphOperation:
    def test_snv_brackets_the_mutated_base(self):
        # A C>T at 1000 keeps 999 and 1001 and puts a T between them.
        assert as_graph_operation(variant()) == (999, 1001, 1, "SNV", "T")

    def test_deletion_of_one_base(self):
        # VCF AG>A at 1000 deletes the G at 1001.
        result = as_graph_operation(variant(ref="AG", alt="A"))
        assert result == (1000, 1002, 1, "DEL", "")

    def test_deletion_spans_every_removed_base(self):
        # CCTGGGCTACTGTGTGTTCAATA>C removes 22 bases starting at 1001.
        ref = "CCTGGGCTACTGTGTGTTCAATA"
        position_1, position_2, size, variant_type, sequence = as_graph_operation(
            variant(ref=ref, alt="C")
        )
        assert (position_1, position_2) == (1000, 1023)
        assert size == len(ref) - 1 == 22
        assert (variant_type, sequence) == ("DEL", "")
        # Exacto derives the size from the gap between the bracketing bases.
        assert position_2 - position_1 - 1 == size

    def test_insertion(self):
        result = as_graph_operation(variant(ref="A", alt="AGGT"))
        assert result == (1000, 1001, 3, "INS", "GGT")

    def test_mnv(self):
        position_1, position_2, size, variant_type, sequence = as_graph_operation(
            variant(ref="CT", alt="GA")
        )
        assert (position_1, position_2, size) == (999, 1002, 2)
        assert (variant_type, sequence) == ("MNV", "GA")

    def test_unencodable_variant_is_loud(self):
        with pytest.raises(SystemExit):
            as_graph_operation(variant(ref="CT", alt="GAC"))


class TestExpectedChange:
    def test_missense_gives_the_predicted_residue(self):
        result = expected_change(variant(protein_change="p.Arg371Lys"))
        assert result["kind"] == "missense"
        assert result["position"] == 371
        assert (result["ref_aa"], result["alt_aa"]) == ("R", "K")

    def test_stop_gained(self):
        assert expected_change(variant(protein_change="p.Gln100Ter"))["alt_aa"] == "*"

    def test_frameshift(self):
        assert expected_change(variant(protein_change="p.Ser775fs"))["kind"] == "frameshift"

    def test_frameshift_without_a_position(self):
        result = expected_change(
            variant(protein_change="fs", consequence="frameshift_variant")
        )
        assert result["kind"] == "frameshift"

    def test_inframe_deletion(self):
        result = expected_change(
            variant(protein_change="p.Ala197_Lys201del", consequence="inframe_deletion")
        )
        assert result["kind"] == "inframe_deletion"

    def test_unannotated(self):
        assert expected_change(variant())["kind"] == "other"


def test_every_vaccine_variant_encodes():
    """The real callset must round-trip — a new variant type would break the run."""
    from pipeline.build_reference import load_variants

    try:
        variants = load_variants()
    except FileNotFoundError:
        pytest.skip("run `python -m pipeline.fetch_osteosarc` first")

    for record in variants:
        position_1, position_2, size, variant_type, sequence = as_graph_operation(record)
        assert position_1 < position_2, record["variant_id"]
        assert size >= 1, record["variant_id"]
        if variant_type == "SNV":
            assert len(sequence) == 1
        if variant_type == "DEL":
            assert sequence == ""
            assert position_2 - position_1 - 1 == size
