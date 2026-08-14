# Copyright 2026 Nathan Fargo
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

import pytest

from geneops import (
    Interval,
    PAM,
    PAMSite,
    Protospacer,
    SequenceContext,
    Spacer,
    TargetSite,
    find_guides,
    get_nuclease,
    reverse_complement,
)

SPACER = "ATGCACGTTAGCTACGATCA"


class TestSpacer:
    def test_is_a_location_free_normalized_guide_sequence(self):
        spacer = Spacer("augc")

        assert spacer.sequence == "ATGC"
        assert spacer.dna_sequence == "ATGC"
        assert spacer.rna_sequence == "AUGC"
        assert spacer.gc_content == pytest.approx(0.5)
        assert len(spacer) == 4

    @pytest.mark.parametrize("sequence", ["", "ACGN", "ACG-"])
    def test_rejects_non_constructible_sequences(self, sequence):
        with pytest.raises(ValueError):
            Spacer(sequence)

    def test_rejects_non_string_sequence(self):
        with pytest.raises(TypeError, match="Spacer sequence"):
            Spacer(123)  # type: ignore[arg-type]


class TestSequenceContext:
    def test_normalizes_and_extracts_global_interval(self):
        context = SequenceContext(
            "acnug",
            Interval(100, 105),
            reference="GRCh38",
            contig="chr7",
        )

        assert context.sequence == "ACNTG"
        assert context.extract(Interval(101, 104)) == "CNT"
        assert context.extract(Interval(101, 104), "-") == "ANG"
        assert context.reference == "GRCh38"
        assert context.contig == "chr7"

    def test_interval_length_must_match_sequence(self):
        with pytest.raises(ValueError, match="interval length"):
            SequenceContext("ACGT", Interval(10, 15))

    def test_interval_must_be_interval_object(self):
        with pytest.raises(TypeError, match="SequenceContext.interval"):
            SequenceContext("AC", (0, 2))  # type: ignore[arg-type]

    def test_extract_rejects_interval_outside_context(self):
        context = SequenceContext("ACGT", Interval(10, 14))
        with pytest.raises(ValueError, match="outside context"):
            context.extract(Interval(9, 11))

    def test_flanks_are_relative_to_requested_strand(self):
        context = SequenceContext("AACCGGTT", Interval(10, 18))

        assert context.extract_flanked(
            Interval(12, 16),
            upstream=2,
            downstream=1,
        ) == "AACCGGT"
        assert context.extract_flanked(
            Interval(12, 16),
            upstream=2,
            downstream=1,
            strand="-",
        ) == "AACCGGT"


class TestProtospacer:
    def test_extracts_sequence_in_pam_strand_orientation(self):
        context = SequenceContext(
            "AAGCTTCA",
            Interval(100, 108),
            reference="GRCh38",
            contig="chr7",
        )

        plus = Protospacer(context, Interval(101, 106), "+")
        minus = Protospacer(context, Interval(101, 106), "-")

        assert plus.sequence == "AGCTT"
        assert minus.sequence == "AAGCT"
        assert plus.target_strand == "-"
        assert minus.target_strand == "+"
        assert plus.context.reference == "GRCh38"
        assert len(plus) == 5

    def test_rejects_unknown_bases(self):
        context = SequenceContext("ACNT", Interval(0, 4))
        with pytest.raises(ValueError, match="unknown bases"):
            Protospacer(context, Interval(0, 4), "+")

    def test_rejects_empty_interval(self):
        context = SequenceContext("ACGT", Interval(0, 4))
        with pytest.raises(ValueError, match="must not be empty"):
            Protospacer(context, Interval(2, 2), "+")


class TestPAMSite:
    def test_validates_observed_sequence_against_pam(self):
        context = SequenceContext("CCAAGGAA", Interval(100, 108))
        site = PAMSite(
            context,
            Interval(103, 106),
            "+",
            PAM("NGG", "3'"),
        )

        assert site.sequence == "AGG"
        assert site.pam.pattern == "NGG"

    def test_reads_reverse_strand_pam_in_its_own_orientation(self):
        context = SequenceContext("AACCAAAA", Interval(100, 108))
        site = PAMSite(
            context,
            Interval(102, 105),
            "-",
            PAM("NGG", "3'"),
        )

        assert site.sequence == "TGG"

    def test_rejects_non_matching_observed_sequence(self):
        context = SequenceContext("CCAAAACC", Interval(100, 108))
        with pytest.raises(ValueError, match="does not match PAM"):
            PAMSite(
                context,
                Interval(102, 105),
                "+",
                PAM("NGG", "3'"),
            )


class TestTargetSite:
    @pytest.fixture
    def valid_arguments(self):
        return {
            "context": SequenceContext(
                "A" * 20 + "AGG",
                Interval(0, 23),
            ),
            "protospacer_interval": Interval(0, 20),
            "pam_interval": Interval(20, 23),
            "pam_strand": "+",
            "nuclease": get_nuclease("SpCas9"),
        }

    @pytest.mark.parametrize(
        ("field", "value", "message"),
        [
            ("context", "not-context", "TargetSite.context"),
            (
                "protospacer_interval",
                (0, 20),
                "TargetSite.protospacer_interval",
            ),
            ("pam_interval", (20, 23), "TargetSite.pam_interval"),
            ("nuclease", "SpCas9", "TargetSite.nuclease"),
        ],
    )
    def test_rejects_wrong_object_types(
        self,
        valid_arguments,
        field,
        value,
        message,
    ):
        arguments = {**valid_arguments, field: value}
        with pytest.raises(TypeError, match=message):
            TargetSite(**arguments)

    def test_find_guides_preserves_reference_coordinates_and_context(self):
        target = "A" * 10 + SPACER + "TGG" + "A" * 10
        context = SequenceContext(
            target,
            Interval(100, 100 + len(target)),
            reference="GRCh38",
            contig="chr7",
        )

        guide = next(
            guide
            for guide in find_guides(
                context,
                get_nuclease("SpCas9"),
                pam_strand="+",
            )
            if guide.pam_interval == Interval(130, 133)
        )

        assert guide.context is context
        assert guide.protospacer_interval == Interval(110, 130)
        assert isinstance(guide.target_site, TargetSite)
        assert isinstance(guide.target_site.spacer, Spacer)
        assert isinstance(guide.target_site.protospacer, Protospacer)
        assert isinstance(guide.target_site.pam_site, PAMSite)
        assert guide.target_site.spacer.sequence == SPACER
        assert guide.target_site.protospacer_sequence == SPACER
        assert guide.target_site.pam_sequence == "TGG"
        assert guide.target_site.sequence_context(
            upstream=4,
            downstream=3,
        ) == "A" * 4 + SPACER + "TGG" + "A" * 3

    def test_reverse_target_context_is_reported_in_pam_orientation(self):
        oriented_site = SPACER + "TGG"
        target = "A" * 10 + reverse_complement(oriented_site) + "A" * 10
        context = SequenceContext(target, Interval(500, 500 + len(target)))

        guide = next(
            guide
            for guide in find_guides(
                context,
                get_nuclease("SpCas9"),
                pam_strand="-",
            )
            if guide.pam_interval == Interval(510, 513)
        )

        assert guide.target_site.protospacer_sequence == SPACER
        assert guide.target_site.pam_sequence == "TGG"
        assert guide.target_site.sequence_context(
            upstream=4,
            downstream=3,
        ) == "T" * 4 + oriented_site + "T" * 3

    def test_context_must_match_guide_sequence(self):
        context = SequenceContext("A" * 20 + "AGG", Interval(0, 23))

        from geneops import Guide

        with pytest.raises(ValueError, match="Guide sequence"):
            Guide(
                sequence="C" * 20,
                nuclease=get_nuclease("SpCas9"),
                pam="AGG",
                pam_interval=Interval(20, 23),
                pam_strand="+",
                context=context,
            )

    def test_direct_guide_gets_minimal_anonymous_target_site(self):
        from geneops import Guide

        guide = Guide(
            sequence=SPACER,
            nuclease=get_nuclease("SpCas9"),
            pam="TGG",
            pam_interval=Interval(20, 23),
            pam_strand="+",
        )

        assert guide.context is None
        assert guide.target_site.context.interval == Interval(0, 23)
        assert guide.target_site.protospacer_sequence == SPACER
