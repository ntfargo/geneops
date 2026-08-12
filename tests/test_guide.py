# Copyright 2026 Nathan Fargo
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest

from geneops.coordinates import Interval
from geneops.guide import (
    Guide,
    find_guides,
    guide_gc_content,
    normalize,
    normalize_target,
    reverse_complement,
    validate_guide_sequence,
)
from geneops.nuclease import Nuclease, get_nuclease


@pytest.fixture
def spcas9() -> Nuclease:
    return get_nuclease("SpCas9")


@pytest.fixture
def ascas12a() -> Nuclease:
    return get_nuclease("AsCas12a")


def make_guide(
    nuclease: Nuclease,
    *,
    sequence: str | None = None,
    pam: str | None = None,
    pam_interval: Interval | None = None,
    pam_strand: str = "+",
) -> Guide:
    if sequence is None:
        sequence = "A" * nuclease.spacer_length
    if pam is None:
        pam = "AGG" if nuclease.name == "SpCas9" else "TTTG"
    if pam_interval is None:
        pam_interval = Interval(30, 30 + len(pam))
    return Guide(
        sequence=sequence,
        nuclease=nuclease,
        pam=pam,
        pam_interval=pam_interval,
        pam_strand=pam_strand,  # type: ignore[arg-type]
    )


class TestNormalize:
    def test_uppercase(self):
        assert normalize("acgt") == "ACGT"

    def test_u_to_t(self):
        assert normalize("AUGC") == "ATGC"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            normalize("")

    @pytest.mark.parametrize("sequence", ["ACGX", "AC GT", "ACG1T"])
    def test_invalid_guide_bases_raise(self, sequence):
        with pytest.raises(ValueError, match="Invalid bases"):
            normalize(sequence)

    def test_target_allows_unknown_bases(self):
        assert normalize_target("acnug") == "ACNTG"

    def test_invalid_target_base_raises(self):
        with pytest.raises(ValueError, match="target sequence"):
            normalize_target("ACGTX")


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("AACG") == "CGTT"

    def test_rna_and_case_are_normalized(self):
        assert reverse_complement("augn") == "NCAT"

    def test_round_trip(self):
        sequence = "GCTANNGCTA"
        assert reverse_complement(reverse_complement(sequence)) == sequence

    def test_invalid_base_raises(self):
        with pytest.raises(ValueError, match="Invalid bases"):
            reverse_complement("ACGX")


class TestGuideConstruction:
    def test_sequence_and_pam_are_normalized(self, spcas9):
        guide = make_guide(
            spcas9,
            sequence="augc" * 5,
            pam="ugg",
        )
        assert guide.sequence == "ATGC" * 5
        assert guide.pam == "TGG"

    def test_invalid_sequence_raises(self, spcas9):
        with pytest.raises(ValueError, match="Invalid bases"):
            make_guide(spcas9, sequence="X" * 20)

    def test_wrong_spacer_length_raises(self, spcas9):
        with pytest.raises(ValueError, match="does not match"):
            make_guide(spcas9, sequence="A" * 19)

    def test_invalid_pam_raises(self, spcas9):
        with pytest.raises(ValueError, match="does not match"):
            make_guide(spcas9, pam="AAA")

    def test_pam_interval_length_must_match(self, spcas9):
        with pytest.raises(ValueError, match="interval length"):
            make_guide(spcas9, pam_interval=Interval(30, 34))

    def test_invalid_pam_strand_raises(self, spcas9):
        with pytest.raises(ValueError, match="pam_strand"):
            make_guide(spcas9, pam_strand="both")

    def test_frozen(self, spcas9):
        guide = make_guide(spcas9)
        with pytest.raises(AttributeError):
            guide.sequence = "T" * 20  # type: ignore[misc]

    @pytest.mark.parametrize(
        ("pam_strand", "target_strand"),
        [("+", "-"), ("-", "+")],
    )
    def test_target_strand_is_opposite_pam_strand(
        self,
        spcas9,
        pam_strand,
        target_strand,
    ):
        guide = make_guide(spcas9, pam_strand=pam_strand)
        assert guide.pam_strand == pam_strand
        assert guide.target_strand == target_strand

    def test_ambiguous_legacy_fields_are_absent(self, spcas9):
        guide = make_guide(spcas9)
        assert not hasattr(guide, "strand")
        assert not hasattr(guide, "pam_position")


class TestProtospacerInterval:
    @pytest.mark.parametrize(
        ("nuclease_name", "pam_interval", "pam_strand", "expected"),
        [
            ("SpCas9", Interval(20, 23), "+", Interval(0, 20)),
            ("SpCas9", Interval(10, 13), "-", Interval(13, 33)),
            ("AsCas12a", Interval(10, 14), "+", Interval(14, 37)),
            ("AsCas12a", Interval(30, 34), "-", Interval(7, 30)),
        ],
    )
    def test_interval_geometry(
        self,
        nuclease_name,
        pam_interval,
        pam_strand,
        expected,
    ):
        nuclease = get_nuclease(nuclease_name)
        guide = make_guide(
            nuclease,
            pam_interval=pam_interval,
            pam_strand=pam_strand,
        )
        assert guide.protospacer_interval == expected
        assert len(guide.protospacer_interval) == nuclease.spacer_length


class TestCutSite:
    @pytest.mark.parametrize(
        ("nuclease_name", "pam_interval", "pam_strand", "expected"),
        [
            ("SpCas9", Interval(20, 23), "+", 17),
            ("SpCas9", Interval(10, 13), "-", 16),
            ("AsCas12a", Interval(10, 14), "+", (32, 37)),
            ("AsCas12a", Interval(30, 34), "-", (7, 12)),
        ],
    )
    def test_cut_boundaries_use_forward_coordinates(
        self,
        nuclease_name,
        pam_interval,
        pam_strand,
        expected,
    ):
        guide = make_guide(
            get_nuclease(nuclease_name),
            pam_interval=pam_interval,
            pam_strand=pam_strand,
        )
        assert guide.cut_site() == expected


class TestFindGuides:
    spacer = "ATGCACGTTAGCTACGATCA"

    @classmethod
    def plus_cas9_target(cls) -> str:
        return "A" * 10 + cls.spacer + "TGG" + "A" * 10

    @classmethod
    def minus_cas9_target(cls) -> str:
        oriented_site = cls.spacer + "TGG"
        return "A" * 10 + reverse_complement(oriented_site) + "A" * 10

    def test_plus_pam_returns_construct_ready_sequence(self, spcas9):
        guides = find_guides(
            self.plus_cas9_target(),
            spcas9,
            pam_strand="+",
        )
        guide = next(g for g in guides if g.pam_interval == Interval(30, 33))

        assert guide.sequence == self.spacer
        assert guide.pam == "TGG"
        assert guide.pam_strand == "+"
        assert guide.target_strand == "-"
        assert guide.protospacer_interval == Interval(10, 30)

    def test_minus_pam_maps_back_to_forward_coordinates(self, spcas9):
        guides = find_guides(
            self.minus_cas9_target(),
            spcas9,
            pam_strand="-",
        )
        guide = next(g for g in guides if g.pam_interval == Interval(10, 13))

        assert guide.sequence == self.spacer
        assert guide.pam == "TGG"
        assert guide.pam_strand == "-"
        assert guide.target_strand == "+"
        assert guide.protospacer_interval == Interval(13, 33)

    def test_reverse_complement_symmetry(self, spcas9):
        target = self.plus_cas9_target()
        reverse_target = reverse_complement(target)
        plus = next(
            g
            for g in find_guides(target, spcas9, pam_strand="+")
            if g.pam_interval == Interval(30, 33)
        )
        minus = next(
            g
            for g in find_guides(reverse_target, spcas9, pam_strand="-")
            if g.pam_interval == Interval(10, 13)
        )

        assert plus.sequence == minus.sequence
        assert plus.pam == minus.pam
        assert plus.pam_strand != minus.pam_strand
        assert plus.pam_interval == Interval(
            len(target) - minus.pam_interval.end,
            len(target) - minus.pam_interval.start,
        )

    def test_cas12a_uses_five_prime_pam_geometry(self, ascas12a):
        spacer = "ACGTACGTACGTACGTACGTACG"
        target = "C" * 10 + "TTTG" + spacer + "C" * 10
        guide = next(
            g
            for g in find_guides(target, ascas12a, pam_strand="+")
            if g.pam_interval == Interval(10, 14)
        )

        assert guide.sequence == spacer
        assert guide.protospacer_interval == Interval(14, 37)

    @pytest.mark.parametrize("pam_strand", ["+", "-", "both"])
    def test_pam_strand_filter(self, spcas9, pam_strand):
        guides = find_guides(
            self.plus_cas9_target(),
            spcas9,
            pam_strand=pam_strand,
        )
        expected = {pam_strand} if pam_strand != "both" else {"+", "-"}
        assert all(guide.pam_strand in expected for guide in guides)

    def test_invalid_pam_strand_filter_raises(self, spcas9):
        with pytest.raises(ValueError, match="pam_strand"):
            find_guides(
                self.plus_cas9_target(),
                spcas9,
                pam_strand="invalid",  # type: ignore[arg-type]
            )

    def test_no_pam_returns_empty(self, spcas9):
        assert find_guides("A" * 50, spcas9, pam_strand="+") == []

    def test_unknown_reference_bases_are_skipped(self, spcas9):
        guides = find_guides("N" * 20 + "AGG", spcas9, pam_strand="+")
        assert guides == []

    def test_invalid_target_base_raises(self, spcas9):
        with pytest.raises(ValueError, match="target sequence"):
            find_guides("A" * 20 + "XGG", spcas9)

    def test_reported_sequences_reconstruct_from_reference(self, spcas9):
        target = self.plus_cas9_target() + self.minus_cas9_target()
        guides = find_guides(target, spcas9)
        assert guides

        for guide in guides:
            protospacer = target[
                guide.protospacer_interval.start : guide.protospacer_interval.end
            ]
            pam = target[guide.pam_interval.start : guide.pam_interval.end]
            if guide.pam_strand == "-":
                protospacer = reverse_complement(protospacer)
                pam = reverse_complement(pam)
            assert guide.sequence == protospacer
            assert guide.pam == pam


class TestValidateGuideSequence:
    @pytest.mark.parametrize("length", [17, 20, 23, 30])
    def test_valid_lengths(self, length):
        validate_guide_sequence("A" * length)

    @pytest.mark.parametrize("length", [0, 16, 31, 35])
    def test_invalid_lengths_raise(self, length):
        with pytest.raises(ValueError):
            validate_guide_sequence("A" * length)

    def test_invalid_alphabet_raises(self):
        with pytest.raises(ValueError):
            validate_guide_sequence("X" * 20)

    def test_rna_input_is_accepted(self):
        validate_guide_sequence("AUGCAUGCAUGCAUGCAUGC")


class TestGuideGCContent:
    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("G" * 20, 1.0),
            ("A" * 20, 0.0),
            ("GCGCGCGCGC" + "A" * 10, 0.5),
        ],
    )
    def test_gc_fraction(self, spcas9, sequence, expected):
        guide = make_guide(spcas9, sequence=sequence)
        assert guide_gc_content(guide) == pytest.approx(expected)
        assert isinstance(guide_gc_content(guide), float)
