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

from geneops.guide import (
    Guide,
    normalize,
    reverse_complement,
    find_guides,
    validate_guide_sequence,
    guide_gc_content,
)
from geneops.nuclease import (
    Nuclease,
    PAM,
    get_nuclease,
    pam_orientation,
    cut_position,
)

class TestNormalize:
    def test_uppercase(self):
        assert normalize("acgt") == "ACGT"

    def test_u_to_t(self):
        assert normalize("AUGC") == "ATGC"

    def test_mixed_case_and_u(self):
        assert normalize("aUgC") == "ATGC"

    def test_already_valid(self):
        assert normalize("ACGT") == "ACGT"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            normalize("")

    def test_invalid_base_raises(self):
        with pytest.raises(ValueError, match="Invalid bases"):
            normalize("ACGX")

    def test_multiple_invalid_bases(self):
        with pytest.raises(ValueError, match="Invalid bases"):
            normalize("XYZACGT")

    def test_whitespace_rejected(self):
        with pytest.raises(ValueError):
            normalize("AC GT")

    def test_digit_rejected(self):
        with pytest.raises(ValueError):
            normalize("ACG1T")

class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ACGT") == "ACGT"  # palindrome

    def test_non_palindrome(self):
        assert reverse_complement("AAAA") == "TTTT"

    def test_single_base(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("G") == "C"

    def test_asymmetric(self):
        assert reverse_complement("AACG") == "CGTT"

    def test_round_trip(self):
        seq = "GCTAGCTA"
        assert reverse_complement(reverse_complement(seq)) == seq

@pytest.fixture
def spcas9() -> Nuclease:
    return get_nuclease("SpCas9")

@pytest.fixture
def ascas12a() -> Nuclease:
    return get_nuclease("AsCas12a")

class TestGuideConstruction:
    def test_sequence_normalized_on_init(self, spcas9):
        g = Guide(
            sequence="acgt" * 5,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert g.sequence == "ACGT" * 5

    def test_rna_input_converted(self, spcas9):
        g = Guide(
            sequence="AUGCAUGCAUGCAUGCAUGC",
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert "U" not in g.sequence
        assert g.sequence == "ATGCATGCATGCATGCATGC"

    def test_invalid_sequence_raises(self, spcas9):
        with pytest.raises(ValueError):
            Guide(
                sequence="XXXX",
                nuclease=spcas9,
                pam="AGG",
                pam_position=20,
                strand="+",
            )

    def test_frozen(self, spcas9):
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        with pytest.raises(AttributeError):
            g.sequence = "T" * 20  # type: ignore[misc]

    def test_strand_values(self, spcas9):
        for s in ("+", "-"):
            g = Guide(
                sequence="A" * 20,
                nuclease=spcas9,
                pam="AGG",
                pam_position=20,
                strand=s,
            )
            assert g.strand == s

class TestBindingInterval:
    """binding_interval returns (start, end) of spacer on target."""

    def test_cas9_forward_strand(self, spcas9):
        # SpCas9: 3' PAM, + strand → spacer is upstream (left) of PAM
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        start, end = g.binding_interval()
        assert start == 0
        assert end == 20
        assert end - start == spcas9.spacer_length

    def test_cas9_reverse_strand(self, spcas9):
        # SpCas9: 3' PAM, - strand → orientation "3'" with "-" means
        # spacer downstream of PAM in fwd coords
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=10,
            strand="-",
        )
        start, end = g.binding_interval()
        assert end - start == spcas9.spacer_length
        # 5' orientation rule: spacer to the right of PAM
        assert start == 10 + 3  # pam_position + pam_len
        assert end == 10 + 3 + 20

    def test_cas12a_forward_strand(self, ascas12a):
        # AsCas12a: 5' PAM, + strand → spacer downstream (right) of PAM
        g = Guide(
            sequence="A" * 23,
            nuclease=ascas12a,
            pam="TTTG",
            pam_position=10,
            strand="+",
        )
        start, end = g.binding_interval()
        assert start == 10 + 4  # pam_pos + pam_len(4)
        assert end == 10 + 4 + 23

    def test_cas12a_reverse_strand(self, ascas12a):
        # AsCas12a: 5' PAM, - strand → spacer upstream (left) of PAM
        g = Guide(
            sequence="A" * 23,
            nuclease=ascas12a,
            pam="TTTG",
            pam_position=30,
            strand="-",
        )
        start, end = g.binding_interval()
        assert end - start == ascas12a.spacer_length
        assert end == 30  # spacer ends at pam_position
        assert start == 30 - 23

    def test_interval_length_matches_spacer(self, spcas9, ascas12a):
        for nuc in (spcas9, ascas12a):
            pam_obj = nuc.get_pam()
            g = Guide(
                sequence="A" * nuc.spacer_length,
                nuclease=nuc,
                pam="A" * len(pam_obj),
                pam_position=30,
                strand="+",
            )
            s, e = g.binding_interval()
            assert e - s == nuc.spacer_length

class TestCutSite:
    def test_cas9_blunt_cut_returns_int(self, spcas9):
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        result = g.cut_site()
        assert isinstance(result, int)

    def test_cas12a_staggered_cut_returns_tuple(self, ascas12a):
        g = Guide(
            sequence="A" * 23,
            nuclease=ascas12a,
            pam="TTTG",
            pam_position=10,
            strand="+",
        )
        result = g.cut_site()
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_cut_site_delegates_to_nuclease(self, spcas9):
        pam_pos = 25
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=pam_pos,
            strand="+",
        )
        expected = cut_position(pam_pos, spcas9, "+")
        assert g.cut_site() == expected

    def test_cut_site_reverse_strand(self, spcas9):
        pam_pos = 15
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=pam_pos,
            strand="-",
        )
        expected = cut_position(pam_pos, spcas9, "-")
        assert g.cut_site() == expected

class TestBindsStrand:
    def test_forward(self, spcas9):
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert g.binds_strand() == "+"

    def test_reverse(self, spcas9):
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="-",
        )
        assert g.binds_strand() == "-"
class TestFindGuides:
    """find_guides scans a target for PAM sites and returns Guide objects."""

    # Helper: build a target with a known PAM at a known position for SpCas9
    # Target layout (+ strand, 3' PAM): [spacer][PAM]...
    @staticmethod
    def _cas9_target() -> str:
        spacer = "ATGCATGCATGCATGCATGC"  # 20 nt
        pam = "TGG"
        flank = "A" * 10
        return flank + spacer + pam + flank

    def test_returns_list(self, spcas9):
        guides = find_guides("AAAAAAAAAAAAAAAAAAAAAAATGGAAAA", spcas9)
        assert isinstance(guides, list)

    def test_finds_cas9_pam_on_plus(self, spcas9):
        target = self._cas9_target()
        guides = find_guides(target, spcas9, strand="+")
        assert len(guides) >= 1
        for g in guides:
            assert g.strand == "+"
            assert isinstance(g, Guide)

    def test_guide_nuclease_matches(self, spcas9):
        target = self._cas9_target()
        guides = find_guides(target, spcas9)
        for g in guides:
            assert g.nuclease is spcas9

    def test_guide_sequence_length(self, spcas9):
        target = self._cas9_target()
        guides = find_guides(target, spcas9)
        for g in guides:
            assert len(g.sequence) == spcas9.spacer_length

    def test_strand_filter_plus(self, spcas9):
        target = self._cas9_target()
        guides = find_guides(target, spcas9, strand="+")
        assert all(g.strand == "+" for g in guides)

    def test_strand_filter_minus(self, spcas9):
        target = self._cas9_target()
        guides = find_guides(target, spcas9, strand="-")
        assert all(g.strand == "-" for g in guides)

    def test_both_strands(self, spcas9):
        target = self._cas9_target()
        guides = find_guides(target, spcas9, strand="both")
        strands = {g.strand for g in guides}
        # Should be able to find guides on at least one strand
        assert len(guides) >= 1

    def test_no_pam_returns_empty(self, spcas9):
        # A target with no NGG at all
        target = "A" * 50
        guides = find_guides(target, spcas9, strand="+")
        assert guides == []

    def test_cas12a_finds_guides(self, ascas12a):
        # Target layout (+ strand, 5' PAM): [PAM][spacer]...
        pam = "TTTG"
        spacer = "A" * 23
        flank = "C" * 10
        target = flank + pam + spacer + flank
        guides = find_guides(target, ascas12a, strand="+")
        assert len(guides) >= 1
        for g in guides:
            assert len(g.sequence) == ascas12a.spacer_length

    def test_rna_target_handled(self, spcas9):
        # Target with U instead of T should still work (normalize_target)
        target = self._cas9_target().replace("T", "U")
        guides = find_guides(target, spcas9)
        # Should find the same guides as DNA target
        assert isinstance(guides, list)

    def test_pam_position_is_nonnegative(self, spcas9):
        target = self._cas9_target()
        guides = find_guides(target, spcas9)
        for g in guides:
            assert g.pam_position >= 0

    def test_guide_pam_field_matches_target(self, spcas9):
        """The .pam field of each returned Guide should actually match the
        nuclease PAM pattern."""
        target = self._cas9_target()
        guides = find_guides(target, spcas9, strand="+")
        for g in guides:
            assert spcas9.matches_pam(g.pam)

class TestValidateGuideSequence:
    def test_valid_20mer(self):
        validate_guide_sequence("A" * 20)  # should not raise

    def test_valid_23mer(self):
        validate_guide_sequence("A" * 23)  # within 17-30

    def test_too_short_raises(self):
        with pytest.raises(ValueError, match="outside the expected range"):
            validate_guide_sequence("A" * 10)

    def test_too_long_raises(self):
        with pytest.raises(ValueError, match="outside the expected range"):
            validate_guide_sequence("A" * 35)

    def test_boundary_17(self):
        validate_guide_sequence("G" * 17)  # min accepted

    def test_boundary_30(self):
        validate_guide_sequence("C" * 30)  # max accepted

    def test_boundary_16_raises(self):
        with pytest.raises(ValueError):
            validate_guide_sequence("G" * 16)

    def test_boundary_31_raises(self):
        with pytest.raises(ValueError):
            validate_guide_sequence("C" * 31)

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            validate_guide_sequence("")

    def test_invalid_alphabet_raises(self):
        with pytest.raises(ValueError):
            validate_guide_sequence("X" * 20)

    def test_rna_input_accepted(self):
        # U→T normalization happens first, so RNA is fine
        validate_guide_sequence("AUGCAUGCAUGCAUGCAUGC")

class TestGuideGCContent:
    def test_all_gc(self, spcas9):
        g = Guide(
            sequence="G" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert guide_gc_content(g) == pytest.approx(1.0)

    def test_no_gc(self, spcas9):
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert guide_gc_content(g) == pytest.approx(0.0)

    def test_half_gc(self, spcas9):
        g = Guide(
            sequence="GCGCGCGCGCAAAAAAAAAA",
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert guide_gc_content(g) == pytest.approx(0.5)

    def test_mixed(self, spcas9):
        seq = "ATGCATGCATGCATGCATGC"  # 10 G/C out of 20
        g = Guide(
            sequence=seq,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert guide_gc_content(g) == pytest.approx(0.5)

    def test_returns_float(self, spcas9):
        g = Guide(
            sequence="A" * 20,
            nuclease=spcas9,
            pam="AGG",
            pam_position=20,
            strand="+",
        )
        assert isinstance(guide_gc_content(g), float)

class TestIntegration:
    """Verify that guides discovered by find_guides have coherent attributes."""

    def test_cas9_guide_cut_inside_target(self, spcas9):
        spacer = "ATGCATGCATGCATGCATGC"
        pam = "TGG"
        flank = "A" * 30
        target = flank + spacer + pam + flank
        guides = find_guides(target, spcas9, strand="+")
        for g in guides:
            cut = g.cut_site()
            assert isinstance(cut, int)
            assert 0 <= cut <= len(target)

    def test_cas12a_guide_cut_inside_target(self, ascas12a):
        pam = "TTTG"
        spacer = "ACGTACGTACGTACGTACGTACG"  # 23 nt
        flank = "C" * 30
        target = flank + pam + spacer + flank
        guides = find_guides(target, ascas12a, strand="+")
        for g in guides:
            cut = g.cut_site()
            assert isinstance(cut, tuple)
            assert all(0 <= c <= len(target) for c in cut)

    def test_binding_interval_inside_target(self, spcas9):
        spacer = "ATGCATGCATGCATGCATGC"
        pam = "TGG"
        flank = "A" * 30
        target = flank + spacer + pam + flank
        guides = find_guides(target, spcas9, strand="+")
        for g in guides:
            s, e = g.binding_interval()
            assert s >= 0
            assert e <= len(target)
            assert e - s == spcas9.spacer_length