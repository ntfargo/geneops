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

"""Tests for cut geometry functions."""

import pytest
from geneops.nuclease import (
    CleavagePattern,
    Nuclease,
    PAM,
    cut_offset,
    produces_blunt_cut,
    overhang_length,
    calculate_cut_site,
    cut_position,
    get_nuclease,
)

class TestCutOffset:
    """Test cut_offset function."""
    
    def test_spcas9_offset(self):
        """SpCas9 cuts 3bp upstream of PAM."""
        assert cut_offset("SpCas9") == -3
    
    def test_cas12a_offset(self):
        """Cas12a cuts downstream of PAM."""
        assert cut_offset("AsCas12a") == 18
        assert cut_offset("LbCas12a") == 18
    
    def test_with_nuclease_object(self):
        """Should work with Nuclease objects."""
        nuclease = get_nuclease("SpCas9")
        assert cut_offset(nuclease) == -3

    def test_name_does_not_select_cut_behavior(self):
        nuclease = Nuclease(
            name="Cas12InNameOnly",
            pam=PAM("NGG", "3'"),
            cleavage=CleavagePattern(17, 17),
        )
        assert cut_offset(nuclease) == -3

class TestProducesBluntCut:
    """Test produces_blunt_cut function."""
    
    def test_cas9_blunt(self):
        """Cas9 family produces blunt cuts."""
        assert produces_blunt_cut("SpCas9") is True
        assert produces_blunt_cut("SaCas9") is True
        assert produces_blunt_cut("SpCas9-HF1") is True
    
    def test_cas12a_staggered(self):
        """Cas12a produces staggered cuts."""
        assert produces_blunt_cut("AsCas12a") is False
        assert produces_blunt_cut("LbCas12a") is False

    def test_explicit_custom_geometry(self):
        nuclease = Nuclease(
            name="Aurora",
            pam=PAM("TTTV", "5'"),
            cleavage=CleavagePattern(18, 23),
            spacer_length=23,
        )
        assert produces_blunt_cut(nuclease) is False

class TestOverhangLength:
    """Test overhang_length function."""
    
    def test_cas9_no_overhang(self):
        """Cas9 has no overhang (blunt cut)."""
        assert overhang_length("SpCas9") == 0
        assert overhang_length("SaCas9") == 0
    
    def test_cas12a_overhang(self):
        """Cas12a has 5nt overhang."""
        assert overhang_length("AsCas12a") == 5
        assert overhang_length("LbCas12a") == 5

    def test_nickase_has_no_double_strand_overhang(self):
        with pytest.raises(ValueError, match="nickase"):
            overhang_length("nCas9")

class TestCutPosition:
    """Test cut_position function."""
    
    def test_spcas9_forward_strand(self):
        """SpCas9 on forward strand."""
        # PAM at position 20, cuts 3bp upstream at position 17
        result = cut_position(20, "SpCas9", "+")
        assert result == 17
    
    def test_spcas9_returns_single_position(self):
        """SpCas9 returns single int (blunt cut)."""
        result = cut_position(20, "SpCas9", "+")
        assert isinstance(result, int)
    
    def test_cas12a_forward_strand(self):
        """Cas12a on forward strand produces staggered cut."""
        # PAM at position 10 (TTTV = 4bp), cuts downstream at ~18bp offset
        result = cut_position(10, "AsCas12a", "+")
        assert isinstance(result, tuple)
        assert len(result) == 2
        # Should be (32, 37) - cut at 10+4+18=32, with 5nt overhang
        assert result == (32, 37)
    
    def test_cas12a_returns_tuple(self):
        """Cas12a returns tuple (staggered cut)."""
        result = cut_position(10, "AsCas12a", "+")
        assert isinstance(result, tuple)
        top, bottom = result
        assert bottom - top == 5  # 5nt overhang

    def test_reverse_strand_uses_explicit_offsets(self):
        assert cut_position(30, "AsCas12a", "-") == (7, 12)

    def test_custom_name_has_no_effect_on_cut_position(self):
        nuclease = Nuclease(
            name="NotCas12DespiteTheName",
            pam=PAM("NGG", "3'"),
            cleavage=CleavagePattern(17, 17),
        )
        assert cut_position(20, nuclease, "+") == 17
    
    def test_different_pam_positions(self):
        """Test with various PAM positions."""
        # Position 0
        assert cut_position(0, "SpCas9", "+") == -3
        
        # Position 100
        assert cut_position(100, "SpCas9", "+") == 97
        
        # Position 50 with Cas12a
        top, bottom = cut_position(50, "AsCas12a", "+")
        assert bottom - top == 5


class TestCalculateCutSite:
    def test_preserves_strand_identity_for_forward_cas12a(self):
        site = calculate_cut_site(10, "AsCas12a", "+")

        assert site.pam_strand_boundary == 32
        assert site.target_strand_boundary == 37
        assert site.plus_strand_boundary == 32
        assert site.minus_strand_boundary == 37
        assert site.boundaries == (32, 37)
        assert site.overhang_length == 5
        assert site.overhang == "5'"

    def test_maps_strand_identity_for_reverse_cas12a(self):
        site = calculate_cut_site(30, "AsCas12a", "-")

        assert site.pam_strand_boundary == 12
        assert site.target_strand_boundary == 7
        assert site.plus_strand_boundary == 7
        assert site.minus_strand_boundary == 12
        assert site.boundaries == (7, 12)
        assert site.overhang == "5'"

    def test_nickase_identifies_cut_reference_strand(self):
        site = calculate_cut_site(20, "nCas9", "+")

        assert site.is_nickase
        assert site.pam_strand_boundary is None
        assert site.target_strand_boundary == 17
        assert site.nicking_strand == "-"
        assert site.overhang_length is None

if __name__ == "__main__":
    pytest.main([__file__, "-v"])