import pytest

from geneops.nuclease import (
    Nuclease,
    RRegistry,
    get_nuclease,
    list_nucleases,
    register_nuclease,
    matches_pam,
    find_pam_sites,
    pam_orientation,
)

# --- Registry tests ---

def test_list_nucleases_contains_defaults():
    names = list_nucleases()
    assert "SpCas9" in names
    assert "SaCas9" in names
    assert "AsCas12a" in names
    assert "nCas9" in names

def test_get_nuclease_case_insensitive():
    sp = get_nuclease("spcas9")
    assert sp.name == "SpCas9"
    assert sp.pam == "NGG"

def test_register_nuclease_adds_to_registry():
    registry = RRegistry()
    registry.register_nuclease(Nuclease(name="AsCas12a", pam="TTTV"))
    assert registry.get_nuclease("AsCas12a").pam == "TTTV"

def test_register_duplicate_raises():
    registry = RRegistry([Nuclease(name="SpCas9", pam="NGG")])
    with pytest.raises(ValueError):
        registry.register_nuclease(Nuclease(name="spcas9", pam="NGG"))

def test_get_unknown_raises():
    registry = RRegistry()
    with pytest.raises(KeyError):
        registry.get_nuclease("Unknown")

# --- PAM matching tests ---

def test_pam_exact_match():
    sp = get_nuclease("SpCas9")
    assert sp.matches_pam("AGG")
    assert sp.matches_pam("TGG")
    assert not sp.matches_pam("AGA")

def test_pam_with_ambiguity_codes():
    spry = get_nuclease("SpRY")  # NRN
    assert spry.matches_pam("AGA")  # A-G-A
    assert spry.matches_pam("GAG")  # G-A-G
    assert spry.matches_pam("TGT")  # T-G-T
    assert not spry.matches_pam("GCG")  # C is not R

def test_cas12a_tt_tv_pam():
    cas12a = get_nuclease("AsCas12a")  # TTTV
    assert cas12a.matches_pam("TTTA")
    assert cas12a.matches_pam("TTTC")
    assert cas12a.matches_pam("TTTG")
    assert not cas12a.matches_pam("TTTT")

def test_pam_length_mismatch():
    sp = get_nuclease("SpCas9")
    assert not sp.matches_pam("GG")
    assert not sp.matches_pam("AGGG")

# --- Integration-style example ---

def find_nuclease_pams(seq: str, nuclease: Nuclease):
    hits = []
    pam_len = nuclease.pam_length()

    for i in range(len(seq) - pam_len + 1):
        if nuclease.matches_pam(seq[i : i + pam_len]):
            hits.append(i)

    return hits

def test_find_pams_with_custom_nuclease():
    seq = "ATGCTGACCTGACCGGATCGTACGTTACGATCGGAGCTTAGGCTAGCTGATCG"

    custom = Nuclease(
        name="MyCas9",
        pam="NGA",
        description="Custom engineered Cas9 variant",
    )
    register_nuclease(custom)

    mycas9 = get_nuclease("mycas9")
    sites = find_nuclease_pams(seq, mycas9)

    # All returned PAMs must truly match NGA
    for pos in sites:
        pam_seq = seq[pos : pos + 3]
        assert mycas9.matches_pam(pam_seq)

# --- PAM logic function tests ---

def test_matches_pam_basic():
    """Test basic PAM matching with exact sequences."""
    assert matches_pam("AGG", "NGG")
    assert matches_pam("TGG", "NGG")
    assert matches_pam("CGG", "NGG")
    assert matches_pam("GGG", "NGG")
    assert not matches_pam("AGC", "NGG")
    assert not matches_pam("AGA", "NGG")

def test_matches_pam_iupac_codes():
    """Test PAM matching with IUPAC ambiguity codes."""
    # R = A or G
    assert matches_pam("AGA", "NRN")
    assert matches_pam("GGT", "NRN")
    assert not matches_pam("GCT", "NRN")  # C is not R
    
    # V = A, C, or G (not T)
    assert matches_pam("TTTA", "TTTV")
    assert matches_pam("TTTC", "TTTV")
    assert matches_pam("TTTG", "TTTV")
    assert not matches_pam("TTTT", "TTTV")
    
    # W = A or T
    assert matches_pam("AAA", "NWN")  # A is W
    assert matches_pam("TAT", "NWN")  # T is W
    assert not matches_pam("GCG", "NWN")  # C is not W

def test_matches_pam_case_insensitive():
    """Test that PAM matching is case-insensitive."""
    assert matches_pam("agg", "NGG")
    assert matches_pam("AGG", "ngg")
    assert matches_pam("aGg", "NgG")

def test_matches_pam_length_mismatch():
    """Test that mismatched lengths return False."""
    assert not matches_pam("AG", "NGG")
    assert not matches_pam("AGGG", "NGG")
    assert not matches_pam("", "NGG")

def test_find_pam_sites_single_pam():
    """Test finding PAM sites with a single match."""
    seq = "ATCAGGTT"
    sites = find_pam_sites(seq, "NGG")
    assert sites == [3]  # AGG at position 3

def test_find_pam_sites_multiple_pams():
    """Test finding multiple PAM sites in a sequence."""
    seq = "ATCGGTAGGATCCGG"
    sites = find_pam_sites(seq, "NGG")
    # Sequence: A T C G G T A G G A T C C G G
    # Index:    0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
    # CGG at index 2, AGG at index 6, CCG at index 11, CGG at index 12
    assert sites == [2, 6, 12]  # CGG at 2, AGG at 6, CGG at 12

def test_find_pam_sites_cas12a():
    """Test finding Cas12a PAM sites (TTTV)."""
    seq = "TTTGATCGTTTAACGTTTC"
    sites = find_pam_sites(seq, "TTTV")
    assert sites == [0, 8, 15]  # TTTG at 0, TTTA at 8, TTTC at 15

def test_find_pam_sites_no_matches():
    """Test finding PAM sites when none exist."""
    seq = "AAAAAAAAAA"
    sites = find_pam_sites(seq, "NGG")
    assert sites == []

def test_find_pam_sites_overlapping():
    """Test finding overlapping PAM sites."""
    seq = "AGGGG"
    sites = find_pam_sites(seq, "NGG")
    # AGG at 0, GGG at 1, GGG at 2
    assert sites == [0, 1, 2]

def test_find_pam_sites_empty_sequence():
    """Test finding PAM sites in an empty sequence."""
    sites = find_pam_sites("", "NGG")
    assert sites == []

def test_find_pam_sites_sequence_shorter_than_pam():
    """Test finding PAM sites when sequence is shorter than PAM."""
    sites = find_pam_sites("AG", "NGG")
    assert sites == []

def test_pam_orientation_cas9():
    """Test that Cas9 variants have 3' PAM."""
    assert pam_orientation("SpCas9") == "3'"
    assert pam_orientation("SaCas9") == "3'"
    assert pam_orientation("SpCas9-NG") == "3'"
    assert pam_orientation("SpG") == "3'"
    assert pam_orientation("SpRY") == "3'"
    assert pam_orientation("nCas9") == "3'"

def test_pam_orientation_cas12a():
    """Test that Cas12a variants have 5' PAM."""
    assert pam_orientation("AsCas12a") == "5'"
    assert pam_orientation("LbCas12a") == "5'"

def test_pam_orientation_with_nuclease_object():
    """Test pam_orientation with a Nuclease object instead of string."""
    cas9 = get_nuclease("SpCas9")
    assert pam_orientation(cas9) == "3'"
    
    cas12a = get_nuclease("AsCas12a")
    assert pam_orientation(cas12a) == "5'"

def test_pam_orientation_custom_nuclease():
    """Test pam_orientation with custom nuclease objects."""
    # Custom Cas9-like nuclease
    custom_cas9 = Nuclease(name="CustomCas9", pam="NGG")
    assert pam_orientation(custom_cas9) == "3'"
    
    # Custom Cas12-like nuclease
    custom_cas12 = Nuclease(name="CustomCas12a", pam="TTTV")
    assert pam_orientation(custom_cas12) == "5'"