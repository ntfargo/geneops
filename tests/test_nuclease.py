import pytest

from geneops.nuclease import (
    Nuclease,
    PAM,
    RRegistry,
    get_nuclease,
    list_nucleases,
    register_nuclease,
    matches_pam,
    find_pam_sites,
    _reverse_complement,
    pam_orientation,
    guide_length,
    is_guide_length_valid,
    is_guide_compatible,
    recognizes_strand,
    find_pam_sites_with_strand,
    binding_strand,
    validate_nuclease,
    is_nickase,
    nicking_offset,
)

@pytest.fixture
def spcas9():
    return get_nuclease("SpCas9")

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

def test_guide_length_cas9():
    """Test that SpCas9 variants require 20nt guides."""
    assert guide_length("SpCas9") == 20
    assert guide_length("nCas9") == 20
    assert guide_length("SpCas9-NG") == 20

def test_guide_length_cas9_orthologs():
    """Test that Cas9 orthologs have their specific spacer lengths."""
    assert guide_length("SaCas9") == 21  # S. aureus Cas9
    assert guide_length("NmCas9") == 24  # N. meningitidis Cas9
    assert guide_length("CjCas9") == 22  # C. jejuni Cas9

def test_guide_length_cas12a():
    """Test that Cas12a variants require 23nt guides."""
    assert guide_length("AsCas12a") == 23
    assert guide_length("LbCas12a") == 23

def test_guide_length_with_nuclease_object():
    """Test guide_length with Nuclease object instead of string."""
    cas9 = get_nuclease("SpCas9")
    assert guide_length(cas9) == 20
    
    cas12a = get_nuclease("AsCas12a")
    assert guide_length(cas12a) == 23

def test_guide_length_custom_nuclease():
    """Test guide_length with custom nuclease with explicit spacer_length."""
    custom = Nuclease(name="CustomNuclease", pam="NGG", spacer_length=18)
    assert guide_length(custom) == 18

def test_is_guide_length_valid_cas9():
    """Test guide length validation for Cas9."""
    valid_guide = "A" * 20
    short_guide = "A" * 19
    long_guide = "A" * 21
    
    assert is_guide_length_valid(valid_guide, "SpCas9")
    assert not is_guide_length_valid(short_guide, "SpCas9")
    assert not is_guide_length_valid(long_guide, "SpCas9")

def test_is_guide_length_valid_cas12a():
    """Test guide length validation for Cas12a."""
    valid_guide = "A" * 23
    short_guide = "A" * 22
    long_guide = "A" * 24
    
    assert is_guide_length_valid(valid_guide, "AsCas12a")
    assert not is_guide_length_valid(short_guide, "AsCas12a")
    assert not is_guide_length_valid(long_guide, "AsCas12a")

def test_is_guide_length_valid_empty_guide():
    """Test guide length validation with empty guide."""
    assert not is_guide_length_valid("", "SpCas9")
    assert not is_guide_length_valid("", "AsCas12a")

def test_is_guide_compatible_cas9_valid():
    """Test guide compatibility for Cas9 with valid guide and target."""
    # Geneops reports the construct-ready guide in PAM-strand orientation.
    guide = "CGATCGATCGATCGATCGAT"
    target = "CGATCGATCGATCGATCGATAGG"
    
    assert is_guide_compatible(guide, target, "SpCas9")

def test_is_guide_compatible_cas9_wrong_length():
    """Test guide compatibility fails with wrong guide length."""
    guide = "CGATCGATCGATCGATCGA"  # 19nt (wrong length)
    target = "CGATCGATCGATCGATCGATAGG"
    
    assert not is_guide_compatible(guide, target, "SpCas9")

def test_is_guide_compatible_cas9_no_pam():
    """Test guide compatibility fails without valid PAM."""
    guide = "CGATCGATCGATCGATCGAT"  # 20nt
    # Target has no valid PAM (NGG)
    target = "CGATCGATCGATCGATCGATAAT"  # No NGG
    
    assert not is_guide_compatible(guide, target, "SpCas9")

def test_is_guide_compatible_cas9_mismatched_guide():
    """Test guide compatibility fails when guide doesn't match binding site."""
    guide = "AAAAAAAAAAAAAAAAAAAA"  # 20nt guide
    target = "CGATCGATCGATCGATCGATAGG"  # Different binding site
    
    assert not is_guide_compatible(guide, target, "SpCas9")

def test_is_guide_compatible_cas12a_valid():
    """Test guide compatibility for Cas12a with valid guide and target."""
    # Cas12a has 5' PAM (TTTV), guide binds downstream of PAM
    # Target: [TTTV][23nt binding site]
    guide = "GATCGATCGATCGATCGATCGAT"  # 23nt guide
    target = "TTTAGATCGATCGATCGATCGATCGAT"  # TTTA + 23nt
    
    assert is_guide_compatible(guide, target, "AsCas12a")

def test_is_guide_compatible_cas12a_wrong_length():
    """Test guide compatibility fails for Cas12a with wrong guide length."""
    guide = "ATCGATCGATCGATCGATCG"  # 20nt (wrong for Cas12a)
    target = "TTTAGATCGATCGATCGATCGATCGAT"
    
    assert not is_guide_compatible(guide, target, "AsCas12a")

def test_is_guide_compatible_cas12a_no_pam():
    """Test guide compatibility fails for Cas12a without valid PAM."""
    guide = "ATCGATCGATCGATCGATCGATC"  # 23nt
    # TTTT is not a valid Cas12a PAM (V = A, C, G, not T)
    target = "TTTTATCGATCGATCGATCGATCGAT"
    
    assert not is_guide_compatible(guide, target, "AsCas12a")

def test_is_guide_compatible_case_insensitive():
    """Test guide compatibility is case-insensitive."""
    guide = "cgatcgatcgatcgatcgat"  # lowercase
    target = "CGATCGATCGATCGATCGATAGG"  # uppercase
    
    assert is_guide_compatible(guide, target, "SpCas9")

def test_is_guide_compatible_with_nuclease_object():
    """Test guide compatibility with Nuclease object instead of string."""
    cas9 = get_nuclease("SpCas9")
    guide = "CGATCGATCGATCGATCGAT"
    target = "CGATCGATCGATCGATCGATAGG"
    
    assert is_guide_compatible(guide, target, cas9)

def test_find_pam_sites_forward_strand(spcas9):
    seq = "AAACCCGGGTTT"  # PAM = GGG at position 6
    pam = spcas9.get_pam()

    sites = list(find_pam_sites_with_strand(seq, pam))

    assert (6, "+") in sites
    assert all(pos >= 0 for pos, _ in sites)

def test_find_pam_sites_reverse_strand_coordinate_mapping(spcas9):
    # PAM exists only on reverse complement
    seq = "CCCAAATTTCCC"
    pam = spcas9.get_pam()

    sites = list(find_pam_sites_with_strand(seq, pam))

    # Ensure reverse-strand PAM is reported in original coordinates
    assert any(strand == "-" for _, strand in sites)
    for pos, strand in sites:
        assert 0 <= pos < len(seq)

def test_binding_strand_pam_on_plus(spcas9):
    """
    PAM on + strand → guide binds − strand
    """
    target = (
        "AAAAA"
        "TTTTTTTTTTTTTTTTTTTT"  # protospacer (20)
        "AGG"                   # PAM
        "AAAAA"
    )

    guide = "TTTTTTTTTTTTTTTTTTTT"

    strand = binding_strand(guide, target, spcas9)

    assert strand == "-"


def test_binding_strand_pam_on_minus(spcas9):
    """A PAM on the reverse strand means the guide binds the plus strand."""
    guide = "ATGCACGTTAGCTACGATCA"
    reverse_oriented_site = _reverse_complement(guide + "TGG")
    target = "A" * 10 + reverse_oriented_site + "A" * 10

    assert binding_strand(guide, target, spcas9) == "+"

def test_validate_nuclease_accepts_valid_cas9():
    """Valid default nucleases pass validation without errors."""
    cas9 = get_nuclease("SpCas9")
    validate_nuclease(cas9)  # should not raise

def test_validate_nuclease_accepts_valid_cas12a():
    cas12a = get_nuclease("AsCas12a")
    validate_nuclease(cas12a)  # should not raise

def test_validate_nuclease_accepts_custom_nuclease():
    custom = Nuclease(name="TestNuc", pam=PAM("NGG", "3'"), spacer_length=20)
    validate_nuclease(custom)

def test_validate_nuclease_accepts_string_pam():
    """Nuclease with a plain string PAM should also pass."""
    nuc = Nuclease(name="PlainPAM", pam="NGG", spacer_length=20)
    validate_nuclease(nuc)

def test_validate_nuclease_rejects_non_nuclease():
    with pytest.raises(TypeError, match="Expected a Nuclease"):
        validate_nuclease("SpCas9")

def test_validate_nuclease_rejects_empty_name():
    nuc = Nuclease(name="", pam="NGG", spacer_length=20)
    with pytest.raises(ValueError, match="name must be a non-empty"):
        validate_nuclease(nuc)

def test_validate_nuclease_rejects_whitespace_name():
    nuc = Nuclease(name="   ", pam="NGG", spacer_length=20)
    with pytest.raises(ValueError, match="name must be a non-empty"):
        validate_nuclease(nuc)

def test_validate_nuclease_rejects_invalid_iupac_pam():
    nuc = Nuclease(name="BadPAM", pam=PAM("XZQ", "3'"), spacer_length=20)
    with pytest.raises(ValueError, match="invalid IUPAC"):
        validate_nuclease(nuc)

def test_validate_nuclease_rejects_empty_pam():
    nuc = Nuclease(name="EmptyPAM", pam="", spacer_length=20)
    with pytest.raises(ValueError, match="PAM pattern must not be empty"):
        validate_nuclease(nuc)

def test_validate_nuclease_rejects_short_spacer():
    nuc = Nuclease(name="TooShort", pam="NGG", spacer_length=10)
    with pytest.raises(ValueError, match="spacer_length"):
        validate_nuclease(nuc)

def test_validate_nuclease_rejects_long_spacer():
    nuc = Nuclease(name="TooLong", pam="NGG", spacer_length=35)
    with pytest.raises(ValueError, match="spacer_length"):
        validate_nuclease(nuc)

def test_validate_nuclease_all_default_nucleases_pass():
    """Every nuclease in the default registry should pass validation."""
    for name in list_nucleases():
        nuc = get_nuclease(name)
        validate_nuclease(nuc)  # should not raise

# --- is_nickase and nicking_offset tests ---

def test_is_nickase_ncas9():
    """nCas9 is a registered nickase."""
    assert is_nickase("nCas9") is True

def test_is_nickase_spcas9():
    """Wild-type SpCas9 is NOT a nickase."""
    assert is_nickase("SpCas9") is False

def test_is_nickase_cas12a():
    """AsCas12a is NOT a nickase."""
    assert is_nickase("AsCas12a") is False

def test_is_nickase_with_nuclease_object():
    nuc = get_nuclease("nCas9")
    assert is_nickase(nuc) is True

    nuc2 = get_nuclease("SpCas9")
    assert is_nickase(nuc2) is False

def test_is_nickase_custom_nickase():
    custom = Nuclease(name="MyNickase", pam="NGG", nickase=True)
    assert is_nickase(custom) is True

def test_is_nickase_custom_full_nuclease():
    custom = Nuclease(name="MyNuclease", pam="NGG", nickase=False)
    assert is_nickase(custom) is False

def test_nicking_offset_ncas9():
    """Cas9 nickase nicks at -3 relative to PAM."""
    assert nicking_offset("nCas9") == -3

def test_nicking_offset_with_nuclease_object():
    nuc = get_nuclease("nCas9")
    assert nicking_offset(nuc) == -3

def test_nicking_offset_custom_cas9_nickase():
    custom = Nuclease(name="D10A-Cas9", pam="NGG", nickase=True)
    assert nicking_offset(custom) == -3

def test_nicking_offset_custom_cas12_nickase():
    custom = Nuclease(name="nCas12a", pam="TTTV", nickase=True)
    assert nicking_offset(custom) == 18

def test_nicking_offset_raises_for_full_nuclease():
    """nicking_offset must raise ValueError for non-nickase nucleases."""
    with pytest.raises(ValueError, match="not a nickase"):
        nicking_offset("SpCas9")

def test_nicking_offset_raises_for_cas12a():
    with pytest.raises(ValueError, match="not a nickase"):
        nicking_offset("AsCas12a")
