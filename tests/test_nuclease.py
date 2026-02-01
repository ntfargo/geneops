import pytest

from geneops.nuclease import (
    Nuclease,
    RRegistry,
    get_nuclease,
    list_nucleases,
    register_nuclease,
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