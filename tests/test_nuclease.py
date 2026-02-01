import pytest

from geneops.nuclease import Nuclease, RRegistry, get_nuclease, list_nucleases, register_nuclease

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

# Example.
print("Supported nucleases:")
for name in list_nucleases():
    print(f"  - {name}")

spcas9 = get_nuclease("spcas9")

print("\nRetrieved nuclease:")
print(f"Name       : {spcas9.name}")
print(f"PAM        : {spcas9.pam}")
print(f"Nickase    : {spcas9.nickase}")
print(f"Description: {spcas9.description}")

custom = Nuclease(
    name="MyCas9",
    pam="NGA",
    description="Custom engineered Cas9 variant",
)

register_nuclease(custom)

mycas9 = get_nuclease("mycas9")
print("\nCustom nuclease registered:")
print(f"{mycas9.name} recognizes PAM {mycas9.pam}")