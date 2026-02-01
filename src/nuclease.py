from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List

IUPAC: Dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "N": {"A", "C", "G", "T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
}

@dataclass(frozen=True, slots=True)
class Nuclease:
    """Nuclease definition.

    Attributes:
        name: Canonical nuclease name (e.g., SpCas9).
        pam: PAM sequence pattern using IUPAC codes (e.g., NGG).
        nickase: Whether the nuclease is a nickase variant.
        description: Optional free-text description.
    """

    name: str
    pam: str
    nickase: bool = False
    description: str | None = None

    def pam_length(self) -> int:
        return len(self.pam)

    def matches_pam(self, seq: str) -> bool:
        """Return True if seq matches this nuclease's PAM."""
        if len(seq) != len(self.pam):
            return False
        return all(
            base in IUPAC[p]
            for base, p in zip(seq.upper(), self.pam.upper())
        )

class RRegistry:
    """Registry for known nucleases."""

    def __init__(self, nucleases: Iterable[Nuclease] | None = None) -> None:
        self._nucleases: Dict[str, Nuclease] = {}
        if nucleases:
            for nuclease in nucleases:
                self.register_nuclease(nuclease)

    @staticmethod
    def _normalize(name: str) -> str:
        return name.strip().lower()

    def register_nuclease(self, nuclease: Nuclease) -> None:
        key = self._normalize(nuclease.name)
        if key in self._nucleases:
            raise ValueError(f"Nuclease already registered: {nuclease.name}")
        self._nucleases[key] = nuclease

    def get_nuclease(self, name: str) -> Nuclease:
        key = self._normalize(name)
        if key not in self._nucleases:
            raise KeyError(f"Unknown nuclease: {name}")
        return self._nucleases[key]

    def list_nucleases(self) -> List[str]:
        return sorted(nuclease.name for nuclease in self._nucleases.values())


# --- Default nuclease definitions ---
_DEFAULT_NUCLEASES = (
    # --- Canonical Cas9 ---
    Nuclease(
        name="SpCas9",
        pam="NGG",
        description="Streptococcus pyogenes Cas9",
    ),
    Nuclease(
        name="nCas9",
        pam="NGG",
        nickase=True,
        description="SpCas9 nickase variant (D10A/H840A)",
    ),

    # --- High-fidelity Cas9 ---
    Nuclease(
        name="SpCas9-HF1",
        pam="NGG",
        description="High-fidelity SpCas9 variant",
    ),
    Nuclease(
        name="eSpCas9",
        pam="NGG",
        description="Enhanced-specificity SpCas9",
    ),

    # --- PAM-expanded Cas9 ---
    Nuclease(
        name="SpCas9-NG",
        pam="NG",
        description="SpCas9 variant recognizing NG PAMs",
    ),
    Nuclease(
        name="SpG",
        pam="NGN",
        description="Engineered SpCas9 with relaxed PAM",
    ),
    Nuclease(
        name="SpRY",
        pam="NRN",
        description="Near-PAMless SpCas9 variant",
    ),

    # --- Smaller Cas9 orthologs ---
    Nuclease(
        name="SaCas9",
        pam="NNGRRT",
        description="Staphylococcus aureus Cas9",
    ),
    Nuclease(
        name="NmCas9",
        pam="NNNNGATT",
        description="Neisseria meningitidis Cas9",
    ),
    Nuclease(
        name="CjCas9",
        pam="NNNVRYAC",
        description="Campylobacter jejuni Cas9",
    ),

    # --- Cas12 family ---
    Nuclease(
        name="AsCas12a",
        pam="TTTV",
        description="Acidaminococcus Cas12a (Cpf1)",
    ),
    Nuclease(
        name="LbCas12a",
        pam="TTTV",
        description="Lachnospiraceae Cas12a (Cpf1)",
    ),
)

DEFAULT_REGISTRY = RRegistry(_DEFAULT_NUCLEASES)
 
def register_nuclease(nuclease: Nuclease) -> None:
    """Register a user-defined nuclease in the default registry."""
    DEFAULT_REGISTRY.register_nuclease(nuclease)

def get_nuclease(name: str) -> Nuclease:
    """Retrieve a known nuclease definition from the default registry."""
    return DEFAULT_REGISTRY.get_nuclease(name)

def list_nucleases() -> List[str]:
    """Enumerate supported nucleases from the default registry."""
    return DEFAULT_REGISTRY.list_nucleases()
