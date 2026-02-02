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

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Literal

# IUPAC nucleotide codes mapping
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

def pam_orientation(nuclease: Nuclease | str) -> Literal["3'", "5'"]:
    """Determine PAM orientation for a nuclease."""
    # Forward reference - will be resolved at runtime
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    
    name_lower = nuclease.name.lower()
    
    # Cas12 family (Cpf1) has 5' PAM
    if 'cas12' in name_lower or 'cpf1' in name_lower:
        return "5'"
    
    # Cas9 family and most others have 3' PAM
    return "3'"


def matches_pam(seq: str, pam: str) -> bool:
    """Check if a sequence matches a PAM pattern (IUPAC-aware)."""
    return PAM(pam, "3'").matches(seq)


def find_pam_sites(seq: str, pam: str) -> List[int]:
    """Locate all PAM site indices in a sequence."""
    return PAM(pam, "3'").find_sites(seq)

@dataclass(frozen=True, slots=True)
class PAM:
    """PAM sequence pattern with IUPAC-aware matching."""
    pattern: str
    orientation: Literal["5'", "3'"]

    def __len__(self) -> int:
        """Return the length of the PAM pattern."""
        return len(self.pattern)
    
    def __eq__(self, other) -> bool:
        """Support comparison with strings (pattern only) or other PAM objects."""
        if isinstance(other, str):
            return self.pattern == other
        return super().__eq__(other)
    
    def matches(self, seq: str) -> bool:
        """Check if a sequence matches this PAM pattern."""
        if len(seq) != len(self.pattern):
            return False
        if any(p not in IUPAC for p in self.pattern.upper()):
            raise ValueError(f"Invalid IUPAC code in PAM: {self.pattern}")
        return all(
            base in IUPAC[p]
            for base, p in zip(seq.upper(), self.pattern.upper())
        )
    
    def find_sites(self, seq: str) -> List[int]:
        """Locate all PAM site indices in a sequence."""
        pam_len = len(self.pattern)
        seq_upper = seq.upper()
        sites = []
        
        for i in range(len(seq) - pam_len + 1):
            candidate = seq_upper[i:i + pam_len]
            if self.matches(candidate):
                sites.append(i)
        return sites

@dataclass(frozen=True, slots=True)
class Nuclease:
    name: str
    pam: PAM | str
    nickase: bool = False
    description: str | None = None

    def get_pam(self) -> PAM:
        """Return a PAM object for this nuclease."""
        if isinstance(self.pam, PAM):
            return self.pam
        # If pam is a string, determine orientation based on nuclease type
        return PAM(self.pam, pam_orientation(self))
    
    def pam_length(self) -> int:
        """Return the length of the PAM sequence."""
        return len(self.pam)

    def matches_pam(self, seq: str) -> bool:
        """Return True if seq matches this nuclease's PAM."""
        return self.get_pam().matches(seq)

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

_DEFAULT_NUCLEASES = (
    # --- Canonical Cas9 ---
    Nuclease(
        name="SpCas9",
        pam=PAM("NGG", "3'"),
        description="Streptococcus pyogenes Cas9",
    ),
    Nuclease(
        name="nCas9",
        pam=PAM("NGG", "3'"),
        nickase=True,
        description="SpCas9 nickase variant (D10A/H840A)",
    ),

    # --- High-fidelity Cas9 ---
    Nuclease(
        name="SpCas9-HF1",
        pam=PAM("NGG", "3'"),
        description="High-fidelity SpCas9 variant",
    ),
    Nuclease(
        name="eSpCas9",
        pam=PAM("NGG", "3'"),
        description="Enhanced-specificity SpCas9",
    ),

    # --- PAM-expanded Cas9 ---
    Nuclease(
        name="SpCas9-NG",
        pam=PAM("NG", "3'"),
        description="SpCas9 variant recognizing NG PAMs",
    ),
    Nuclease(
        name="SpG",
        pam=PAM("NGN", "3'"),
        description="Engineered SpCas9 with relaxed PAM",
    ),
    Nuclease(
        name="SpRY",
        pam=PAM("NRN", "3'"),
        description="Near-PAMless SpCas9 variant (context-dependent efficiency)"
    ),

    # --- Smaller Cas9 orthologs ---
    Nuclease(
        name="SaCas9",
        pam=PAM("NNGRRT", "3'"),
        description="Staphylococcus aureus Cas9",
    ),
    Nuclease(
        name="NmCas9",
        pam=PAM("NNNNGATT", "3'"),
        description="Neisseria meningitidis Cas9",
    ),
    Nuclease(
        name="CjCas9",
        pam=PAM("NNNNRYAC", "3'"),
        description="Campylobacter jejuni Cas9",
    ),

    # --- Cas12 family ---
    Nuclease(
        name="AsCas12a",
        pam=PAM("TTTV", "5'"),
        description="Acidaminococcus Cas12a (Cpf1)",
    ),
    Nuclease(
        name="LbCas12a",
        pam=PAM("TTTV", "5'"),
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
