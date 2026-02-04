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

def matches_pam(seq: str, pam: PAM | str, orientation: Literal["3'", "5'"] = "3'") -> bool:
    if isinstance(pam, str):
        pam = PAM(pam, orientation)
    return pam.matches(seq)

def find_pam_sites(seq: str, pam: PAM | str, orientation: Literal["3'", "5'"] = "3'") -> List[int]:
    if isinstance(pam, str):
        pam = PAM(pam, orientation)
    return pam.find_sites(seq)

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
    spacer_length: int = 20
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
        spacer_length=20,
        description="Streptococcus pyogenes Cas9",
    ),
    Nuclease(
        name="nCas9",
        pam=PAM("NGG", "3'"),
        spacer_length=20,
        nickase=True,
        description="SpCas9 nickase variant (D10A/H840A)",
    ),

    # --- High-fidelity Cas9 ---
    Nuclease(
        name="SpCas9-HF1",
        pam=PAM("NGG", "3'"),
        spacer_length=20,
        description="High-fidelity SpCas9 variant",
    ),
    Nuclease(
        name="eSpCas9",
        pam=PAM("NGG", "3'"),
        spacer_length=20,
        description="Enhanced-specificity SpCas9",
    ),

    # --- PAM-expanded Cas9 ---
    Nuclease(
        name="SpCas9-NG",
        pam=PAM("NG", "3'"),
        spacer_length=20,
        description="SpCas9 variant recognizing NG PAMs",
    ),
    Nuclease(
        name="SpG",
        pam=PAM("NGN", "3'"),
        spacer_length=20,
        description="Engineered SpCas9 with relaxed PAM",
    ),
    Nuclease(
        name="SpRY",
        pam=PAM("NRN", "3'"),
        spacer_length=20,
        description="Near-PAMless SpCas9 variant (context-dependent efficiency)"
    ),

    # --- Smaller Cas9 orthologs ---
    Nuclease(
        name="SaCas9",
        pam=PAM("NNGRRT", "3'"),
        spacer_length=21,
        description="Staphylococcus aureus Cas9",
    ),
    Nuclease(
        name="NmCas9",
        pam=PAM("NNNNGATT", "3'"),
        spacer_length=24,
        description="Neisseria meningitidis Cas9",
    ),
    Nuclease(
        name="CjCas9",
        pam=PAM("NNNNRYAC", "3'"),
        spacer_length=22,
        description="Campylobacter jejuni Cas9",
    ),

    # --- Cas12 family ---
    Nuclease(
        name="AsCas12a",
        pam=PAM("TTTV", "5'"),
        spacer_length=23,
        description="Acidaminococcus Cas12a (Cpf1)",
    ),
    Nuclease(
        name="LbCas12a",
        pam=PAM("TTTV", "5'"),
        spacer_length=23,
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

def cut_offset(nuclease: Nuclease | str) -> int:
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    
    name_lower = nuclease.name.lower()
    
    # Cas12 family cuts downstream of PAM (in the protospacer)
    if 'cas12' in name_lower or 'cpf1' in name_lower:
        return 18  # Typical Cas12a cut position
    
    # Cas9 family cuts 3bp upstream of PAM (between protospacer and PAM)
    return -3


def produces_blunt_cut(nuclease: Nuclease | str) -> bool:
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    
    name_lower = nuclease.name.lower()
    
    # Cas12 family produces staggered cuts
    if 'cas12' in name_lower or 'cpf1' in name_lower:
        return False
    
    # Cas9 family produces blunt cuts
    return True


def overhang_length(nuclease: Nuclease | str) -> int:
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    
    name_lower = nuclease.name.lower()
    
    # Cas12 family produces 5nt 5' overhangs
    if 'cas12' in name_lower or 'cpf1' in name_lower:
        return 5
    
    # Cas9 produces blunt cuts (0 overhang)
    return 0


def cut_position(
    pam_index: int, 
    nuclease: Nuclease | str, 
    strand: Literal["+", "-"] = "+"
) -> int | tuple[int, int]:
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    
    offset = cut_offset(nuclease)
    orientation = pam_orientation(nuclease)
    
    # For 3' PAM (like Cas9): PAM is downstream, cut is upstream
    if orientation == "3'":
        if strand == "+":
            cut_pos = pam_index + offset
        else:
            # Reverse strand: PAM position is at the end, adjust accordingly
            pam_len = nuclease.pam_length()
            cut_pos = pam_index + pam_len - offset
    
    # For 5' PAM (like Cas12a): PAM is upstream, cut is downstream
    else:
        if strand == "+":
            pam_len = nuclease.pam_length()
            cut_pos = pam_index + pam_len + offset
        else:
            cut_pos = pam_index - offset
    
    # Return single position or tuple depending on cut type
    if produces_blunt_cut(nuclease):
        return cut_pos
    else:
        overhang = overhang_length(nuclease)
        return (cut_pos, cut_pos + overhang)


def guide_length(nuclease: Nuclease | str) -> int:
    """Return the required guide length for a nuclease."""
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    
    return nuclease.spacer_length

def is_guide_length_valid(guide: str, nuclease: Nuclease | str) -> bool:
    """Check whether a guide sequence matches the required length for a nuclease."""
    required_length = guide_length(nuclease)
    return len(guide) == required_length

def is_guide_compatible(
    guide: str,
    target: str,
    nuclease: Nuclease | str
) -> bool:
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    
    # Condition 1: Guide length must match nuclease requirement
    if not is_guide_length_valid(guide, nuclease):
        return False
    
    pam = nuclease.get_pam()
    orientation = pam_orientation(nuclease)
    pam_len = len(pam)
    guide_len = len(guide)
    target_upper = target.upper()
    guide_upper = guide.upper()
    
    # For 3' PAM (Cas9): guide binds upstream of PAM on the target strand
    # Target structure: [guide binding site][PAM]
    if orientation == "3'":
        # Check all positions where a PAM could exist
        for i in range(len(target) - pam_len + 1):
            pam_candidate = target_upper[i:i + pam_len]
            if pam.matches(pam_candidate):
                # Guide binding site is upstream of PAM
                guide_start = i - guide_len
                if guide_start >= 0:
                    binding_site = target_upper[guide_start:i]
                    # Guide binds to complementary strand, so check reverse complement
                    if _reverse_complement(guide_upper) == binding_site:
                        return True
    
    # For 5' PAM (Cas12): guide binds downstream of PAM on the target strand
    # Target structure: [PAM][guide binding site]
    else:
        for i in range(len(target) - pam_len + 1):
            pam_candidate = target_upper[i:i + pam_len]
            if pam.matches(pam_candidate):
                # Guide binding site is downstream of PAM
                guide_start = i + pam_len
                guide_end = guide_start + guide_len
                if guide_end <= len(target):
                    binding_site = target_upper[guide_start:guide_end]
                    # Guide binds to complementary strand, so check reverse complement
                    if _reverse_complement(guide_upper) == binding_site:
                        return True
    
    return False

def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))