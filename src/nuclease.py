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

def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                   '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq = list(seq)
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]

def pam_orientation(nuclease: Nuclease | str) -> Literal["3'", "5'"]:
    """Return the explicitly configured PAM orientation for a nuclease."""
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)
    return nuclease.get_pam().orientation

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
class CleavagePattern:
    """Strand-specific cut offsets along the PAM-strand protospacer.

    Both offsets are measured from the protospacer's 5' end in PAM-strand
    orientation.  An offset is a boundary between bases: ``0`` is immediately
    before the protospacer and ``spacer_length`` immediately after it.
    ``None`` means that strand is not cut, allowing nickases to be represented
    without a separate behavior flag.
    """

    pam_strand_offset: int | None
    target_strand_offset: int | None

    def __post_init__(self) -> None:
        offsets = (self.pam_strand_offset, self.target_strand_offset)
        if offsets == (None, None):
            raise ValueError("CleavagePattern must cut at least one strand")
        if any(
            offset is not None
            and (not isinstance(offset, int) or isinstance(offset, bool))
            for offset in offsets
        ):
            raise TypeError("Cleavage offsets must be integers or None")

    @property
    def is_nickase(self) -> bool:
        """Return whether exactly one DNA strand is cut."""
        return (self.pam_strand_offset is None) != (
            self.target_strand_offset is None
        )

    @property
    def is_blunt(self) -> bool:
        """Return whether both strands are cut at the same boundary."""
        return (
            self.pam_strand_offset is not None
            and self.target_strand_offset is not None
            and self.pam_strand_offset == self.target_strand_offset
        )

    @property
    def overhang_length(self) -> int | None:
        """Return staggered-cut length, or ``None`` for a nickase."""
        if self.is_nickase:
            return None
        assert self.pam_strand_offset is not None
        assert self.target_strand_offset is not None
        return abs(self.pam_strand_offset - self.target_strand_offset)

    @property
    def nicking_strand(self) -> Literal["pam", "target"] | None:
        """Return the cut strand for a nickase, otherwise ``None``."""
        if not self.is_nickase:
            return None
        return "pam" if self.pam_strand_offset is not None else "target"


@dataclass(frozen=True, slots=True)
class Nuclease:
    """A nuclease with explicit targeting and strand-cleavage behavior."""

    name: str
    pam: PAM
    cleavage: CleavagePattern
    spacer_length: int = 20
    description: str | None = None

    def get_pam(self) -> PAM:
        """Return this nuclease's explicit PAM definition."""
        if not isinstance(self.pam, PAM):
            raise TypeError("Nuclease.pam must be a PAM instance")
        return self.pam
    
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
        cleavage=CleavagePattern(17, 17),
        spacer_length=20,
        description="Streptococcus pyogenes Cas9",
    ),
    Nuclease(
        name="nCas9",
        pam=PAM("NGG", "3'"),
        cleavage=CleavagePattern(None, 17),
        spacer_length=20,
        description="SpCas9 D10A nickase; cuts the target strand",
    ),

    # --- High-fidelity Cas9 ---
    Nuclease(
        name="SpCas9-HF1",
        pam=PAM("NGG", "3'"),
        cleavage=CleavagePattern(17, 17),
        spacer_length=20,
        description="High-fidelity SpCas9 variant",
    ),
    Nuclease(
        name="eSpCas9",
        pam=PAM("NGG", "3'"),
        cleavage=CleavagePattern(17, 17),
        spacer_length=20,
        description="Enhanced-specificity SpCas9",
    ),

    # --- PAM-expanded Cas9 ---
    Nuclease(
        name="SpCas9-NG",
        pam=PAM("NG", "3'"),
        cleavage=CleavagePattern(17, 17),
        spacer_length=20,
        description="SpCas9 variant recognizing NG PAMs",
    ),
    Nuclease(
        name="SpG",
        pam=PAM("NGN", "3'"),
        cleavage=CleavagePattern(17, 17),
        spacer_length=20,
        description="Engineered SpCas9 with relaxed PAM",
    ),
    Nuclease(
        name="SpRY",
        pam=PAM("NRN", "3'"),
        cleavage=CleavagePattern(17, 17),
        spacer_length=20,
        description="Near-PAMless SpCas9 variant (context-dependent efficiency)"
    ),

    # --- Smaller Cas9 orthologs ---
    Nuclease(
        name="SaCas9",
        pam=PAM("NNGRRT", "3'"),
        cleavage=CleavagePattern(18, 18),
        spacer_length=21,
        description="Staphylococcus aureus Cas9",
    ),
    Nuclease(
        name="NmCas9",
        pam=PAM("NNNNGATT", "3'"),
        cleavage=CleavagePattern(21, 21),
        spacer_length=24,
        description="Neisseria meningitidis Cas9",
    ),
    Nuclease(
        name="CjCas9",
        pam=PAM("NNNNRYAC", "3'"),
        cleavage=CleavagePattern(19, 19),
        spacer_length=22,
        description="Campylobacter jejuni Cas9",
    ),

    # --- Cas12 family ---
    Nuclease(
        name="AsCas12a",
        pam=PAM("TTTV", "5'"),
        cleavage=CleavagePattern(18, 23),
        spacer_length=23,
        description="Acidaminococcus Cas12a (Cpf1)",
    ),
    Nuclease(
        name="LbCas12a",
        pam=PAM("TTTV", "5'"),
        cleavage=CleavagePattern(18, 23),
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

def _resolve_nuclease(nuclease: Nuclease | str) -> Nuclease:
    """Resolve a registry name while leaving explicit objects intact."""
    if isinstance(nuclease, str):
        return DEFAULT_REGISTRY.get_nuclease(nuclease)
    return nuclease


def _offset_from_pam(nuclease: Nuclease, protospacer_offset: int) -> int:
    """Convert a protospacer offset to the legacy PAM-relative form."""
    if pam_orientation(nuclease) == "3'":
        return protospacer_offset - nuclease.spacer_length
    return protospacer_offset


def is_nickase(nuclease: Nuclease | str) -> bool:
    """Return whether ``nuclease`` cuts exactly one DNA strand."""
    nuclease = _resolve_nuclease(nuclease)
    return nuclease.cleavage.is_nickase

def nicking_offset(nuclease: Nuclease | str) -> int:
    """Return a nickase's cut offset relative to the PAM boundary."""
    nuclease = _resolve_nuclease(nuclease)
    if not is_nickase(nuclease):
        raise ValueError(
            f"{nuclease.name} is not a nickase; use cut_offset() for "
            f"double-strand nucleases"
        )

    offset = nuclease.cleavage.pam_strand_offset
    if offset is None:
        offset = nuclease.cleavage.target_strand_offset
    assert offset is not None
    return _offset_from_pam(nuclease, offset)

def cut_offset(nuclease: Nuclease | str) -> int:
    """Return the PAM-strand cut offset relative to the PAM boundary.

    This compatibility helper projects the explicit cleavage pattern into the
    older PAM-relative convention. Use ``Nuclease.cleavage`` when strand
    identity matters.
    """
    nuclease = _resolve_nuclease(nuclease)
    if is_nickase(nuclease):
        raise ValueError(
            f"{nuclease.name} is a nickase; use nicking_offset() instead"
        )

    offset = nuclease.cleavage.pam_strand_offset
    assert offset is not None
    return _offset_from_pam(nuclease, offset)

def produces_blunt_cut(nuclease: Nuclease | str) -> bool:
    """Return whether both strand cuts occur at the same boundary."""
    nuclease = _resolve_nuclease(nuclease)
    return nuclease.cleavage.is_blunt

def overhang_length(nuclease: Nuclease | str) -> int:
    """Return the distance between the two strand cuts."""
    nuclease = _resolve_nuclease(nuclease)
    length = nuclease.cleavage.overhang_length
    if length is None:
        raise ValueError(f"{nuclease.name} is a nickase and has no overhang")
    return length

def cut_position(
    pam_index: int,
    nuclease: Nuclease | str,
    pam_strand: Literal["+", "-"] = "+",
) -> int | tuple[int, int]:
    """Return cut boundaries in zero-based forward-reference coordinates.

    ``pam_index`` is the start of the PAM's half-open interval on the forward
    reference, and ``pam_strand`` is the strand containing the recognized PAM.
    A blunt cut returns one boundary; a staggered cut returns two boundaries.
    """
    nuclease = _resolve_nuclease(nuclease)
    if pam_strand not in {"+", "-"}:
        raise ValueError(
            f"pam_strand must be '+' or '-', got {pam_strand!r}"
        )

    orientation = pam_orientation(nuclease)
    pam_len = nuclease.pam_length()

    if pam_strand == "+":
        protospacer_start = (
            pam_index - nuclease.spacer_length
            if orientation == "3'"
            else pam_index + pam_len
        )

        def position(offset: int) -> int:
            return protospacer_start + offset

    else:
        protospacer_start = (
            pam_index + pam_len + nuclease.spacer_length
            if orientation == "3'"
            else pam_index
        )

        def position(offset: int) -> int:
            return protospacer_start - offset

    positions = [
        position(offset)
        for offset in (
            nuclease.cleavage.pam_strand_offset,
            nuclease.cleavage.target_strand_offset,
        )
        if offset is not None
    ]
    if len(positions) == 1 or positions[0] == positions[1]:
        return positions[0]
    return tuple(sorted(positions))


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
    target_upper = target.upper().replace("U", "T")
    guide_upper = guide.upper().replace("U", "T")
    
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
                    # Guide sequences are reported in PAM-strand orientation.
                    if guide_upper == binding_site:
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
                    # Guide sequences are reported in PAM-strand orientation.
                    if guide_upper == binding_site:
                        return True
    
    return False

def validate_nuclease(nuclease: Nuclease) -> None:
    if not isinstance(nuclease, Nuclease):
        raise TypeError(
            f"Expected a Nuclease instance, got {type(nuclease).__name__}"
        )

    if not isinstance(nuclease.name, str) or not nuclease.name.strip():
        raise ValueError("Nuclease name must be a non-empty string")

    if not isinstance(nuclease.pam, PAM):
        raise TypeError(
            f"PAM must be a PAM instance, got {type(nuclease.pam).__name__}"
        )
    pam_pattern = nuclease.pam.pattern

    if not pam_pattern:
        raise ValueError("PAM pattern must not be empty")

    invalid_chars = {
        ch for ch in pam_pattern.upper() if ch not in IUPAC
    }
    if invalid_chars:
        raise ValueError(
            f"PAM pattern contains invalid IUPAC character(s): "
            f"{sorted(invalid_chars)}"
        )
    if nuclease.pam.orientation not in {"3'", "5'"}:
        raise ValueError(
            "PAM orientation must be \"3'\" or \"5'\", "
            f"got {nuclease.pam.orientation!r}"
        )
    if (
        not isinstance(nuclease.spacer_length, int)
        or isinstance(nuclease.spacer_length, bool)
        or nuclease.spacer_length < 17
        or nuclease.spacer_length > 30
    ):
        raise ValueError(
            f"spacer_length must be an integer in 17–30, "
            f"got {nuclease.spacer_length}"
        )
    if not isinstance(nuclease.cleavage, CleavagePattern):
        raise TypeError(
            "cleavage must be a CleavagePattern instance, "
            f"got {type(nuclease.cleavage).__name__}"
        )

def recognizes_strand(nuclease: Nuclease | str) -> Literal["+", "-", "both"]:
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)

    # All standard Cas9/Cas12 nucleases recognise PAM on both strands
    return "both"

def find_pam_sites_with_strand(seq, pam):
    pam_len = len(pam)

    for i in range(len(seq) - pam_len + 1):
        if pam.matches(seq[i:i + pam_len]):
            yield i, "+"

    rc = _reverse_complement(seq)
    for i in range(len(rc) - pam_len + 1):
        if pam.matches(rc[i:i + pam_len]):
            orig_pos = len(seq) - (i + pam_len)
            yield orig_pos, "-"

def binding_strand(
    guide: str,
    target: str,
    nuclease: Nuclease | str,
) -> Literal["+", "-"] | None:
    if isinstance(nuclease, str):
        nuclease = DEFAULT_REGISTRY.get_nuclease(nuclease)

    pam = nuclease.get_pam()
    orientation = pam_orientation(nuclease)
    guide = guide.upper().replace("U", "T")
    target = target.upper().replace("U", "T")

    for pam_pos, pam_strand in find_pam_sites_with_strand(target, pam):
        bind_strand = "-" if pam_strand == "+" else "+"

        if pam_strand == "+" and orientation == "3'":
            start, end = pam_pos - len(guide), pam_pos
        elif pam_strand == "+" and orientation == "5'":
            start = pam_pos + len(pam)
            end = start + len(guide)
        elif pam_strand == "-" and orientation == "3'":
            start = pam_pos + len(pam)
            end = start + len(guide)
        else:
            start, end = pam_pos - len(guide), pam_pos

        if start < 0 or end > len(target):
            continue

        protospacer = target[start:end]

        if pam_strand == "-":
            protospacer = _reverse_complement(protospacer)

        if protospacer == guide:
            return bind_strand

    return None
