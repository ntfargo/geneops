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
from typing import Literal

from .nuclease import (
    Nuclease,
    cut_position,
    is_guide_compatible,
    pam_orientation,
)
 
_VALID_DNA_BASES = frozenset("ACGT")

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

def normalize(seq: str) -> str:
    seq = seq.upper().replace("U", "T")
    if not seq:
        raise ValueError("Sequence must not be empty")
    invalid = set(seq) - _VALID_DNA_BASES
    if invalid:
        raise ValueError(f"Invalid bases in sequence: {sorted(invalid)}")
    return seq

def reverse_complement(seq: str) -> str:
    return "".join(_COMPLEMENT[b] for b in reversed(seq))

@dataclass(frozen=True)
class Guide:
    sequence: str
    nuclease: Nuclease
    pam: str
    pam_position: int
    strand: Literal["+", "-"]

    def __post_init__(self) -> None:
        object.__setattr__(self, "sequence", normalize(self.sequence))

    def binding_interval(self) -> tuple[int, int]:
        orientation = pam_orientation(self.nuclease)
        pam_len = len(self.pam)
        spacer_len = self.nuclease.spacer_length

        # Spacer sits upstream of PAM on the strand with the PAM:
        #   (+, 3') or (-, 5')  →  spacer is left of PAM in fwd coords
        # Spacer sits downstream of PAM:
        #   (+, 5') or (-, 3')  →  spacer is right of PAM in fwd coords
        spacer_upstream = (
            (self.strand == "+" and orientation == "3'")
            or (self.strand == "-" and orientation == "5'")
        )

        if spacer_upstream:
            return (self.pam_position - spacer_len, self.pam_position)
        return (
            self.pam_position + pam_len,
            self.pam_position + pam_len + spacer_len,
        )

    def cut_site(self) -> int | tuple[int, int]:
        return cut_position(self.pam_position, self.nuclease, self.strand)

    def binds_strand(self) -> Literal["+", "-"]:
        return self.strand

def normalize_target(seq: str) -> str:
    return seq.upper().replace("U", "T")

def find_guides(
    target: str,
    nuclease: Nuclease,
    strand: Literal["+", "-", "both"] = "both",
) -> list[Guide]:
    target_norm = normalize_target(target)
    pam_obj = nuclease.get_pam()
    orientation = pam_orientation(nuclease)
    pam_len = len(pam_obj)
    spacer_len = nuclease.spacer_length
    seq_len = len(target_norm)

    guides: list[Guide] = []
    if strand in ("+", "both"):
        for p in pam_obj.find_sites(target_norm):
            if orientation == "3'":
                sp_start, sp_end = p - spacer_len, p
            else:
                sp_start = p + pam_len
                sp_end = sp_start + spacer_len

            if sp_start < 0 or sp_end > seq_len:
                continue

            protospacer = target_norm[sp_start:sp_end]
            guide_seq = reverse_complement(protospacer)

            guides.append(
                Guide(
                    sequence=guide_seq,
                    nuclease=nuclease,
                    pam=target_norm[p : p + pam_len],
                    pam_position=p,
                    strand="+",
                )
            )

    if strand in ("-", "both"):
        rev_target = reverse_complement(target_norm)
        for j in pam_obj.find_sites(rev_target):
            if orientation == "3'":
                sp_start, sp_end = j - spacer_len, j
            else:
                sp_start = j + pam_len
                sp_end = sp_start + spacer_len

            if sp_start < 0 or sp_end > seq_len:
                continue

            protospacer = rev_target[sp_start:sp_end]
            guide_seq = reverse_complement(protospacer)

            # Map PAM position back to forward-strand coordinates
            pam_pos_fwd = seq_len - j - pam_len

            guides.append(
                Guide(
                    sequence=guide_seq,
                    nuclease=nuclease,
                    pam=target_norm[pam_pos_fwd : pam_pos_fwd + pam_len],
                    pam_position=pam_pos_fwd,
                    strand="-",
                )
            )

    return guides

def validate_guide_sequence(seq: str) -> None:
    normalized = normalize(seq)  # handles empty / bad alphabet
    if len(normalized) < 17 or len(normalized) > 30:
        raise ValueError(
            f"Guide length {len(normalized)} is outside the expected "
            f"range (17–30 nt)"
        )

def guide_gc_content(guide: Guide) -> float:
    seq = guide.sequence
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq)