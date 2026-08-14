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

from .coordinates import CutSite, Interval, Strand, opposite_strand

from .nuclease import (
    Nuclease,
    pam_orientation,
)
from .target import PAMSite, Protospacer, SequenceContext, Spacer, TargetSite

_VALID_DNA_BASES = frozenset("ACGT")
_VALID_TARGET_BASES = frozenset("ACGTN")

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

def normalize(seq: str) -> str:
    if not isinstance(seq, str):
        raise TypeError("Sequence must be a string")
    seq = seq.upper().replace("U", "T")
    if not seq:
        raise ValueError("Sequence must not be empty")
    invalid = set(seq) - _VALID_DNA_BASES
    if invalid:
        raise ValueError(f"Invalid bases in sequence: {sorted(invalid)}")
    return seq

def reverse_complement(seq: str) -> str:
    if not isinstance(seq, str):
        raise TypeError("Sequence must be a string")
    seq = seq.upper().replace("U", "T")
    invalid = set(seq) - _VALID_TARGET_BASES
    if invalid:
        raise ValueError(f"Invalid bases in sequence: {sorted(invalid)}")
    return "".join(_COMPLEMENT[b] for b in reversed(seq))


@dataclass(frozen=True, slots=True)
class Guide:
    """A guide candidate associated with a genomic target site.

    ``sequence`` is the DNA spelling of the guide RNA spacer, stored 5′ to 3′.
    It is identical to the DNA protospacer sequence for the DNA-targeting
    nucleases currently supported. It does not include a guide RNA scaffold.
    All intervals use forward-reference coordinates.
    """

    sequence: str
    nuclease: Nuclease
    pam: str
    pam_interval: Interval
    pam_strand: Strand
    context: SequenceContext | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.nuclease, Nuclease):
            raise TypeError("Guide.nuclease must be a Nuclease")
        if not isinstance(self.pam_interval, Interval):
            raise TypeError("Guide.pam_interval must be an Interval")
        if self.context is not None and not isinstance(
            self.context,
            SequenceContext,
        ):
            raise TypeError(
                "Guide.context must be a SequenceContext or None"
            )
        object.__setattr__(self, "sequence", normalize(self.sequence))
        object.__setattr__(self, "pam", normalize(self.pam))

        if self.pam_strand not in ("+", "-"):
            raise ValueError(
                f"pam_strand must be '+' or '-', got {self.pam_strand!r}"
            )
        if len(self.sequence) != self.nuclease.spacer_length:
            raise ValueError(
                f"Guide length {len(self.sequence)} does not match "
                f"{self.nuclease.name} spacer length "
                f"{self.nuclease.spacer_length}"
            )
        if len(self.pam_interval) != len(self.pam):
            raise ValueError(
                "PAM interval length must equal the observed PAM length"
            )
        if not self.nuclease.matches_pam(self.pam):
            raise ValueError(
                f"PAM {self.pam!r} does not match "
                f"{self.nuclease.name} pattern {self.nuclease.get_pam().pattern!r}"
            )
        if self.context is not None:
            site = TargetSite(
                context=self.context,
                protospacer_interval=self.protospacer_interval,
                pam_interval=self.pam_interval,
                pam_strand=self.pam_strand,
                nuclease=self.nuclease,
            )
            if site.protospacer_sequence != self.sequence:
                raise ValueError(
                    "Guide sequence does not match its SequenceContext"
                )
            if site.pam_sequence != self.pam:
                raise ValueError("Guide PAM does not match its SequenceContext")

    @property
    def target_strand(self) -> Strand:
        """The DNA strand hybridized to the guide RNA."""
        return opposite_strand(self.pam_strand)

    @property
    def spacer(self) -> Spacer:
        """Return the targeting portion of the guide RNA construct."""
        return Spacer(self.sequence)

    @property
    def protospacer_interval(self) -> Interval:
        """Return the protospacer as a forward-reference interval."""
        orientation = pam_orientation(self.nuclease)
        spacer_len = self.nuclease.spacer_length

        # Relative positions are mapped onto forward-reference coordinates.
        spacer_upstream = (
            (self.pam_strand == "+" and orientation == "3'")
            or (self.pam_strand == "-" and orientation == "5'")
        )

        if spacer_upstream:
            return Interval(
                self.pam_interval.start - spacer_len,
                self.pam_interval.start,
            )
        return Interval(
            self.pam_interval.end,
            self.pam_interval.end + spacer_len,
        )

    @property
    def target_site(self) -> TargetSite:
        """Return this guide bound to its validated sequence context."""
        context = self.context
        if context is None:
            orientation = pam_orientation(self.nuclease)
            oriented_site = (
                self.sequence + self.pam
                if orientation == "3'"
                else self.pam + self.sequence
            )
            forward_site = (
                oriented_site
                if self.pam_strand == "+"
                else reverse_complement(oriented_site)
            )
            site_interval = Interval(
                min(self.protospacer_interval.start, self.pam_interval.start),
                max(self.protospacer_interval.end, self.pam_interval.end),
            )
            context = SequenceContext(forward_site, site_interval)

        return TargetSite(
            context=context,
            protospacer_interval=self.protospacer_interval,
            pam_interval=self.pam_interval,
            pam_strand=self.pam_strand,
            nuclease=self.nuclease,
        )

    @property
    def protospacer(self) -> Protospacer:
        """Return the reference-located DNA sequence targeted by the spacer."""
        return self.target_site.protospacer

    @property
    def pam_site(self) -> PAMSite:
        """Return the reference-located PAM match for this guide."""
        return self.target_site.pam_site

    def cut_site(self) -> CutSite:
        """Return the configured strand-aware cut site."""
        return self.target_site.cut_site()

    def cut_position(self) -> int | tuple[int, int]:
        """Return the cut site in the legacy integer-or-tuple form."""
        return self.cut_site().as_position()

def normalize_target(seq: str) -> str:
    if not isinstance(seq, str):
        raise TypeError("Target sequence must be a string")
    seq = seq.upper().replace("U", "T")
    if not seq:
        raise ValueError("Target sequence must not be empty")
    invalid = set(seq) - _VALID_TARGET_BASES
    if invalid:
        raise ValueError(f"Invalid bases in target sequence: {sorted(invalid)}")
    return seq

def find_guides(
    target: str | SequenceContext,
    nuclease: Nuclease,
    pam_strand: Literal["+", "-", "both"] = "both",
) -> list[Guide]:
    """Find guides, reporting all coordinates on the forward reference.

    ``pam_strand`` filters the strand containing the recognized PAM.  Guide
    sequences are returned 5' to 3' in that strand's orientation.
    """
    if pam_strand not in ("+", "-", "both"):
        raise ValueError(
            "pam_strand must be '+', '-', or 'both', "
            f"got {pam_strand!r}"
        )
    if not isinstance(nuclease, Nuclease):
        raise TypeError("nuclease must be a Nuclease instance")

    if isinstance(target, SequenceContext):
        context = target
        target_norm = context.sequence
    else:
        target_norm = normalize_target(target)
        context = SequenceContext(
            target_norm,
            Interval(0, len(target_norm)),
        )
    coordinate_offset = context.interval.start
    pam_obj = nuclease.get_pam()
    orientation = pam_orientation(nuclease)
    pam_len = len(pam_obj)
    spacer_len = nuclease.spacer_length
    seq_len = len(target_norm)

    guides: list[Guide] = []
    if pam_strand in ("+", "both"):
        for p in pam_obj.find_sites(target_norm):
            if orientation == "3'":
                sp_start, sp_end = p - spacer_len, p
            else:
                sp_start = p + pam_len
                sp_end = sp_start + spacer_len

            if sp_start < 0 or sp_end > seq_len:
                continue

            protospacer = target_norm[sp_start:sp_end]
            if set(protospacer) - _VALID_DNA_BASES:
                continue

            guides.append(
                Guide(
                    sequence=protospacer,
                    nuclease=nuclease,
                    pam=target_norm[p : p + pam_len],
                    pam_interval=Interval(
                        coordinate_offset + p,
                        coordinate_offset + p + pam_len,
                    ),
                    pam_strand="+",
                    context=context,
                )
            )

    if pam_strand in ("-", "both"):
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
            if set(protospacer) - _VALID_DNA_BASES:
                continue

            # Map the reverse-oriented PAM back to a forward-reference span.
            pam_start_fwd = seq_len - j - pam_len

            guides.append(
                Guide(
                    sequence=protospacer,
                    nuclease=nuclease,
                    pam=rev_target[j : j + pam_len],
                    pam_interval=Interval(
                        coordinate_offset + pam_start_fwd,
                        coordinate_offset + pam_start_fwd + pam_len,
                    ),
                    pam_strand="-",
                    context=context,
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
    if not isinstance(guide, Guide):
        raise TypeError("guide must be a Guide instance")
    return guide.spacer.gc_content
