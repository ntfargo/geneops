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

"""
Coordinate primitives shared by genome-editing workflows.
Geneops uses zero-based, half-open intervals throughout.  An interval
``[start, end)`` contains ``start`` but not ``end``; its length is therefore
``end - start``.  Coordinates always refer to the forward reference sequence,
even when a feature is found on the reverse strand.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Literal

Strand = Literal["+", "-"]
Overhang = Literal["5'", "3'"]

def opposite_strand(strand: Strand) -> Strand:
    """Return the DNA strand complementary to *strand*."""
    if strand == "+":
        return "-"
    if strand == "-":
        return "+"
    raise ValueError(f"Strand must be '+' or '-', got {strand!r}")

@dataclass(frozen=True, slots=True, order=True)
class Interval:
    """A zero-based, half-open interval on the forward reference sequence."""

    start: int
    end: int

    def __post_init__(self) -> None:
        if (
            not isinstance(self.start, int)
            or isinstance(self.start, bool)
            or not isinstance(self.end, int)
            or isinstance(self.end, bool)
        ):
            raise TypeError("Interval coordinates must be integers")
        if self.start < 0:
            raise ValueError("Interval start must be non-negative")
        if self.end < self.start:
            raise ValueError("Interval end must be greater than or equal to start")

    def __len__(self) -> int:
        return self.end - self.start

    def __iter__(self) -> Iterator[int]:
        yield self.start
        yield self.end

    def contains(self, position: int) -> bool:
        """Return whether a zero-based base position lies in this interval."""
        return self.start <= position < self.end

    def contains_interval(self, interval: Interval) -> bool:
        """Return whether *interval* is fully contained in this interval."""
        if not isinstance(interval, Interval):
            raise TypeError("interval must be an Interval instance")
        return self.start <= interval.start and interval.end <= self.end


@dataclass(frozen=True, slots=True)
class CutSite:
    """Configured strand-aware cut boundaries on a forward reference.

    Boundaries are zero-based positions between bases. ``pam_strand_boundary``
    and ``target_strand_boundary`` retain the biological identity of each cut;
    ``None`` means that strand is not cut. ``pam_strand`` maps those roles onto
    the forward (``+``) and reverse (``-``) reference strands. These positions
    represent the nuclease's configured nominal geometry, not a claim that
    every molecule is cleaved at exactly the same boundary.
    """

    pam_strand_boundary: int | None
    target_strand_boundary: int | None
    pam_strand: Strand

    def __post_init__(self) -> None:
        boundaries = (
            self.pam_strand_boundary,
            self.target_strand_boundary,
        )
        if boundaries == (None, None):
            raise ValueError("CutSite must contain at least one cut boundary")
        if any(
            boundary is not None
            and (
                not isinstance(boundary, int)
                or isinstance(boundary, bool)
            )
            for boundary in boundaries
        ):
            raise TypeError("Cut boundaries must be integers or None")
        if any(
            boundary is not None and boundary < 0
            for boundary in boundaries
        ):
            raise ValueError("Cut boundaries must be non-negative")
        if self.pam_strand not in ("+", "-"):
            raise ValueError(
                f"pam_strand must be '+' or '-', got {self.pam_strand!r}"
            )

    @property
    def target_strand(self) -> Strand:
        """Return the strand complementary to the PAM strand."""
        return opposite_strand(self.pam_strand)

    @property
    def plus_strand_boundary(self) -> int | None:
        """Return the cut boundary on the forward reference strand."""
        if self.pam_strand == "+":
            return self.pam_strand_boundary
        return self.target_strand_boundary

    @property
    def minus_strand_boundary(self) -> int | None:
        """Return the cut boundary on the reverse reference strand."""
        if self.pam_strand == "-":
            return self.pam_strand_boundary
        return self.target_strand_boundary

    @property
    def is_nickase(self) -> bool:
        """Return whether exactly one DNA strand is cut."""
        return (self.pam_strand_boundary is None) != (
            self.target_strand_boundary is None
        )

    @property
    def is_blunt(self) -> bool:
        """Return whether both strands are cut at the same boundary."""
        return (
            self.pam_strand_boundary is not None
            and self.target_strand_boundary is not None
            and self.pam_strand_boundary == self.target_strand_boundary
        )

    @property
    def nicking_strand(self) -> Strand | None:
        """Return the reference strand cut by a nickase, otherwise ``None``."""
        if not self.is_nickase:
            return None
        if self.plus_strand_boundary is not None:
            return "+"
        return "-"

    @property
    def overhang_length(self) -> int | None:
        """Return staggered-cut length, or ``None`` for a nickase."""
        if self.is_nickase:
            return None
        assert self.pam_strand_boundary is not None
        assert self.target_strand_boundary is not None
        return abs(
            self.pam_strand_boundary - self.target_strand_boundary
        )

    @property
    def overhang(self) -> Overhang | None:
        """Return overhang polarity, or ``None`` for blunt cuts and nicks."""
        if self.is_nickase or self.is_blunt:
            return None
        assert self.plus_strand_boundary is not None
        assert self.minus_strand_boundary is not None
        if self.plus_strand_boundary < self.minus_strand_boundary:
            return "5'"
        return "3'"

    @property
    def boundaries(self) -> tuple[int, ...]:
        """Return unique cut boundaries in ascending reference order."""
        return tuple(
            sorted(
                {
                    boundary
                    for boundary in (
                        self.pam_strand_boundary,
                        self.target_strand_boundary,
                    )
                    if boundary is not None
                }
            )
        )

    def as_position(self) -> int | tuple[int, int]:
        """Project this site into the legacy integer-or-tuple representation."""
        if len(self.boundaries) == 1:
            return self.boundaries[0]
        first, second = self.boundaries
        return first, second
