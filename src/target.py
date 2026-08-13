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

"""Reference-aware sequence and genome-editing target objects."""

from __future__ import annotations

from dataclasses import dataclass

from .coordinates import CutSite, Interval, Strand, opposite_strand
from .nuclease import Nuclease, calculate_cut_site, pam_orientation

_VALID_CONTEXT_BASES = frozenset("ACGTN")
_VALID_PROTOSPACER_BASES = frozenset("ACGT")
_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _normalize_context_sequence(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise TypeError("SequenceContext sequence must be a string")
    normalized = sequence.upper().replace("U", "T")
    if not normalized:
        raise ValueError("SequenceContext sequence must not be empty")
    invalid = set(normalized) - _VALID_CONTEXT_BASES
    if invalid:
        raise ValueError(
            f"Invalid bases in SequenceContext: {sorted(invalid)}"
        )
    return normalized


def _reverse_complement(sequence: str) -> str:
    return "".join(_COMPLEMENT[base] for base in reversed(sequence))


@dataclass(frozen=True, slots=True)
class SequenceContext:
    """A forward-reference sequence slice with stable coordinates.

    ``interval`` identifies where ``sequence`` lies on its reference and must
    have the same length. ``reference`` and ``contig`` are optional identifiers
    so small, anonymous sequences remain useful while genome-backed workflows
    can preserve their assembly or accession and contig identity.
    """

    sequence: str
    interval: Interval
    reference: str | None = None
    contig: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.interval, Interval):
            raise TypeError("SequenceContext.interval must be an Interval")
        normalized = _normalize_context_sequence(self.sequence)
        object.__setattr__(self, "sequence", normalized)
        if len(normalized) != len(self.interval):
            raise ValueError(
                "SequenceContext interval length must equal sequence length"
            )
        for name, value in (
            ("reference", self.reference),
            ("contig", self.contig),
        ):
            if value is not None and (
                not isinstance(value, str) or not value.strip()
            ):
                raise ValueError(f"{name} must be a non-empty string or None")

    def contains(self, interval: Interval) -> bool:
        """Return whether *interval* lies fully inside this context."""
        return self.interval.contains_interval(interval)

    def extract(
        self,
        interval: Interval,
        strand: Strand = "+",
    ) -> str:
        """Extract an interval 5′ to 3′ on the requested strand."""
        if strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {strand!r}")
        if not self.contains(interval):
            raise ValueError(
                f"Interval {interval!r} is outside context {self.interval!r}"
            )
        local_start = interval.start - self.interval.start
        local_end = interval.end - self.interval.start
        sequence = self.sequence[local_start:local_end]
        if strand == "-":
            return _reverse_complement(sequence)
        return sequence

    def extract_flanked(
        self,
        interval: Interval,
        *,
        upstream: int = 0,
        downstream: int = 0,
        strand: Strand = "+",
    ) -> str:
        """Extract a feature plus strand-relative upstream/downstream bases."""
        for name, value in (("upstream", upstream), ("downstream", downstream)):
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 0
            ):
                raise ValueError(f"{name} must be a non-negative integer")
        if strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {strand!r}")

        if strand == "+":
            start = interval.start - upstream
            end = interval.end + downstream
        else:
            start = interval.start - downstream
            end = interval.end + upstream
        if start < self.interval.start or end > self.interval.end:
            raise ValueError("Requested flanks extend outside SequenceContext")
        return self.extract(Interval(start, end), strand)


@dataclass(frozen=True, slots=True)
class TargetSite:
    """A validated protospacer/PAM site bound to sequence context.

    Intervals always use forward-reference coordinates. Sequence properties
    are returned 5′ to 3′ in the PAM-strand orientation, matching
    :class:`geneops.guide.Guide`.
    """

    context: SequenceContext
    protospacer_interval: Interval
    pam_interval: Interval
    pam_strand: Strand
    nuclease: Nuclease

    def __post_init__(self) -> None:
        if not isinstance(self.context, SequenceContext):
            raise TypeError("TargetSite.context must be a SequenceContext")
        if not isinstance(self.protospacer_interval, Interval):
            raise TypeError(
                "TargetSite.protospacer_interval must be an Interval"
            )
        if not isinstance(self.pam_interval, Interval):
            raise TypeError("TargetSite.pam_interval must be an Interval")
        if not isinstance(self.nuclease, Nuclease):
            raise TypeError("TargetSite.nuclease must be a Nuclease")
        if self.pam_strand not in ("+", "-"):
            raise ValueError(
                f"pam_strand must be '+' or '-', got {self.pam_strand!r}"
            )
        if len(self.protospacer_interval) != self.nuclease.spacer_length:
            raise ValueError(
                "Protospacer interval length does not match nuclease spacer "
                "length"
            )
        if len(self.pam_interval) != self.nuclease.pam_length():
            raise ValueError("PAM interval length does not match nuclease PAM")
        if not self.context.contains(self.protospacer_interval):
            raise ValueError("Protospacer interval is outside SequenceContext")
        if not self.context.contains(self.pam_interval):
            raise ValueError("PAM interval is outside SequenceContext")

        spacer_upstream = (
            (self.pam_strand == "+" and pam_orientation(self.nuclease) == "3'")
            or (
                self.pam_strand == "-"
                and pam_orientation(self.nuclease) == "5'"
            )
        )
        if spacer_upstream:
            adjacent = self.protospacer_interval.end == self.pam_interval.start
        else:
            adjacent = self.pam_interval.end == self.protospacer_interval.start
        if not adjacent:
            raise ValueError(
                "Protospacer and PAM intervals do not match nuclease geometry"
            )

        if set(self.protospacer_sequence) - _VALID_PROTOSPACER_BASES:
            raise ValueError("TargetSite protospacer contains unknown bases")
        if not self.nuclease.matches_pam(self.pam_sequence):
            raise ValueError(
                f"Observed PAM {self.pam_sequence!r} does not match "
                f"{self.nuclease.name}"
            )

    @property
    def target_strand(self) -> Strand:
        """Return the DNA strand hybridized by the guide RNA."""
        return opposite_strand(self.pam_strand)

    @property
    def site_interval(self) -> Interval:
        """Return the smallest interval containing protospacer and PAM."""
        return Interval(
            min(self.protospacer_interval.start, self.pam_interval.start),
            max(self.protospacer_interval.end, self.pam_interval.end),
        )

    @property
    def protospacer_sequence(self) -> str:
        """Return the protospacer 5′ to 3′ on the PAM strand."""
        return self.context.extract(
            self.protospacer_interval,
            self.pam_strand,
        )

    @property
    def pam_sequence(self) -> str:
        """Return the observed PAM 5′ to 3′ on the PAM strand."""
        return self.context.extract(self.pam_interval, self.pam_strand)

    def sequence_context(
        self,
        *,
        upstream: int = 0,
        downstream: int = 0,
    ) -> str:
        """Return the oriented site plus requested sequence flanks."""
        return self.context.extract_flanked(
            self.site_interval,
            upstream=upstream,
            downstream=downstream,
            strand=self.pam_strand,
        )

    def cut_site(self) -> CutSite:
        """Return configured strand-aware cut coordinates for this target."""
        return calculate_cut_site(
            self.pam_interval.start,
            self.nuclease,
            self.pam_strand,
        )
