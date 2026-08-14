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
from .nuclease import PAM, Nuclease, calculate_cut_site, pam_orientation

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
class Spacer:
    """The variable targeting sequence used in a guide RNA construct.

    The stored value uses the DNA alphabet so it can be compared directly with
    DNA protospacers. ``rna_sequence`` provides the equivalent RNA spelling.
    A spacer has no genomic location and may target more than one protospacer.
    """

    sequence: str

    def __post_init__(self) -> None:
        if not isinstance(self.sequence, str):
            raise TypeError("Spacer sequence must be a string")
        normalized = self.sequence.upper().replace("U", "T")
        if not normalized:
            raise ValueError("Spacer sequence must not be empty")
        invalid = set(normalized) - _VALID_PROTOSPACER_BASES
        if invalid:
            raise ValueError(f"Invalid bases in spacer: {sorted(invalid)}")
        object.__setattr__(self, "sequence", normalized)

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def dna_sequence(self) -> str:
        """Return the normalized DNA spelling of the spacer."""
        return self.sequence

    @property
    def rna_sequence(self) -> str:
        """Return the spacer 5′ to 3′ using the RNA alphabet."""
        return self.sequence.replace("T", "U")

    @property
    def gc_content(self) -> float:
        """Return the fraction of spacer bases that are G or C."""
        return sum(base in "GC" for base in self.sequence) / len(self)


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
class Protospacer:
    """A reference-located DNA sequence targeted by a spacer.

    ``pam_strand`` identifies the strand whose 5′ to 3′ sequence matches the
    DNA spelling of the guide spacer. Coordinates always remain on the forward
    reference.
    """

    context: SequenceContext
    interval: Interval
    pam_strand: Strand

    def __post_init__(self) -> None:
        if not isinstance(self.context, SequenceContext):
            raise TypeError("Protospacer.context must be a SequenceContext")
        if not isinstance(self.interval, Interval):
            raise TypeError("Protospacer.interval must be an Interval")
        if self.pam_strand not in ("+", "-"):
            raise ValueError(
                f"pam_strand must be '+' or '-', got {self.pam_strand!r}"
            )
        if not self.context.contains(self.interval):
            raise ValueError("Protospacer interval is outside SequenceContext")
        if not self.interval:
            raise ValueError("Protospacer interval must not be empty")
        if set(self.sequence) - _VALID_PROTOSPACER_BASES:
            raise ValueError("Protospacer contains unknown bases")

    def __len__(self) -> int:
        return len(self.interval)

    @property
    def sequence(self) -> str:
        """Return the protospacer 5′ to 3′ on the PAM strand."""
        return self.context.extract(self.interval, self.pam_strand)

    @property
    def target_strand(self) -> Strand:
        """Return the strand complementary to the guide spacer."""
        return opposite_strand(self.pam_strand)


@dataclass(frozen=True, slots=True)
class PAMSite:
    """An observed, reference-located sequence matching a PAM definition."""

    context: SequenceContext
    interval: Interval
    pam_strand: Strand
    pam: PAM

    def __post_init__(self) -> None:
        if not isinstance(self.context, SequenceContext):
            raise TypeError("PAMSite.context must be a SequenceContext")
        if not isinstance(self.interval, Interval):
            raise TypeError("PAMSite.interval must be an Interval")
        if not isinstance(self.pam, PAM):
            raise TypeError("PAMSite.pam must be a PAM")
        if self.pam_strand not in ("+", "-"):
            raise ValueError(
                f"pam_strand must be '+' or '-', got {self.pam_strand!r}"
            )
        if not self.context.contains(self.interval):
            raise ValueError("PAM interval is outside SequenceContext")
        if len(self.interval) != len(self.pam):
            raise ValueError("PAM interval length does not match PAM pattern")
        if not self.pam.matches(self.sequence):
            raise ValueError(
                f"Observed sequence {self.sequence!r} does not match "
                f"PAM pattern {self.pam.pattern!r}"
            )

    @property
    def sequence(self) -> str:
        """Return the observed PAM 5′ to 3′ on the PAM strand."""
        return self.context.extract(self.interval, self.pam_strand)


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
        protospacer = self.protospacer
        pam_site = self.pam_site
        if len(protospacer) != self.nuclease.spacer_length:
            raise ValueError(
                "Protospacer interval length does not match nuclease spacer "
                "length"
            )

        spacer_upstream = (
            (self.pam_strand == "+" and pam_orientation(self.nuclease) == "3'")
            or (
                self.pam_strand == "-"
                and pam_orientation(self.nuclease) == "5'"
            )
        )
        if spacer_upstream:
            adjacent = protospacer.interval.end == pam_site.interval.start
        else:
            adjacent = pam_site.interval.end == protospacer.interval.start
        if not adjacent:
            raise ValueError(
                "Protospacer and PAM intervals do not match nuclease geometry"
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
    def spacer(self) -> Spacer:
        """Return the guide spacer that targets this protospacer."""
        return Spacer(self.protospacer.sequence)

    @property
    def protospacer(self) -> Protospacer:
        """Return the reference-located protospacer object."""
        return Protospacer(
            context=self.context,
            interval=self.protospacer_interval,
            pam_strand=self.pam_strand,
        )

    @property
    def pam_site(self) -> PAMSite:
        """Return the observed, reference-located PAM match."""
        return PAMSite(
            context=self.context,
            interval=self.pam_interval,
            pam_strand=self.pam_strand,
            pam=self.nuclease.get_pam(),
        )

    @property
    def protospacer_sequence(self) -> str:
        """Return the protospacer 5′ to 3′ on the PAM strand."""
        return self.protospacer.sequence

    @property
    def pam_sequence(self) -> str:
        """Return the observed PAM 5′ to 3′ on the PAM strand."""
        return self.pam_site.sequence

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
