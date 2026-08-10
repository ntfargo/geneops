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
        if not isinstance(self.start, int) or not isinstance(self.end, int):
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