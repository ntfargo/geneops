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

"""geneops src base"""

from ._version import __version__

from .coordinates import CutSite, Interval, Overhang, Strand, opposite_strand

from .nuclease import (
    CleavagePattern,
    Nuclease,
    PAM,
    RRegistry,
    get_nuclease,
    list_nucleases,
    register_nuclease,
    is_nickase,
    nicking_offset,
    calculate_cut_site,
    cut_position,
    cut_offset,
    produces_blunt_cut,
    overhang_length,
    validate_nuclease,
)

from .target import SequenceContext, TargetSite

from .guide import (
    Guide,
    normalize,
    reverse_complement,
    find_guides,
    validate_guide_sequence,
    guide_gc_content,
)

__all__ = [
    # coordinates
    "Interval",
    "CutSite",
    "Overhang",
    "Strand",
    "opposite_strand",
    # target
    "SequenceContext",
    "TargetSite",
    # nuclease
    "CleavagePattern",
    "Nuclease",
    "PAM",
    "RRegistry",
    "get_nuclease",
    "list_nucleases",
    "register_nuclease",
    "is_nickase",
    "nicking_offset",
    "calculate_cut_site",
    "cut_position",
    "cut_offset",
    "produces_blunt_cut",
    "overhang_length",
    "validate_nuclease",
    # guide
    "Guide",
    "normalize",
    "reverse_complement",
    "find_guides",
    "validate_guide_sequence",
    "guide_gc_content",
]
