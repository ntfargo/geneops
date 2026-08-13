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

import pytest

from geneops.coordinates import Interval, opposite_strand
class TestInterval:
    def test_zero_based_half_open_length(self):
        interval = Interval(10, 20)
        assert len(interval) == 10
        assert tuple(interval) == (10, 20)

    def test_contains_excludes_end(self):
        interval = Interval(10, 20)
        assert interval.contains(10)
        assert interval.contains(19)
        assert not interval.contains(20)

    def test_empty_interval_is_valid(self):
        assert len(Interval(5, 5)) == 0

    @pytest.mark.parametrize("start,end", [(-1, 2), (5, 4)])
    def test_invalid_bounds_raise(self, start, end):
        with pytest.raises(ValueError):
            Interval(start, end)

    def test_non_integer_bounds_raise(self):
        with pytest.raises(TypeError):
            Interval(1.5, 2)  # type: ignore[arg-type]

    @pytest.mark.parametrize("start,end", [(True, 2), (1, False)])
    def test_boolean_bounds_raise(self, start, end):
        with pytest.raises(TypeError, match="integers"):
            Interval(start, end)
class TestOppositeStrand:
    def test_opposites(self):
        assert opposite_strand("+") == "-"
        assert opposite_strand("-") == "+"

    def test_invalid_strand_raises(self):
        with pytest.raises(ValueError, match="Strand"):
            opposite_strand("both")  # type: ignore[arg-type]
