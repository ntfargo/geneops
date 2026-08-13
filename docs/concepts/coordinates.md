# Coordinates and strands

Genome-editing software becomes difficult to compose when each component makes different assumptions about indices and strand orientation. Geneops therefore uses one coordinate contract throughout its public core.

## Zero-based, half-open intervals

An interval `[start, end)` includes `start` and excludes `end`. Its length is always `end - start`.

<div class="coordinate-strip">
  <span>10 · A</span><span>11 · C</span><span>12 · G</span><span>13 · T</span>
</div>

The four bases above occupy `Interval(10, 14)`.

```python
from geneops import Interval

region = Interval(10, 14)

len(region)          # 4
region.contains(10)  # True
region.contains(13)  # True
region.contains(14)  # False
```

This convention matches Python slicing and BED intervals, making conversions predictable.

## Always use the forward reference

All Geneops intervals refer to the supplied sequence in its forward orientation. This remains true when a PAM is discovered on the reverse strand.

For a [`Guide`][geneops.guide.Guide]:

- `pam_interval` is the PAM span on the forward reference.
- `protospacer_interval` is the protospacer span on the forward reference.
- `pam_strand` is the strand that contains the recognized PAM.
- `target_strand` is the complementary strand hybridized by the guide RNA.
- `sequence` is stored 5′ to 3′ in the PAM-strand orientation.

!!! important
    `pam_strand` and `target_strand` are opposites. Avoid using an unqualified field named only `strand` in downstream data structures.

## Cut sites are strand-aware boundaries

Cut positions are zero-based boundaries between bases, not base indices. A [`CutSite`][geneops.coordinates.CutSite] keeps the PAM-strand and target-strand boundaries separate, which matters for staggered cuts and nickases.

```python
from geneops import calculate_cut_site

site = calculate_cut_site(10, "AsCas12a", "+")

site.pam_strand_boundary     # 32
site.target_strand_boundary  # 37
site.plus_strand_boundary    # 32
site.minus_strand_boundary   # 37
site.overhang                # "5'"
```

The first argument is the start of the PAM's half-open interval on the forward reference. Use `guide.cut_site()` when you already have a guide. Existing coordinate-only code can continue to use `cut_position()`, which returns an integer for one unique boundary or a sorted tuple for two.

!!! note "Changed in 0.0.6"
    `Guide.cut_site()` now returns a strand-aware `CutSite` instead of an integer or tuple. Replace numeric uses with `guide.cut_position()`, or use properties such as `guide.cut_site().pam_strand_boundary` when strand identity matters.

Cut coordinates are the nominal boundaries configured by the selected nuclease. Real cleavage can be heterogeneous across molecules, target sequences, and experimental conditions.

## Converting one-based coordinates

For a one-based, fully closed interval `[start, end]`, subtract one from the start and leave the end unchanged:

```python
from geneops import Interval

one_based_start = 11
one_based_end = 14
interval = Interval(one_based_start - 1, one_based_end)
assert interval == Interval(10, 14)
```
