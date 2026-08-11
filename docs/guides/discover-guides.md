# Discover guides

[`find_guides`][geneops.guide.find_guides] searches a target sequence for recognized PAMs and returns valid, construct-ready [`Guide`][geneops.guide.Guide] objects.

## Cas9: a 3′ PAM

For SpCas9, the protospacer is upstream of the PAM when read in the PAM strand's orientation:

```text
5′  [-------- 20 nt protospacer --------][NGG]  3′
```

```python
from geneops import find_guides, get_nuclease

spacer = "ATGCACGTTAGCTACGATCA"
target = "A" * 10 + spacer + "TGG" + "A" * 10

guide = find_guides(
    target,
    get_nuclease("SpCas9"),
    pam_strand="+",
)[0]

assert guide.sequence == spacer
assert guide.pam == "TGG"
assert tuple(guide.pam_interval) == (30, 33)
assert tuple(guide.protospacer_interval) == (10, 30)
assert guide.target_strand == "-"
assert guide.cut_site() == 27
```

## Cas12a: a 5′ PAM

Cas12a reverses the local geometry: its PAM is upstream of the protospacer.

```text
5′  [TTTV][---------- 23 nt protospacer ----------]  3′
```

```python
spacer = "ACGTACGTACGTACGTACGTACG"
target = "C" * 10 + "TTTG" + spacer + "C" * 10

guide = find_guides(
    target,
    get_nuclease("AsCas12a"),
    pam_strand="+",
)[0]

assert guide.sequence == spacer
assert tuple(guide.pam_interval) == (10, 14)
assert tuple(guide.protospacer_interval) == (14, 37)
assert guide.cut_site() == (32, 37)
```

## Handle ambiguous target bases

Target sequences may contain `N`, which is useful for masked or uncertain bases. A candidate is not returned when its protospacer contains an `N`, because `Guide.sequence` is construct-ready and accepts only `A`, `C`, `G`, and `T`.

```python
from geneops.guide import normalize_target

normalize_target("acnug")  # ACNTG
```

## Rank or filter results

Geneops currently provides GC content as a small composable primitive:

```python
from geneops import guide_gc_content

filtered = [
    guide
    for guide in find_guides(target, get_nuclease("AsCas12a"))
    if 0.35 <= guide_gc_content(guide) <= 0.70
]
```

More advanced activity and off-target scores belong to later roadmap phases. Keeping discovery and scoring separate lets a workflow choose the model appropriate for its assay.
