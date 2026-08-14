# Discover guides

[`find_guides`][geneops.guide.find_guides] searches a target sequence for recognized PAMs and returns validated [`Guide`][geneops.guide.Guide] candidates. A candidate contains the targeting spacer and its reference-located protospacer and PAM it does not include a complete guide RNA scaffold.

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
assert guide.spacer.rna_sequence == spacer.replace("T", "U")
assert guide.pam == "TGG"
assert guide.pam_site.sequence == "TGG"
assert tuple(guide.pam_interval) == (30, 33)
assert tuple(guide.protospacer_interval) == (10, 30)
assert guide.target_strand == "-"
assert guide.cut_site().boundaries == (27,)
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
assert guide.cut_site().boundaries == (32, 37)
assert guide.cut_site().overhang == "5'"
```

## Handle ambiguous target bases

Target sequences may contain `N`, which is useful for masked or uncertain bases. A candidate is not returned when its protospacer contains an `N`, because a concrete `Spacer` and `Protospacer` accept only `A`, `C`, `G`, and `T`.

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

## Keep genome coordinates

Wrap a reference slice in [`SequenceContext`][geneops.target.SequenceContext] to preserve its assembly, contig, and coordinate offset:

```python
from geneops import Interval, SequenceContext

context = SequenceContext(
    target,
    Interval(1_000, 1_000 + len(target)),
    reference="GRCh38",
    contig="chr7",
)

guides = find_guides(context, get_nuclease("AsCas12a"))
```

Every returned guide retains this context through `guide.context` and exposes a validated `guide.target_site`.
