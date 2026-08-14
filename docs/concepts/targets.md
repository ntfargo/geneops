# Sequence context and target sites

A guide sequence alone does not identify a genome-editing target. The same sequence can occur on different contigs, assemblies, or strands, and many scoring and outcome models also require bases flanking the protospacer and PAM.

GeneOps separates the related biological objects instead of treating every
sequence as a generic string:

- [`Spacer`][geneops.target.Spacer] is the location-free targeting portion of
  a guide RNA construct.
- [`Protospacer`][geneops.target.Protospacer] is the reference-located DNA
  sequence targeted by a spacer.
- [`PAMSite`][geneops.target.PAMSite] is an observed, located sequence that
  matches a nuclease's PAM definition.
- [`SequenceContext`][geneops.target.SequenceContext] is a forward-reference sequence slice with an interval and optional reference and contig identifiers.
- [`TargetSite`][geneops.target.TargetSite] binds validated protospacer and PAM intervals to that context and a nuclease.
- [`Guide`][geneops.guide.Guide] is a guide candidate connecting a spacer to a
  target site. It does not represent the full scaffold-containing guide RNA.

## Spacer is not protospacer

The terms are intentionally distinct. A spacer belongs to the guide construct
and has no genomic coordinate. A protospacer belongs to a reference sequence
the same spacer can match multiple genomic protospacers.

```python
guide.spacer.sequence        # DNA spelling used for comparison
guide.spacer.rna_sequence    # RNA spelling used in the guide
guide.protospacer.sequence   # DNA sequence at this target
guide.protospacer.interval   # forward-reference location
guide.pam_site.sequence      # observed PAM, such as TGG
guide.pam_site.pam.pattern   # accepted pattern, such as NGG
```

For the DNA-targeting nucleases currently supported, spacer and protospacer
have the same nucleotide spelling. They remain separate objects because one is
a reagent component and the other is a reference-located biological feature.
This follows the terminology used by the
[crisprVerse domain model](https://github.com/crisprVerse/crisprBase) and keeps
location separate from literal sequence as recommended by
[GA4GH VRS](https://vrs.ga4gh.org/en/latest/concepts/LocationAndReference/SequenceLocation.html).

## Preserve reference coordinates

```python
from geneops import Interval, SequenceContext, find_guides, get_nuclease

spacer = "ATGCACGTTAGCTACGATCA"
sequence = "A" * 10 + spacer + "TGG" + "A" * 10

context = SequenceContext(
    sequence=sequence,
    interval=Interval(1_000, 1_000 + len(sequence)),
    reference="GRCh38",
    contig="chr7",
)

guide = find_guides(
    context,
    get_nuclease("SpCas9"),
    pam_strand="+",
)[0]

guide.pam_interval          # Interval(start=1030, end=1033)
guide.protospacer_interval  # Interval(start=1010, end=1030)
guide.context.reference     # "GRCh38"
guide.spacer.rna_sequence   # "AUGCACGUUAGCUACGAUCA"
```

Passing a plain string remains supported. GeneOps then creates an anonymous `SequenceContext` spanning `[0, len(sequence))`.

## Request sequence windows for models

Sequence models use different input windows, so `SequenceContext` does not hard-code one length. Request exactly the strand-oriented flanks a model needs:

```python
site = guide.target_site
window = site.sequence_context(upstream=4, downstream=3)

assert len(window) == 30  # 4 + 20 nt spacer + 3 nt PAM + 3
```

This 30-nucleotide layout is used by models such as [Rule Set 3](https://www.nature.com/articles/s41467-022-33024-2) and [CRISPRon](https://www.nature.com/articles/s41467-021-23576-0). The sequence is only one model input: for example, Rule Set 3 also uses tracrRNA identity, while CRISPRon includes binding energy. Repair-outcome models such as [inDelphi](https://www.nature.com/articles/s41586-018-0686-x) also depend on local target sequence. Keeping the window operation explicit lets future scoring and outcome modules declare all of their requirements without changing target identity.

## One coordinate contract

The context sequence is always stored in forward-reference orientation. `TargetSite.protospacer_sequence`, `TargetSite.pam_sequence`, and `TargetSite.sequence_context()` return sequence 5′ to 3′ in the PAM-strand orientation. Intervals never change orientation.

Reference and contig identifiers are optional strings rather than assumptions about a particular genome provider. This leaves room for future integration with accession- and digest-based identifiers while keeping small synthetic targets easy to use.
