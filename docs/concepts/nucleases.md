# Nucleases and PAMs

A [`Nuclease`][geneops.nuclease.Nuclease] combines the constraints needed to interpret a candidate editing site: a name, explicit PAM definition, spacer length, and explicit strand-cleavage pattern.

## Built-in registry

Geneops currently includes these definitions:

| Family | Nuclease | PAM | Spacer | Orientation |
| --- | --- | --- | ---: | :---: |
| Cas9 | SpCas9 | `NGG` | 20 nt | 3′ |
| Cas9 | nCas9 | `NGG` | 20 nt | 3′ |
| Cas9 | SpCas9-HF1 | `NGG` | 20 nt | 3′ |
| Cas9 | eSpCas9 | `NGG` | 20 nt | 3′ |
| Cas9 | SpCas9-NG | `NG` | 20 nt | 3′ |
| Cas9 | SpG | `NGN` | 20 nt | 3′ |
| Cas9 | SpRY | `NRN` | 20 nt | 3′ |
| Cas9 | SaCas9 | `NNGRRT` | 21 nt | 3′ |
| Cas9 | NmCas9 | `NNNNGATT` | 24 nt | 3′ |
| Cas9 | CjCas9 | `NNNNRYAC` | 22 nt | 3′ |
| Cas12a | AsCas12a | `TTTV` | 23 nt | 5′ |
| Cas12a | LbCas12a | `TTTV` | 23 nt | 5′ |

Lookups are case-insensitive:

```python
from geneops import get_nuclease

assert get_nuclease("spcas9").name == "SpCas9"
```

## IUPAC-aware PAM matching

PAM definitions use the IUPAC DNA alphabet. For example, `N` accepts any canonical DNA base, `R` accepts A or G, and `V` accepts A, C, or G.

```python
spcas9 = get_nuclease("SpCas9")

spcas9.matches_pam("AGG")  # True
spcas9.matches_pam("TGG")  # True
spcas9.matches_pam("AGA")  # False
```

## Explicit nominal cleavage geometry

[`CleavagePattern`][geneops.nuclease.CleavagePattern] stores independent cut
configured offsets for the PAM strand and target strand. Both are boundary offsets measured
from the protospacer's 5′ end in PAM-strand orientation.

- Equal offsets describe a blunt cut.
- Different offsets describe a staggered cut their absolute difference is the overhang length.
- `None` on one strand describes a nickase.

The built-in SpCas9 pattern is the canonical `(17, 17)` geometry for its 20-nt
protospacer. The built-in Cas12a patterns use canonical `(18, 23)` boundaries.
These values live in each registry entry Geneops does not infer them from names
such as `Cas9`, `Cas12`, or `Cpf1`.