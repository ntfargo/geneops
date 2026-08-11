# Nucleases and PAMs

A [`Nuclease`][geneops.nuclease.Nuclease] combines the constraints needed to interpret a candidate editing site: a name, PAM pattern, spacer length, PAM orientation, and whether the enzyme is a nickase.

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

## Geometry assumptions

The current cut model provides conventional family-level geometry:

- Cas9-family definitions use a 3′ PAM and a blunt cut three bases upstream of the PAM.
- Cas12-family definitions use a 5′ PAM and return staggered boundaries with a five-nucleotide overhang.
- Nickases report single-strand behavior through `is_nickase` and `nicking_offset`.

!!! warning "Model scope"
    These are compact workflow defaults, not experimental-condition-specific cleavage models. Validate assumptions against the enzyme variant and assay used in your experiment.
