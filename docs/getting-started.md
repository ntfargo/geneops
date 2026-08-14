# Get started

Geneops requires Python 3.12 or newer. The project is currently installed from source while its public API is being established.

## Install from the repository

```bash
git clone https://github.com/ntfargo/geneops.git
cd geneops
python3 -m venv .venv
.venv/bin/python -m pip install -e .
```

On Windows, use `.venv\Scripts\python` in place of `.venv/bin/python`.

## Find your first guide

The example target contains a 20-nucleotide SpCas9 protospacer followed by a `TGG` PAM.

```python
from geneops import find_guides, get_nuclease

target = "AAAAAAAAAAATGCACGTTAGCTACGATCATGGAAAAAAAAAA"
nuclease = get_nuclease("SpCas9")

for guide in find_guides(target, nuclease):
    print(
        guide.spacer.rna_sequence,
        guide.pam_site.sequence,
        guide.pam_strand,
        guide.protospacer.interval,
        guide.cut_site().boundaries,
    )
```

`find_guides` returns [`Guide`][geneops.guide.Guide] candidates rather than bare strings. Each result distinguishes the guide [`Spacer`][geneops.target.Spacer], located [`Protospacer`][geneops.target.Protospacer], observed [`PAMSite`][geneops.target.PAMSite], nuclease, and configured nominal cut geometry needed by later workflow stages.

## Filter by PAM strand

By default both strands are searched. Restrict the search when your experiment or upstream data requires it:

```python
plus_guides = find_guides(target, nuclease, pam_strand="+")
minus_guides = find_guides(target, nuclease, pam_strand="-")
all_guides = find_guides(target, nuclease, pam_strand="both")
```

`pam_strand` identifies the strand containing the recognized PAM. It is not the strand hybridized by the guide RNA that complementary strand is available as `guide.target_strand`.

## Inspect built-in nucleases

```python
from geneops import get_nuclease, list_nucleases

print(list_nucleases())

cas12a = get_nuclease("AsCas12a")
print(cas12a.get_pam().pattern)       # TTTV
print(cas12a.get_pam().orientation)   # 5'
print(cas12a.spacer_length)           # 23
```

Continue with [guide discovery](guides/discover-guides.md), learn how to preserve [sequence context and target identity](concepts/targets.md), or review the [coordinate conventions](concepts/coordinates.md) before integrating Geneops with genomic files.
