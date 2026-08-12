# API reference

The API reference is generated from Geneops source code and docstrings by `mkdocstrings`. It therefore stays synchronized with the installed package at documentation build time.

## Public package surface

Most workflows can import these objects directly from `geneops`:

```python
from geneops import (
    Guide,
    Interval,
    Nuclease,
    find_guides,
    get_nuclease,
)
```

Module-specific helpers are documented in the following sections:

- [Coordinates](coordinates.md) intervals and strand operations
- [Guides](guide.md) sequence normalization, discovery, and guide objects
- [Nucleases](nuclease.md) PAM definitions, registries, and cut geometry
