# Define a custom nuclease

Use a custom [`Nuclease`][geneops.nuclease.Nuclease] when a variant is not present in the built-in registry.

## Construct and validate

```python
from geneops import CleavagePattern, Nuclease, PAM, validate_nuclease

custom = Nuclease(
    name="MyEditor",
    pam=PAM("NGA", "3'"),
    cleavage=CleavagePattern(
        pam_strand_offset=17,
        target_strand_offset=17,
    ),
    spacer_length=20,
    description="An experimental editing nuclease",
)

validate_nuclease(custom)
```

Both targeting and nominal cleavage geometry are explicit. Cleavage offsets are
zero-based boundaries measured from the protospacer's 5′ end in PAM-strand
orientation. Use `None` for a strand that is not cut. Choose offsets from an
appropriate experimental source and document when they are canonical rather
than exact measurements for the target being analyzed.

For example, a target-strand nickase uses a one-strand pattern:

```python
nickase = Nuclease(
    name="MyNickase",
    pam=PAM("NGG", "3'"),
    cleavage=CleavagePattern(
        pam_strand_offset=None,
        target_strand_offset=17,
    ),
)
```

## Use it directly

Registration is optional. Any API accepting a `Nuclease` can use the object directly:

```python
from geneops import find_guides

guides = find_guides(target_sequence, custom)
```

## Add it to a registry

Use an independent registry when a workflow needs its own controlled catalog:

```python
from geneops import RRegistry

registry = RRegistry([custom])
assert registry.get_nuclease("mycas9") is custom
```

Or register it globally for the current Python process:

```python
from geneops import get_nuclease, register_nuclease

register_nuclease(custom)
assert get_nuclease("MyCas9") is custom
```

Duplicate normalized names raise `ValueError`, and unknown names raise `KeyError`.
