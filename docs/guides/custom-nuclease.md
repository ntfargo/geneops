# Define a custom nuclease

Use a custom [`Nuclease`][geneops.nuclease.Nuclease] when a variant is not present in the built-in registry.

## Construct and validate

```python
from geneops import Nuclease, validate_nuclease
from geneops.nuclease import PAM

custom = Nuclease(
    name="MyCas9",
    pam=PAM("NGA", "3'"),
    spacer_length=20,
    description="An experimental Cas9-like variant",
)

validate_nuclease(custom)
```

Using a `PAM` object makes the orientation explicit. A plain string is also accepted Geneops currently infers 5′ orientation for names containing `Cas12` or `Cpf1`, and 3′ orientation otherwise.

!!! tip
    Prefer an explicit `PAM` object for custom enzymes. Name-based inference is convenient for common families but should not encode novel enzyme biology.

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
