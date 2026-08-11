# GeneOps

Unified Python toolkit for genome editing workflows: CRISPR, base editing, and prime editing.
GeneOps is building a consistent foundation for guide design, scoring, outcome modeling, and analysis.

---

## What is GeneOps?

GeneOps is an open-source Python library that provides shared core primitives and pipelines for
genome editing tasks, including:

- CRISPR guide RNA design
- Off-target discovery and scoring
- Nuclease and PAM constraint modeling
- Edit outcome prediction (indels, base edits, prime edits)
- QC and downstream analysis
 
## Install

```bash
# Create virtual env
python3 -m venv .venv

# Install GeneOps
.venv/bin/python -m pip install -e .
```

## Documentation

The documentation source lives in [`docs/`](docs/)

```bash
.venv/bin/python -m pip install -e ".[docs]"
.venv/bin/python -m mkdocs serve
```

Public running docs: https://geneops.linearfox.com/

Contributions are welcome! Please open issues or pull requests on the GitHub repository.

Contact: ntfargo@proton.me