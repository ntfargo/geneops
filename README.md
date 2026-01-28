# Geneops

Unified Python toolkit for genome editing workflows: CRISPR, base editing, and prime editing.
Design guides, score off-targets, model edit outcomes, and analyze results through a consistent.

---

## What is Geneops?

Geneops is an open-source Python library that provides shared core primitives and pipelines for
genome editing tasks, including:

- CRISPR guide RNA design
- Off-target discovery and scoring
- Nuclease and PAM constraint modeling
- Edit outcome prediction (indels, base edits, prime edits)
- QC and downstream analysis

The goal is to streamline work across laboratories and software projects by providing well-tested core primitives (targets, guides, nuclease constraints, genomic intervals, and edit outcomes), interoperable file I/O, and pluggable computational backends.

# Install
```
# Create virtual env
python -m venv venv
source venv/bin/activate

# Install geneops
pip install -e .
```

Contact: ntfargo@proton.me