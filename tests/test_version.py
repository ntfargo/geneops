"""Tests for the package version contract."""

import re

from geneops import __version__
from geneops._version import __version__ as canonical_version


def test_public_version_uses_canonical_source() -> None:
    assert __version__ == canonical_version


def test_version_is_a_pep_440_release_version() -> None:
    assert re.fullmatch(r"\d+\.\d+\.\d+", __version__)
