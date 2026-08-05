"""Make the workflow's scripts importable and expose the fixture helpers.

`workflow` has no __init__.py -- it is an implicit namespace package, which is
how the Snakefile imports it (`from workflow.scripts.common import *`). Putting
the repo root on sys.path here reproduces that, so the tests import the same
modules the workflow does.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tests import fixtures as _fixtures  # noqa: E402


@pytest.fixture
def fx():
    """The fixture-building helpers, as a module."""
    return _fixtures


@pytest.fixture
def repo_root() -> Path:
    return REPO_ROOT
