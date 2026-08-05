"""The vendored parser must stay identical to its upstream original.

Two copies of a parser whose failure mode is silently returning the wrong field
would drift invisibly -- both keep working, each on a slightly different set of
headers. Comparing the sources object by object turns that into a failing test.

Skipped where upstream is absent, which is every ported deployment. Point
WESINGEST_PATH at the file to run it elsewhere.
"""

from __future__ import annotations

import importlib.util
import inspect
import os
import sys
from pathlib import Path

import pytest

from workflow.scripts import fastq_header

UPSTREAM = Path(os.environ.get(
    "WESINGEST_PATH", "/mnt/data/NGS/WES/scripts/wesingest/header.py"))

VENDORED = ["ReadHeader", "parse_header", "read_first_header"]


@pytest.fixture(scope="module")
def upstream():
    if not UPSTREAM.is_file():
        pytest.skip(f"upstream parser not present at {UPSTREAM}")
    spec = importlib.util.spec_from_file_location("_wesingest_header", UPSTREAM)
    module = importlib.util.module_from_spec(spec)
    # inspect.getsource resolves a class through sys.modules, so the dynamically
    # loaded module has to be registered before its classes can be read back.
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        del sys.modules[spec.name]


@pytest.mark.parametrize("name", VENDORED)
def test_vendored_object_matches_upstream(upstream, name):
    ours = inspect.getsource(getattr(fastq_header, name))
    theirs = inspect.getsource(getattr(upstream, name))
    assert ours == theirs, (
        f"{name} has diverged from {UPSTREAM}. The upstream copy is canonical: "
        f"re-copy it here rather than editing this file."
    )


def test_header_regex_matches_upstream(upstream):
    assert fastq_header.HEADER_RE.pattern == upstream.HEADER_RE.pattern


def test_pipeline_local_additions_are_not_upstream(upstream):
    """parse_index is ours until docs/TODO.md item 2 lands upstream."""
    assert not hasattr(upstream, "parse_index")
