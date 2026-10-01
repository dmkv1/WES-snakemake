"""mutect2_keep_labels / mutect2_include_expr: the FilterMutectCalls label set kept in the
final VCF and the bcftools expression built from it."""

from __future__ import annotations

import pytest

from workflow.scripts.common import (
    MUTECT2_FILTER_LABELS,
    mutect2_include_expr,
    mutect2_keep_labels,
)


def test_default_keeps_nothing_but_pass():
    assert mutect2_keep_labels({}) == []
    expr = mutect2_include_expr([])
    assert expr.count("FILTER!~") == len(MUTECT2_FILTER_LABELS)


def test_deprecated_switch_adds_strand_bias():
    labels = mutect2_keep_labels({"keep_strand_bias_calls": True, "keep_filter_labels": ["haplotype"]})
    assert labels == ["haplotype", "strand_bias"]


def test_listed_labels_are_not_excluded():
    expr = mutect2_include_expr(["clustered_events", "haplotype"])
    assert 'FILTER!~"clustered_events"' not in expr
    assert 'FILTER!~"haplotype"' not in expr
    assert 'FILTER!~"strand_bias"' in expr
    assert 'FILTER!~"germline"' in expr


def test_germline_rescue_is_bounded_by_popaf():
    expr = mutect2_include_expr(["strand_bias"], germline_min_popaf=5.0)
    assert 'FILTER!~"germline" &&' not in expr
    assert expr.endswith('(FILTER!~"germline" || INFO/POPAF[0]>=5.000000)')


@pytest.mark.parametrize("labels", [["clustred_events"], ["germline"]])
def test_rejected_labels(labels):
    with pytest.raises(ValueError):
        mutect2_keep_labels({"keep_filter_labels": labels})
