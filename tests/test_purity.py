"""resolve_purity: source priority and the purity_confidence QC flag."""

from __future__ import annotations

from workflow.scripts.common import PURITY_CONFIDENCE_THRESHOLD, resolve_purity


def test_known_wins_over_purecn():
    purity, source, confidence = resolve_purity(
        known=0.9, purecn_purity="0.5", purecn_failed=False, purecn_available=True
    )
    assert (purity, source, confidence) == ("0.9", "known", "high")


def test_flagged_purecn_estimate_is_used_not_rejected():
    """The 2.3.0 fix: a numeric PureCN fit is used regardless of what flagged
    it (LOW PURITY, NON-ABERRANT, ...) -- flag content only feeds
    purity_confidence, it no longer falls back to assumed_pure."""
    purity, source, confidence = resolve_purity(
        known=None, purecn_purity="0.18", purecn_failed=False, purecn_available=True
    )
    assert (purity, source) == ("0.18", "purecn")
    assert confidence == "low_purity"


def test_purecn_at_or_above_threshold_is_high_confidence():
    purity, source, confidence = resolve_purity(
        known=None,
        purecn_purity=str(PURITY_CONFIDENCE_THRESHOLD),
        purecn_failed=False,
        purecn_available=True,
    )
    assert (purity, source, confidence) == ("0.3", "purecn", "high")


def test_purecn_just_below_threshold_is_low_purity():
    purity, source, confidence = resolve_purity(
        known=None, purecn_purity="0.29", purecn_failed=False, purecn_available=True
    )
    assert (purity, source, confidence) == ("0.29", "purecn", "low_purity")


def test_p139_regression_non_aberrant_flag_no_longer_zeroes_purity():
    """P139 (WES-MCL-II): NON-ABERRANT-flagged fit at 0.31, above PureCN's own
    low-purity cutoff. Before the fix this was force-set to assumed_pure=1."""
    purity, source, confidence = resolve_purity(
        known=None, purecn_purity="0.31", purecn_failed=False, purecn_available=True
    )
    assert (purity, source, confidence) == ("0.31", "purecn", "high")


def test_purecn_failed_falls_back_to_assumed_pure():
    purity, source, confidence = resolve_purity(
        known=None, purecn_purity="0.4", purecn_failed=True, purecn_available=True
    )
    assert (purity, source, confidence) == ("1", "assumed_pure", "unknown")


def test_purecn_unavailable_falls_back_to_assumed_pure():
    purity, source, confidence = resolve_purity(
        known=None, purecn_purity="", purecn_failed=False, purecn_available=False
    )
    assert (purity, source, confidence) == ("1", "assumed_pure", "unknown")


def test_purecn_no_numeric_purity_falls_back_to_assumed_pure():
    purity, source, confidence = resolve_purity(
        known=None, purecn_purity="NA", purecn_failed=False, purecn_available=True
    )
    assert (purity, source, confidence) == ("1", "assumed_pure", "unknown")


def test_known_below_threshold_is_still_flagged_low_purity():
    purity, source, confidence = resolve_purity(
        known=0.2, purecn_purity="0.5", purecn_failed=False, purecn_available=True
    )
    assert (purity, source, confidence) == ("0.2", "known", "low_purity")
