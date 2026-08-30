"""Claim-10 display rounding: min-threshold drop, renormalize, round to whole %.

Given a composition vector (nonneg, sums to 1), apply the patent's display rule:
  1. Zero out any type with mass < min_threshold.
  2. Renormalize the survivors so they sum to 1.
  3. Round to integer percentages that sum to exactly 100 (largest-remainder
     method — the "delta rounding" referenced in claim 10).
"""

from __future__ import annotations

import numpy as np

DEFAULT_MIN_THRESHOLD = 0.05  # 5%


def compose_percentages(
    composition: np.ndarray,
    min_threshold: float = DEFAULT_MIN_THRESHOLD,
) -> np.ndarray:
    """Return integer percentages summing to exactly 100."""
    if composition.ndim != 1:
        raise ValueError(f"composition must be 1D; got shape {composition.shape}")
    if (composition < 0).any():
        raise ValueError("composition entries must be nonneg")
    total = composition.sum()
    if total <= 0:
        raise ValueError("composition sum must be positive")

    # Normalize to a proper distribution first.
    p = composition / total
    # Drop sub-threshold types.
    survivors_mask = p >= min_threshold
    if not survivors_mask.any():
        # Everything is below threshold; fall back to keeping the argmax.
        survivors_mask = np.zeros_like(p, dtype=bool)
        survivors_mask[int(p.argmax())] = True

    # Renormalize survivors to sum to 1.
    survivors = np.where(survivors_mask, p, 0.0)
    survivors = survivors / survivors.sum()

    # Largest-remainder rounding to 100.
    raw = survivors * 100.0
    floors = np.floor(raw).astype(np.int64)
    remainders = raw - floors
    short = 100 - int(floors.sum())
    if short > 0:
        # Give the extra points to the largest remainders.
        # Tie-break on original composition value (stable, deterministic).
        order = np.argsort(-remainders, kind="stable")
        for k in order[:short]:
            floors[k] += 1
    elif short < 0:
        # Overshoot (rare with integer arithmetic here, but be defensive).
        order = np.argsort(remainders, kind="stable")
        for k in order[:-short]:
            floors[k] -= 1
    return floors
