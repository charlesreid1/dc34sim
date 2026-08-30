"""Test double for the JS popsim getter and end-to-end HMM harness.

In production, `window.popsim.snapshot()` returns `{ t, f_T }` describing the
current live population. Here, we hand-set both values.

Also exposes `decode_composition()` which chains stage-1 → stage-2 → compose
into a single call — the same shape the JS `badgecestry.loadDiploid` wrapper
will have.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from hmm.compose import compose_percentages
from hmm.founder_panel import (
    ACTIVE_TYPE_MASK, FOUNDER_PANEL, INACTIVE_TYPES, N_ACTIVE_TYPES,
    N_TYPES, TYPES,
)
from hmm.full_stage import decode_haplotype
from hmm.mutation_kernel import (
    T_GRID, build_emission_tensor, snap_t_to_grid,
)
from hmm.simplified_stage import decode_diploid


@dataclass
class PopsimSnapshot:
    t: int                 # generations since founding
    f_T: np.ndarray        # (n_types,) marginal type fractions


def uniform_snapshot(t: int = 0) -> PopsimSnapshot:
    """Uniform over the ACTIVE types only (None is excluded from ancestry)."""
    f = np.where(ACTIVE_TYPE_MASK, 1.0 / N_ACTIVE_TYPES, 0.0)
    return PopsimSnapshot(t=t, f_T=f)


def biased_snapshot(fractions: dict, t: int = 0) -> PopsimSnapshot:
    """Build a snapshot from a {type_name: fraction} dict.

    Missing types get zero; the fractions are normalized to sum to 1.
    Inactive types (None) are silently dropped — the ancestry decoder does
    not model them.
    """
    f = np.zeros(N_TYPES)
    for name, val in fractions.items():
        if name in INACTIVE_TYPES:
            continue
        f[TYPES.index(name)] = val
    s = f.sum()
    if s <= 0:
        raise ValueError("fractions must sum to something positive on active types")
    f = f / s
    return PopsimSnapshot(t=t, f_T=f)


@dataclass
class DecodeReport:
    percentages: np.ndarray            # (n_types,) integer percentages summing to 100
    composition: np.ndarray            # (n_types,) real-valued composition
    per_window_h0: np.ndarray          # (n_windows, k) per-haplotype label posteriors
    per_window_h1: np.ndarray          # (n_windows, k)
    C0: np.ndarray                     # (k,) stage-1 candidates for haplo0
    C1: np.ndarray
    joint_posteriors: np.ndarray       # (n_windows, k, k, 2) stage-2 posteriors
    stage1_h0_posteriors: np.ndarray   # (n_windows, n_types) full-type per-window
    stage1_h1_posteriors: np.ndarray
    stage1_h0_alpha: np.ndarray
    stage1_h0_beta: np.ndarray
    stage1_h1_alpha: np.ndarray
    stage1_h1_beta: np.ndarray
    t_used: int                        # grid-snapped t
    rate_name: str


def decode_composition(
    diploid: np.ndarray,
    snapshot: PopsimSnapshot,
    rate_name: str = "Baseline",
    top_k: int = 3,
    min_threshold: float = 0.05,
    emission_tensor: Optional[np.ndarray] = None,
) -> DecodeReport:
    """End-to-end two-stage decode with claim-10 display rounding.

    diploid          : (2, 9) uint8
    snapshot         : PopsimSnapshot with (t, f_T)
    rate_name        : mutation rate the population has been evolving under
    top_k            : stage-1 → stage-2 pruning size (default 3, per §5.5)
    min_threshold    : claim-10 min mass to survive display rounding
    emission_tensor  : optionally precomputed (n_grid, n_types, n_loci, n_alleles);
                       one is built if not provided (rebuilding is O(seconds)
                       but this lets callers cache across many diploids).
    """
    if emission_tensor is None:
        emission_tensor = build_emission_tensor(FOUNDER_PANEL, rate_name)
    t_idx = snap_t_to_grid(snapshot.t)
    E_TLA = emission_tensor[t_idx]

    r0 = decode_haplotype(diploid[0], E_TLA, snapshot.f_T, top_k=top_k)
    r1 = decode_haplotype(diploid[1], E_TLA, snapshot.f_T, top_k=top_k)

    res2 = decode_diploid(
        diploid, E_TLA, snapshot.f_T,
        C0=r0.top_k_indices, C1=r1.top_k_indices,
    )

    return DecodeReport(
        percentages=compose_percentages(res2.composition, min_threshold),
        composition=res2.composition,
        per_window_h0=res2.joint_posteriors.sum(axis=(2, 3)),
        per_window_h1=res2.joint_posteriors.sum(axis=(1, 3)),
        C0=r0.top_k_indices, C1=r1.top_k_indices,
        joint_posteriors=res2.joint_posteriors,
        stage1_h0_posteriors=r0.posteriors,
        stage1_h1_posteriors=r1.posteriors,
        stage1_h0_alpha=r0.alpha, stage1_h0_beta=r0.beta,
        stage1_h1_alpha=r1.alpha, stage1_h1_beta=r1.beta,
        t_used=T_GRID[t_idx],
        rate_name=rate_name,
    )
