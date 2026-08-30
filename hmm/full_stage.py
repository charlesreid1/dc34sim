"""Stage-1 (full-type) per-haplotype HMM.

Runs forward-backward on the 8-state (all badge types) HMM across the four
windows on one haploid observation. Returns α, β, per-window posteriors, and
the top-k candidate types by summed posterior mass (§5.5 of plan).

Transition matrix is rank-1: P(T' | T) = f_{T'}(t) (§5.4). Emission per
(T, W) factorizes as ∏_{ℓ ∈ W} E_{T, ℓ, t}(a_ℓ).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hmm.founder_panel import ACTIVE_TYPE_MASK, LOCUS, N_TYPES

# The four windows of §5.1. W_DE merges linkage groups D and E.
WINDOWS = (
    ("W_A", (LOCUS["cd_period"], LOCUS["cd_rate"], LOCUS["cd_dir"])),
    ("W_B", (LOCUS["sat"],)),
    ("W_C", (LOCUS["hue_ratedir"], LOCUS["hue_base"], LOCUS["hue_bound"])),
    ("W_DE", (LOCUS["chaser"], LOCUS["nonlin"])),
)
WINDOW_NAMES = tuple(w[0] for w in WINDOWS)
WINDOW_LOCI = tuple(w[1] for w in WINDOWS)
N_WINDOWS = len(WINDOWS)


@dataclass
class FullStageResult:
    """Everything the internals panel needs from the stage-1 decode."""
    emissions: np.ndarray            # (n_windows, n_types)  E_T(W_i)
    alpha: np.ndarray                # (n_windows, n_types)  forward messages
    beta: np.ndarray                 # (n_windows, n_types)  backward messages
    posteriors: np.ndarray           # (n_windows, n_types)  γ_i(T)
    summed_posterior: np.ndarray     # (n_types,)           Σ_i γ_i(T)
    top_k_indices: np.ndarray        # (k,) type indices, sorted by mass desc
    log_likelihood: float            # log P(observation | model)


def per_window_emission(haploid: np.ndarray, E_TLA: np.ndarray) -> np.ndarray:
    """Compute E_T(W_i) = ∏_{ℓ ∈ W_i} E[T, ℓ, a_ℓ] for one haploid.

    haploid : (9,) uint8 alleles.
    E_TLA   : (n_types, n_loci, n_alleles) emission tensor for the current t.
    returns : (n_windows, n_types) emission likelihoods.
    """
    out = np.ones((N_WINDOWS, N_TYPES), dtype=np.float64)
    for w_idx, loci in enumerate(WINDOW_LOCI):
        for ell in loci:
            a = int(haploid[ell])
            out[w_idx] *= E_TLA[:, ell, a]
    return out


def forward_backward(emissions: np.ndarray, f_T: np.ndarray) -> tuple:
    """Run α/β on the rank-1 HMM.

    emissions : (n_windows, n_types)
    f_T       : (n_types,) marginal prior per type

    Because transitions are rank-1 (P(T' | T) = f_{T'} for all T), forward
    scales without mixing across states in the usual way: the message after
    a rank-1 transition step is `f_T * (α · 1) = f_T * total_mass`. That
    makes α_i(T) = f_T · E_i(T) · (Π_{j < i} <α_j, 1>), with the message
    total absorbed into the normalizer. Same for β symmetrically.

    Rather than special-case rank-1, we just run the general forward-backward
    with the explicit T×T transition matrix `A = 1 ⊗ f_T`, in scaled form.
    This is O(n_windows · n_types^2), which at 4 × 64 = 256 mults per pass
    is completely negligible and keeps the code shape identical to what a
    later non-rank-1 transition would need.

    Returns (alpha, beta, log_likelihood) with alpha/beta scaled so each row
    sums to 1 (Rabiner-style scaling); log_likelihood is the sum of log
    scaling factors from the forward pass.
    """
    n_windows, n_types = emissions.shape
    A = np.broadcast_to(f_T, (n_types, n_types)).copy()  # A[T, T'] = f_{T'}

    # Zero-emission fallback: if a window has no type that can generate the
    # observation (common at t=0 when the observation falls outside every
    # founder support), that window's emission is uninformative — fall back
    # to a uniform row so forward-backward stays defined and the observation
    # simply contributes no likelihood ratio to the decode. This matches the
    # generative semantics "this window is unexplained under every type."
    em = emissions.copy()
    for i in range(n_windows):
        if em[i].sum() <= 0:
            em[i] = 1.0

    alpha = np.zeros((n_windows, n_types), dtype=np.float64)
    scales = np.zeros(n_windows, dtype=np.float64)

    # Initial: α_0(T) ∝ f_T · E_0(T).
    a0 = f_T * em[0]
    scales[0] = a0.sum()
    alpha[0] = a0 / scales[0]

    for i in range(1, n_windows):
        # a_i(T') = E_i(T') · Σ_T α_{i-1}(T) · A[T, T']
        prev = alpha[i - 1] @ A         # (n_types,)
        ai = em[i] * prev
        s = ai.sum()
        scales[i] = s
        alpha[i] = ai / s

    beta = np.zeros_like(alpha)
    beta[-1] = 1.0
    # Scale beta[-1] like the others so shapes match.
    beta[-1] /= beta[-1].sum()
    for i in range(n_windows - 2, -1, -1):
        rhs = em[i + 1] * beta[i + 1]
        bi = A @ rhs
        s = bi.sum()
        beta[i] = bi / s

    log_likelihood = float(np.log(scales).sum())
    return alpha, beta, log_likelihood


def decode_haplotype(
    haploid: np.ndarray,
    E_TLA: np.ndarray,
    f_T: np.ndarray,
    top_k: int = 3,
) -> FullStageResult:
    """Full stage-1 decode of one haploid.

    haploid : (9,) uint8
    E_TLA   : (n_types, n_loci, n_alleles) emissions at current t
    f_T     : (n_types,) type marginals
    top_k   : keep this many types by summed posterior mass (default 3, §5.5)
    """
    if f_T.shape != (N_TYPES,):
        raise ValueError(f"f_T must have shape ({N_TYPES},); got {f_T.shape}")
    if not np.isclose(f_T.sum(), 1.0):
        raise ValueError(f"f_T must sum to 1; got {f_T.sum()}")

    # Enforce the ancestry-decoder's exclusion of inactive badge types (None).
    # Zeroing the prior propagates through forward-backward: any state with
    # zero prior mass stays zero everywhere.
    f_T = np.where(ACTIVE_TYPE_MASK, f_T, 0.0)
    s = f_T.sum()
    if s <= 0:
        raise ValueError("f_T has zero mass on active types")
    f_T = f_T / s

    emissions = per_window_emission(haploid, E_TLA)
    alpha, beta, log_L = forward_backward(emissions, f_T)

    posteriors = alpha * beta
    posteriors /= posteriors.sum(axis=1, keepdims=True)

    summed = posteriors.sum(axis=0)
    # Restrict top-k selection to active types only; inactive types get -inf
    # score so they never appear in candidate sets.
    ranked = np.where(ACTIVE_TYPE_MASK, summed, -np.inf)
    top_k_indices = np.argsort(-ranked)[:top_k]

    return FullStageResult(
        emissions=emissions,
        alpha=alpha,
        beta=beta,
        posteriors=posteriors,
        summed_posterior=summed,
        top_k_indices=top_k_indices,
        log_likelihood=log_L,
    )
