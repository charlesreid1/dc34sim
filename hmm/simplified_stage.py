"""Stage-2 (simplified, joint) diploid HMM over C0 × C1 × {switch}.

Runs forward-backward on the pruned joint state space (§5.5). Emission per
window is the product of per-locus factors on whichever haploid is currently
assigned to each parent-of-origin label (§5.3). For W_DE, an additional hard
constraint tying `phenotype.nonlin = sat_add(haplo0.chaser, haplo1.nonlin)`
is baked into the emission via a delta factor (§5.6).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hmm.founder_panel import LOCUS
from hmm.full_stage import N_WINDOWS, WINDOW_LOCI, WINDOW_NAMES

# Index of W_DE (the merged D+E window that carries the nonlin phenotype).
W_DE_INDEX = WINDOW_NAMES.index("W_DE")
CHASER = LOCUS["chaser"]
NONLIN = LOCUS["nonlin"]


def sat_add(x: int, y: int) -> int:
    return min(255, int(x) + int(y))


def phenotype_nonlin(diploid: np.ndarray) -> int:
    """Compute the shipped `phenotype().nonlin` per §4 of pop-gen doc."""
    return sat_add(int(diploid[0, CHASER]), int(diploid[1, NONLIN]))


@dataclass
class SimplifiedStageResult:
    """Everything the internals panel needs from the stage-2 decode."""
    C0: np.ndarray                  # (k,) type indices for haplo0 candidates
    C1: np.ndarray                  # (k,) type indices for haplo1 candidates
    emissions: np.ndarray           # (n_windows, k, k, 2)  E^{(2)}(W)
    alpha: np.ndarray               # (n_windows, k, k, 2)
    beta: np.ndarray                # (n_windows, k, k, 2)
    joint_posteriors: np.ndarray    # (n_windows, k, k, 2)  γ_i(T0, T1, switch)
    composition: np.ndarray         # (n_types,) summed mass over windows & haploids
    log_likelihood: float


def _stage2_emission(
    diploid: np.ndarray,
    E_TLA: np.ndarray,
    C0: np.ndarray,
    C1: np.ndarray,
    use_nonlin_bug: bool = True,
) -> np.ndarray:
    """Compute the (n_windows, |C0|, |C1|, 2) joint-state emission tensor.

    diploid          : (2, 9) uint8
    E_TLA            : (n_types, n_loci, n_alleles)
    C0, C1           : candidate type indices for each haplotype
    use_nonlin_bug   : if True, apply the shipped δ(phenotype.nonlin = ...)
                       factor on W_DE (§5.6). If False, drop it (used by tests
                       that exercise the bug's asymmetry).
    """
    k0, k1 = len(C0), len(C1)
    E = np.zeros((N_WINDOWS, k0, k1, 2), dtype=np.float64)

    # Precompute per-window emissions per (type, switch) for each haplotype's
    # assignment. Under switch=0: T0 emits haplo0, T1 emits haplo1.
    # Under switch=1: T0 emits haplo1, T1 emits haplo0.
    # For efficiency we compute e0[T, w] = ∏ over ℓ ∈ w of E[T, ℓ, haplo0[ℓ]]
    # and same for e1, then combine.
    n_types = E_TLA.shape[0]
    e0 = np.ones((n_types, N_WINDOWS), dtype=np.float64)
    e1 = np.ones((n_types, N_WINDOWS), dtype=np.float64)
    for w_idx, loci in enumerate(WINDOW_LOCI):
        for ell in loci:
            e0[:, w_idx] *= E_TLA[:, ell, int(diploid[0, ell])]
            e1[:, w_idx] *= E_TLA[:, ell, int(diploid[1, ell])]

    for w in range(N_WINDOWS):
        for i, T0 in enumerate(C0):
            for j, T1 in enumerate(C1):
                # switch=0: T0 sees haplo0, T1 sees haplo1
                E[w, i, j, 0] = e0[T0, w] * e1[T1, w]
                # switch=1: T0 sees haplo1, T1 sees haplo0
                E[w, i, j, 1] = e1[T0, w] * e0[T1, w]

    if use_nonlin_bug:
        # W_DE δ factor: observed phenotype.nonlin must equal what the state's
        # switch-assigned slots produce under the shipped `sat_add(slot0.chaser,
        # slot1.nonlin)`. switch=0 always matches by construction (it's how
        # the phenotype was computed). switch=1 matches iff swapping the two
        # implied slots preserves the value — which is what makes the bug
        # asymmetric and stage-2 informative.
        observed_pnonlin = phenotype_nonlin(diploid)
        a0_chaser = int(diploid[0, CHASER])
        a1_chaser = int(diploid[1, CHASER])
        a0_nonlin = int(diploid[0, NONLIN])
        a1_nonlin = int(diploid[1, NONLIN])
        implied_switch0 = sat_add(a0_chaser, a1_nonlin)
        implied_switch1 = sat_add(a1_chaser, a0_nonlin)
        if implied_switch0 != observed_pnonlin:
            E[W_DE_INDEX, :, :, 0] = 0.0
        if implied_switch1 != observed_pnonlin:
            E[W_DE_INDEX, :, :, 1] = 0.0

    return E


def _joint_forward_backward(
    emissions: np.ndarray,
    f_T: np.ndarray,
    C0: np.ndarray,
    C1: np.ndarray,
) -> tuple:
    """Forward-backward on the joint (T0, T1, switch) state space.

    emissions : (n_windows, k, k, 2)
    f_T       : (n_types,) type prior
    C0, C1    : candidate index arrays

    Transitions factor: T0' | T0 ~ f_{T0'} restricted to C0 and renormalized;
    same for T1'; switch flips with P(switch)=1/2. So the joint transition
    is the outer product of three independent 1D distributions.
    """
    n_windows = emissions.shape[0]
    k0, k1 = len(C0), len(C1)

    # Restricted, renormalized type priors on candidate sets.
    p_T0 = f_T[C0].astype(np.float64)
    if p_T0.sum() <= 0:
        p_T0 = np.ones(k0) / k0
    else:
        p_T0 = p_T0 / p_T0.sum()

    p_T1 = f_T[C1].astype(np.float64)
    if p_T1.sum() <= 0:
        p_T1 = np.ones(k1) / k1
    else:
        p_T1 = p_T1 / p_T1.sum()

    p_switch = np.array([0.5, 0.5])

    # Initial prior over (i, j, s): p_T0[i] * p_T1[j] * p_switch[s]
    prior = np.einsum("i,j,s->ijs", p_T0, p_T1, p_switch)

    em = emissions.copy()
    for w in range(n_windows):
        if em[w].sum() <= 0:
            em[w] = 1.0  # zero-emission fallback (same rationale as stage-1)

    alpha = np.zeros_like(em)
    beta = np.zeros_like(em)
    scales = np.zeros(n_windows, dtype=np.float64)

    a0 = prior * em[0]
    scales[0] = a0.sum()
    alpha[0] = a0 / scales[0]

    for w in range(1, n_windows):
        # Because the joint transition factors, message propagation reduces to
        #   next[i', j', s'] = p_T0[i'] * p_T1[j'] * p_switch[s'] * total_mass
        # where total_mass = sum over (i, j, s) of alpha[w-1, i, j, s] = 1 (scaled).
        # So the "prior-only" contribution is `prior` itself, and the emission
        # multiplies in. That's a factored transition of the same shape as the
        # initial step, which is what we'd expect from the rank-1 stage-1.
        aw = prior * em[w]
        s = aw.sum()
        scales[w] = s
        alpha[w] = aw / s

    beta[-1] = 1.0
    beta[-1] /= beta[-1].sum()
    for w in range(n_windows - 2, -1, -1):
        # β_w = (Σ over next state of prior_next * em_{w+1} * β_{w+1}) then
        # projected back — but since the transition is factored and rank-1
        # in each factor, the sum over next state (i', j', s') is just
        # sum(prior * em[w+1] * beta[w+1]) which is a scalar. So beta[w] is
        # uniform over states (proportional to 1), scaled.
        rhs = prior * em[w + 1] * beta[w + 1]
        total = rhs.sum()
        beta[w] = np.ones_like(beta[w]) * total
        beta[w] /= beta[w].sum()

    log_L = float(np.log(scales).sum())
    return alpha, beta, log_L


def decode_diploid(
    diploid: np.ndarray,
    E_TLA: np.ndarray,
    f_T: np.ndarray,
    C0: np.ndarray,
    C1: np.ndarray,
    use_nonlin_bug: bool = True,
) -> SimplifiedStageResult:
    """Full stage-2 decode on a diploid, restricted to the given candidate sets.

    diploid : (2, 9) uint8
    E_TLA   : (n_types, n_loci, n_alleles) emissions at current t
    f_T     : (n_types,) type marginals
    C0, C1  : (k,) candidate type indices from stage-1 (§5.5)
    """
    if diploid.shape != (2, 9):
        raise ValueError(f"diploid must have shape (2, 9); got {diploid.shape}")

    emissions = _stage2_emission(diploid, E_TLA, C0, C1, use_nonlin_bug)
    alpha, beta, log_L = _joint_forward_backward(emissions, f_T, C0, C1)

    joint = alpha * beta
    sums = joint.sum(axis=(1, 2, 3), keepdims=True)
    joint = np.where(sums > 0, joint / sums, joint)

    # Composition: for each window, marginalize over switch and the "other"
    # haplotype's label to get per-haplotype label posteriors, sum across
    # windows and haplotypes.
    n_types = E_TLA.shape[0]
    composition = np.zeros(n_types, dtype=np.float64)
    # γ_w(T0 = t) summed over (j, s) -> P(haplo0 lineage = t at window w).
    p_h0 = joint.sum(axis=(2, 3))  # (n_windows, k0)
    p_h1 = joint.sum(axis=(1, 3))  # (n_windows, k1)
    for w in range(joint.shape[0]):
        for i, T in enumerate(C0):
            composition[T] += p_h0[w, i]
        for j, T in enumerate(C1):
            composition[T] += p_h1[w, j]
    composition /= composition.sum()

    return SimplifiedStageResult(
        C0=C0, C1=C1,
        emissions=emissions,
        alpha=alpha, beta=beta,
        joint_posteriors=joint,
        composition=composition,
        log_likelihood=log_L,
    )
