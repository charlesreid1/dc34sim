import numpy as np
import pytest

from hmm.founder_panel import FOUNDER_PANEL, N_TYPES, TYPE_INDEX
from hmm.full_stage import (
    N_WINDOWS, WINDOW_LOCI,
    decode_haplotype, forward_backward, per_window_emission,
)
from hmm.mutation_kernel import build_emission_tensor, T_GRID


def _uniform_f_T():
    """Uniform over the 7 active types (None excluded). Matches what the
    decoder does internally when handed a naive-uniform input."""
    from hmm.founder_panel import ACTIVE_TYPE_MASK, N_ACTIVE_TYPES
    return np.where(ACTIVE_TYPE_MASK, 1.0 / N_ACTIVE_TYPES, 0.0)


def _sample_haploid_from_type(T: str, seed: int = 0) -> np.ndarray:
    """Draw a founder haploid from `Haploid::from_type(T)` distribution."""
    rng = np.random.default_rng(seed)
    ti = TYPE_INDEX[T]
    pi = FOUNDER_PANEL[ti]  # (9, 256)
    out = np.zeros(9, dtype=np.uint8)
    for ell in range(9):
        out[ell] = rng.choice(256, p=pi[ell])
    return out


def test_per_window_emission_shape():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = np.zeros(9, dtype=np.uint8)
    em = per_window_emission(haploid, E[0])
    assert em.shape == (N_WINDOWS, N_TYPES)


def test_forward_backward_alpha_beta_shape_and_normalization():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = _sample_haploid_from_type("Uber", seed=1)
    em = per_window_emission(haploid, E[0])
    alpha, beta, log_L = forward_backward(em, _uniform_f_T())
    assert alpha.shape == (N_WINDOWS, N_TYPES)
    assert beta.shape == (N_WINDOWS, N_TYPES)
    # Rows sum to 1 by scaling convention.
    assert np.allclose(alpha.sum(axis=1), 1.0)
    # Beta[-1] initialized to 1 (unscaled), then scaled by row sum.
    assert np.allclose(beta.sum(axis=1), 1.0)
    assert np.isfinite(log_L)


@pytest.mark.parametrize("T", ["Uber", "Goon", "Village", "CtfContest"])
def test_founder_haploid_top1_is_own_type_at_t0(T):
    """A haploid sampled from Haploid::from_type(T) at t=0 should decode to T
    as the top-posterior state under a uniform f_T prior. Skipped for Human
    and None because their founder distributions are near-identical (differ
    only in cd_period_max), so which one comes out on top depends on the
    per-locus draw."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = _sample_haploid_from_type(T, seed=42)
    res = decode_haplotype(haploid, E[0], _uniform_f_T())
    assert res.top_k_indices[0] == TYPE_INDEX[T], (
        f"expected top1={T}, got summed={res.summed_posterior}, "
        f"top_k={res.top_k_indices}"
    )


def test_none_never_appears_in_candidates():
    """None is inactive in the ancestry decoder — it must never be in the
    top-k candidate set for either haplotype, regardless of the observation."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    for seed in range(20):
        # Sample a None-drawn haploid (worst case for exclusion): even then,
        # top-k must not contain None.
        haploid = _sample_haploid_from_type("None", seed=seed)
        # top_k=7 = all active types; None still must not appear.
        res = decode_haplotype(haploid, E[0], _uniform_f_T(), top_k=7)
        assert TYPE_INDEX["None"] not in res.top_k_indices, (
            f"seed={seed}: None leaked into candidates: {res.top_k_indices}"
        )


def test_none_gets_zero_posterior_mass():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = _sample_haploid_from_type("Human", seed=1)
    res = decode_haplotype(haploid, E[0], _uniform_f_T())
    # None column should be identically zero across every window.
    assert (res.posteriors[:, TYPE_INDEX["None"]] == 0.0).all()
    assert res.summed_posterior[TYPE_INDEX["None"]] == 0.0


def test_human_wins_over_none_by_construction():
    """With None excluded, a founder Human haploid must decode to Human as
    the top-1 candidate; the Human/None near-tie collapses to Human."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = _sample_haploid_from_type("Human", seed=42)
    res = decode_haplotype(haploid, E[0], _uniform_f_T())
    assert res.top_k_indices[0] == TYPE_INDEX["Human"]


def test_uniform_f_T_input_survives_none_zeroing():
    """A uniform f_T that includes None must be renormalized to uniform over
    active types, not rejected."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = _sample_haploid_from_type("Uber", seed=0)
    naive_uniform = np.ones(N_TYPES) / N_TYPES  # includes None
    res = decode_haplotype(haploid, E[0], naive_uniform)
    assert res.top_k_indices[0] == TYPE_INDEX["Uber"]
    assert TYPE_INDEX["None"] not in res.top_k_indices


def test_all_none_f_T_rejected():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = np.zeros(9, dtype=np.uint8)
    f_T = np.zeros(N_TYPES)
    f_T[TYPE_INDEX["None"]] = 1.0
    with pytest.raises(ValueError):
        decode_haplotype(haploid, E[0], f_T)


def test_large_t_posterior_approaches_prior():
    """At the largest grid t, emissions are near-uniform, so per-window
    posterior should be close to f_T."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    # Random haploid
    rng = np.random.default_rng(7)
    haploid = rng.integers(0, 256, size=9, dtype=np.uint8)
    f_T = _uniform_f_T()
    res = decode_haplotype(haploid, E[-1], f_T)
    # Per-window marginal should be close to uniform.
    for w in range(N_WINDOWS):
        assert np.allclose(res.posteriors[w], f_T, atol=0.05)


def test_posteriors_sum_to_one_per_window():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = _sample_haploid_from_type("Community", seed=3)
    res = decode_haplotype(haploid, E[0], _uniform_f_T())
    assert np.allclose(res.posteriors.sum(axis=1), 1.0)


def test_top_k_size():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = _sample_haploid_from_type("Uber", seed=0)
    res3 = decode_haplotype(haploid, E[0], _uniform_f_T(), top_k=3)
    res5 = decode_haplotype(haploid, E[0], _uniform_f_T(), top_k=5)
    assert res3.top_k_indices.shape == (3,)
    assert res5.top_k_indices.shape == (5,)
    # Top 3 of top_k=5 must equal top_k=3.
    assert np.array_equal(res5.top_k_indices[:3], res3.top_k_indices)


def test_bad_f_T_rejected():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    haploid = np.zeros(9, dtype=np.uint8)
    with pytest.raises(ValueError):
        decode_haplotype(haploid, E[0], np.array([1.0]))  # wrong shape
    with pytest.raises(ValueError):
        decode_haplotype(haploid, E[0], np.zeros(N_TYPES))  # doesn't sum to 1


def test_nonuniform_f_T_shifts_ambiguous_posterior():
    """If f_T heavily favors one type, an emission-ambiguous haploid should
    decode to that type."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    # Fully-uniform haploid values → emission likelihoods are close for many types.
    haploid = np.array([3, 128, 128, 200, 128, 128, 200, 200, 128], dtype=np.uint8)
    f_T = np.zeros(N_TYPES)
    f_T[TYPE_INDEX["Uber"]] = 0.99
    f_T += 0.01 / N_TYPES
    f_T /= f_T.sum()
    res = decode_haplotype(haploid, E[0], f_T)
    # With f_T that concentrated, Uber should still float to the top on
    # observations that don't strongly rule it out.
    # (This is a sanity test that f_T flows through the transition.)
    # For the top type at each window, we compare to the uniform case.
    res_uniform = decode_haplotype(haploid, E[0], _uniform_f_T())
    # The Uber-summed posterior mass should grow when we bias f_T toward Uber.
    assert res.summed_posterior[TYPE_INDEX["Uber"]] > res_uniform.summed_posterior[TYPE_INDEX["Uber"]]
