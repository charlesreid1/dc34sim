import numpy as np
import pytest

from hmm.founder_panel import FOUNDER_PANEL, N_TYPES, TYPE_INDEX
from hmm.full_stage import N_WINDOWS, decode_haplotype
from hmm.mutation_kernel import build_emission_tensor
from hmm.simplified_stage import (
    W_DE_INDEX, decode_diploid, phenotype_nonlin,
)


def _uniform_f_T():
    return np.ones(N_TYPES) / N_TYPES


def _sample_haploid_from_type(T: str, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ti = TYPE_INDEX[T]
    pi = FOUNDER_PANEL[ti]
    out = np.zeros(9, dtype=np.uint8)
    for ell in range(9):
        out[ell] = rng.choice(256, p=pi[ell])
    return out


def _make_diploid(T0: str, T1: str, seed: int = 0) -> np.ndarray:
    dip = np.zeros((2, 9), dtype=np.uint8)
    dip[0] = _sample_haploid_from_type(T0, seed=seed)
    dip[1] = _sample_haploid_from_type(T1, seed=seed + 1)
    return dip


def _stage1_candidates(diploid, E_TLA, f_T, k=3):
    r0 = decode_haplotype(diploid[0], E_TLA, f_T, top_k=k)
    r1 = decode_haplotype(diploid[1], E_TLA, f_T, top_k=k)
    return r0.top_k_indices, r1.top_k_indices


def test_decode_diploid_shapes():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Uber", "Human", seed=1)
    C0, C1 = _stage1_candidates(dip, E[0], _uniform_f_T())
    res = decode_diploid(dip, E[0], _uniform_f_T(), C0, C1)
    assert res.joint_posteriors.shape == (N_WINDOWS, len(C0), len(C1), 2)
    assert res.composition.shape == (N_TYPES,)


def test_joint_posteriors_sum_to_one_per_window():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Uber", "Goon", seed=2)
    C0, C1 = _stage1_candidates(dip, E[0], _uniform_f_T())
    res = decode_diploid(dip, E[0], _uniform_f_T(), C0, C1)
    per_window = res.joint_posteriors.sum(axis=(1, 2, 3))
    assert np.allclose(per_window, 1.0)


def test_composition_sums_to_one():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Uber", "Human", seed=3)
    C0, C1 = _stage1_candidates(dip, E[0], _uniform_f_T())
    res = decode_diploid(dip, E[0], _uniform_f_T(), C0, C1)
    assert np.isclose(res.composition.sum(), 1.0)


def test_pure_uber_diploid_composition_favors_uber():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Uber", "Uber", seed=4)
    C0, C1 = _stage1_candidates(dip, E[0], _uniform_f_T())
    res = decode_diploid(dip, E[0], _uniform_f_T(), C0, C1)
    assert res.composition.argmax() == TYPE_INDEX["Uber"]


def test_uber_goon_composition_carries_both_types():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Uber", "Goon", seed=5)
    C0, C1 = _stage1_candidates(dip, E[0], _uniform_f_T())
    res = decode_diploid(dip, E[0], _uniform_f_T(), C0, C1)
    # Both Uber and Goon should carry meaningful mass in the composition.
    assert res.composition[TYPE_INDEX["Uber"]] > 0.1
    assert res.composition[TYPE_INDEX["Goon"]] > 0.1


def test_haplotype_reversal_flips_per_haplotype_posteriors_with_bug():
    """The nonlin bug's asymmetry lives in the per-haplotype (not aggregate)
    label posteriors on W_DE. The aggregate composition sums over haplotypes,
    so it is invariant under a full (haplo swap, C swap, switch swap) — that
    is a genuine symmetry of the model. The per-haplotype posteriors are not
    invariant, and this test pins that."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    # Find an Uber/Human pair whose phenotype.nonlin actually differs under swap.
    dip = None
    for seed in range(30):
        cand = _make_diploid("Uber", "Human", seed=seed)
        if phenotype_nonlin(cand) != phenotype_nonlin(cand[::-1]):
            dip = cand
            break
    assert dip is not None
    reversed_dip = dip[::-1].copy()
    r0 = decode_haplotype(dip[0], E[0], _uniform_f_T(), top_k=3)
    r1 = decode_haplotype(dip[1], E[0], _uniform_f_T(), top_k=3)
    union = np.array(sorted(set(list(r0.top_k_indices) + list(r1.top_k_indices))))
    res = decode_diploid(dip, E[0], _uniform_f_T(), union, union, use_nonlin_bug=True)
    res_rev = decode_diploid(
        reversed_dip, E[0], _uniform_f_T(), union, union, use_nonlin_bug=True,
    )
    # Per-haplotype 0 posterior on W_DE under original vs reversed diploid.
    p_h0_orig = res.joint_posteriors[W_DE_INDEX].sum(axis=(1, 2))  # (k,)
    p_h0_rev = res_rev.joint_posteriors[W_DE_INDEX].sum(axis=(1, 2))
    assert not np.allclose(p_h0_orig, p_h0_rev), (
        "haplo0 label posterior on W_DE was invariant under reversal — "
        "expected the nonlin bug to break that"
    )


def test_without_bug_composition_is_symmetric_under_reversal():
    """When the δ factor is dropped and we use the same candidate sets on both
    orderings, reversing haplo0 ↔ haplo1 leaves the composition invariant."""
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Uber", "Human", seed=6)
    reversed_dip = dip[::-1].copy()
    # Force identical candidate sets on both orderings by unioning both stage-1
    # calls, so the only difference in the two decodes is the observation order.
    r0 = decode_haplotype(dip[0], E[0], _uniform_f_T(), top_k=3)
    r1 = decode_haplotype(dip[1], E[0], _uniform_f_T(), top_k=3)
    union = np.array(sorted(set(list(r0.top_k_indices) + list(r1.top_k_indices))))
    res = decode_diploid(dip, E[0], _uniform_f_T(), union, union, use_nonlin_bug=False)
    res_rev = decode_diploid(
        reversed_dip, E[0], _uniform_f_T(), union, union, use_nonlin_bug=False,
    )
    # In the sum-over-haplotypes composition, the two orderings must agree
    # because we've made the emission model symmetric under slot exchange.
    assert np.allclose(res.composition, res_rev.composition, atol=1e-10)


def test_marginalization_over_switch_matches_per_haplotype():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Uber", "Community", seed=8)
    C0, C1 = _stage1_candidates(dip, E[0], _uniform_f_T())
    res = decode_diploid(dip, E[0], _uniform_f_T(), C0, C1)
    # Marginalizing joint over switch and over C1 gives per-window per-haplo0
    # posteriors; those must sum to 1 across the C0 axis.
    p_h0 = res.joint_posteriors.sum(axis=(2, 3))  # (n_windows, k0)
    assert np.allclose(p_h0.sum(axis=1), 1.0)


def test_bad_diploid_shape_rejected():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    with pytest.raises(ValueError):
        decode_diploid(
            np.zeros(9, dtype=np.uint8), E[0], _uniform_f_T(),
            np.array([0]), np.array([0]),
        )


def test_W_DE_switch1_zeroed_when_swap_would_change_phenotype():
    """Direct check: for a diploid where slot swap changes phenotype.nonlin,
    the δ factor zeroes W_DE's switch=1 emission slice."""
    from hmm.simplified_stage import _stage2_emission
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    # Use realistic founder-supported values so the base emission is nonzero.
    # slot0 = Uber-drawn, slot1 = Human-drawn (both in-support).
    dip = np.zeros((2, 9), dtype=np.uint8)
    dip[0] = _sample_haploid_from_type("Uber", seed=10)
    dip[1] = _sample_haploid_from_type("Human", seed=11)
    # Force a phenotype asymmetry: make slot0.chaser != slot1.chaser.
    dip[0, 7] = 10
    dip[1, 7] = 40
    dip[0, 8] = 50
    dip[1, 8] = 100
    assert phenotype_nonlin(dip) != phenotype_nonlin(dip[::-1])
    # Candidates include Uber and Human so emissions have support.
    C = np.array([TYPE_INDEX["Uber"], TYPE_INDEX["Human"], TYPE_INDEX["Community"]])
    emissions = _stage2_emission(dip, E[-1], C, C, use_nonlin_bug=True)
    assert emissions[W_DE_INDEX, :, :, 1].sum() == 0.0
    assert emissions[W_DE_INDEX, :, :, 0].sum() > 0.0
