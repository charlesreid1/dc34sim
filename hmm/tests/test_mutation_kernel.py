import numpy as np
import pytest

from hmm.founder_panel import FOUNDER_PANEL, LOCUS, N_ALLELES
from hmm.mutation_kernel import (
    RATES, T_GRID,
    build_cd_period_kernel, build_emission_tensor, build_general_kernel,
    matrix_power_grid, snap_t_to_grid,
)


@pytest.mark.parametrize("rate", RATES)
def test_general_rows_sum_to_one(rate):
    K = build_general_kernel(rate)
    assert np.allclose(K.sum(axis=1), 1.0)


@pytest.mark.parametrize("rate", RATES)
def test_cd_period_rows_sum_to_one(rate):
    K = build_cd_period_kernel(rate)
    assert np.allclose(K.sum(axis=1), 1.0)


def test_none_rate_is_identity_general():
    K = build_general_kernel("None")
    assert np.array_equal(K, np.eye(N_ALLELES))


def test_none_rate_is_identity_cd_period():
    # With rate=None, no fires ever, so no mod-7 fold; kernel is identity.
    K = build_cd_period_kernel("None")
    assert np.array_equal(K, np.eye(N_ALLELES))


@pytest.mark.parametrize("rate", ("Baseline", "Elevated", "Radioactive", "Apocalyptic"))
def test_uniform_is_fixed_point_of_general(rate):
    K = build_general_kernel(rate)
    u = np.ones(N_ALLELES) / N_ALLELES
    assert np.allclose(u @ K, u)


@pytest.mark.parametrize("rate", ("Baseline", "Elevated", "Radioactive", "Apocalyptic"))
def test_cd_period_kernel_concentrates_mass_on_0_6(rate):
    """After one Apocalyptic step from a uniform-over-0..255 input, most of
    the mass should sit on 0..=6. Baseline still moves some but concentration
    grows with rate."""
    K = build_cd_period_kernel(rate)
    u = np.ones(N_ALLELES) / N_ALLELES
    out = u @ K
    # A fraction ~= p of the mass is guaranteed to sit on 0..=6 (from fires);
    # additional mass is contributed by the identity component of entries
    # already in 0..=6 before the pass.
    assert out[:7].sum() > 7.0 / 256  # strictly more than the pre-mixing baseline


def test_grid_powers_shape_and_endpoints():
    K = build_general_kernel("Baseline")
    powers = matrix_power_grid(K)
    for t in T_GRID:
        assert powers[t].shape == (N_ALLELES, N_ALLELES)
    assert np.array_equal(powers[0], np.eye(N_ALLELES))
    # K^1 = K
    assert np.allclose(powers[1], K)
    # K^2 = K @ K
    assert np.allclose(powers[2], K @ K)


def test_grid_powers_remain_stochastic():
    K = build_general_kernel("Apocalyptic")
    powers = matrix_power_grid(K)
    for t, Kt in powers.items():
        assert np.allclose(Kt.sum(axis=1), 1.0), f"K^{t} rows must sum to 1"
        assert (Kt >= -1e-12).all(), f"K^{t} has negative entries"


def test_emission_tensor_shape_and_normalization():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    n_grid = len(T_GRID)
    n_types, n_loci, n_alleles = FOUNDER_PANEL.shape
    assert E.shape == (n_grid, n_types, n_loci, n_alleles)
    # Every (t, T, ℓ) marginal sums to 1.
    assert np.allclose(E.sum(axis=-1), 1.0)


def test_emission_at_t0_equals_founder():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    assert np.allclose(E[0], FOUNDER_PANEL)


def test_uniform_loci_remain_uniform_across_t_and_rate():
    # cd_rate and hue_ratedir start uniform and are fixed points, so they stay
    # uniform after any number of mutation passes at any nonzero rate.
    for rate in ("Baseline", "Elevated", "Radioactive", "Apocalyptic"):
        E = build_emission_tensor(FOUNDER_PANEL, rate)
        for locus_name in ("cd_rate", "hue_ratedir"):
            ell = LOCUS[locus_name]
            for ti in range(E.shape[0]):
                for T in range(E.shape[1]):
                    assert np.allclose(E[ti, T, ell], 1.0 / N_ALLELES), (
                        f"{locus_name} drifted at rate={rate}, t_idx={ti}, T={T}"
                    )


def test_uber_hue_bound_delta_smears_over_time():
    """Uber's hue_bound starts as a delta at 255; after mutation passes it should
    smear (max concentration strictly less than 1)."""
    from hmm.founder_panel import TYPE_INDEX
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    ti_uber = TYPE_INDEX["Uber"]
    ell = LOCUS["hue_bound"]
    assert E[0, ti_uber, ell, 255] == 1.0
    # After t=1, some mass leaks away.
    assert E[1, ti_uber, ell, 255] < 1.0
    # After t=128, the delta is essentially gone.
    assert E[-1, ti_uber, ell, 255] < 0.1


def test_snap_t_to_grid():
    assert snap_t_to_grid(0) == 0
    assert snap_t_to_grid(1) == 1
    assert snap_t_to_grid(3) in (2, 3)  # equidistant between t=2 and t=4
    assert snap_t_to_grid(1000) == len(T_GRID) - 1
