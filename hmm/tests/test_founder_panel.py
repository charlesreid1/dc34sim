import numpy as np

from hmm.founder_panel import (
    FOUNDER_PANEL, LOCUS, N_ALLELES, N_LOCI, N_TYPES, TYPE_INDEX,
    build_founder_panel,
)


def test_shape():
    assert FOUNDER_PANEL.shape == (N_TYPES, N_LOCI, N_ALLELES)


def test_normalization():
    sums = FOUNDER_PANEL.sum(axis=-1)
    assert np.allclose(sums, 1.0), "every (T, ℓ) marginal must sum to 1"


def test_uniform_loci_are_type_invariant():
    # cd_rate and hue_ratedir are unconstrained u8 for every type per §2.
    for locus_name in ("cd_rate", "hue_ratedir"):
        ell = LOCUS[locus_name]
        slice_ = FOUNDER_PANEL[:, ell, :]
        assert np.allclose(slice_, slice_[0]), f"{locus_name} must be type-invariant"
        assert np.allclose(slice_[0], 1.0 / N_ALLELES)


def test_goon_hue_base_is_delta_zero():
    ti = TYPE_INDEX["Goon"]
    pmf = FOUNDER_PANEL[ti, LOCUS["hue_base"]]
    assert pmf[0] == 1.0
    assert pmf[1:].sum() == 0.0


def test_uber_hue_bound_is_delta_255():
    ti = TYPE_INDEX["Uber"]
    pmf = FOUNDER_PANEL[ti, LOCUS["hue_bound"]]
    assert pmf[255] == 1.0
    assert pmf[:255].sum() == 0.0


def test_uber_chaser_supported_on_0_45():
    ti = TYPE_INDEX["Uber"]
    pmf = FOUNDER_PANEL[ti, LOCUS["chaser"]]
    assert np.isclose(pmf[0:46].sum(), 1.0)
    assert pmf[46:].sum() == 0.0


def test_human_hue_base_supported_on_128_160():
    ti = TYPE_INDEX["Human"]
    pmf = FOUNDER_PANEL[ti, LOCUS["hue_base"]]
    assert np.isclose(pmf[128:161].sum(), 1.0)
    assert pmf[:128].sum() == 0.0 and pmf[161:].sum() == 0.0


def test_other_nonlin_supported_on_0_90():
    ti = TYPE_INDEX["Other"]
    pmf = FOUNDER_PANEL[ti, LOCUS["nonlin"]]
    assert np.isclose(pmf[0:91].sum(), 1.0)
    assert pmf[91:].sum() == 0.0


def test_hue_bound_geq_hue_base_in_support():
    # For non-override types, hue_bound support is contained in [hue_lo, hue_hi],
    # same as hue_base's support, so every mass point of bound has some base <= it.
    for T, ti in TYPE_INDEX.items():
        if T in ("Goon", "Uber"):
            continue
        base = FOUNDER_PANEL[ti, LOCUS["hue_base"]]
        bound = FOUNDER_PANEL[ti, LOCUS["hue_bound"]]
        base_support = np.nonzero(base)[0]
        bound_support = np.nonzero(bound)[0]
        assert bound_support.min() >= base_support.min()
        assert bound_support.max() <= base_support.max()


def test_builder_is_pure():
    # Rebuilding gives an identical tensor.
    a = build_founder_panel()
    b = build_founder_panel()
    assert np.array_equal(a, b)
