import numpy as np
import pytest

from hmm.founder_panel import FOUNDER_PANEL, TYPE_INDEX
from hmm.mutation_kernel import build_emission_tensor
from hmm.popsim_bridge import (
    PopsimSnapshot, biased_snapshot, decode_composition, uniform_snapshot,
)


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


def test_uniform_snapshot():
    """Uniform over the 7 ACTIVE types; None (index 7) is excluded."""
    snap = uniform_snapshot(t=0)
    assert snap.t == 0
    assert np.isclose(snap.f_T.sum(), 1.0)
    # 7 active types, each 1/7.
    assert snap.f_T[TYPE_INDEX["None"]] == 0.0
    active_vals = np.delete(snap.f_T, TYPE_INDEX["None"])
    assert np.allclose(active_vals, 1 / 7)


def test_biased_snapshot():
    snap = biased_snapshot({"Human": 0.7, "Uber": 0.01, "Goon": 0.29}, t=5)
    assert snap.t == 5
    assert np.isclose(snap.f_T.sum(), 1.0)
    assert np.isclose(snap.f_T[TYPE_INDEX["Human"]], 0.7 / 1.0)


def test_biased_snapshot_missing_type_gets_zero():
    snap = biased_snapshot({"Human": 1.0})
    assert snap.f_T[TYPE_INDEX["Uber"]] == 0.0
    assert snap.f_T[TYPE_INDEX["Human"]] == 1.0


def test_biased_snapshot_rejects_zero_sum():
    with pytest.raises(ValueError):
        biased_snapshot({})


def test_biased_snapshot_silently_drops_none():
    """The ancestry decoder never models None; requesting it in a bias dict
    is a no-op rather than an error, so callers reading population fractions
    from popsim don't need to filter."""
    snap = biased_snapshot({"Human": 0.7, "None": 0.3})
    assert snap.f_T[TYPE_INDEX["None"]] == 0.0
    assert snap.f_T[TYPE_INDEX["Human"]] == 1.0  # 0.7 / 0.7 after None dropped


def test_biased_snapshot_only_none_rejected():
    with pytest.raises(ValueError):
        biased_snapshot({"None": 1.0})


def test_decode_composition_end_to_end():
    dip = _make_diploid("Uber", "Human", seed=1)
    snap = uniform_snapshot(t=0)
    rep = decode_composition(dip, snap)
    assert rep.percentages.sum() == 100
    assert np.isclose(rep.composition.sum(), 1.0)
    assert rep.C0.shape == (3,)
    assert rep.C1.shape == (3,)
    assert rep.joint_posteriors.shape == (4, 3, 3, 2)
    assert rep.stage1_h0_posteriors.shape == (4, 8)
    assert rep.stage1_h1_posteriors.shape == (4, 8)


def test_pure_uber_diploid_reports_uber_dominant():
    dip = _make_diploid("Uber", "Uber", seed=2)
    snap = uniform_snapshot(t=0)
    rep = decode_composition(dip, snap)
    assert rep.percentages[TYPE_INDEX["Uber"]] >= 50
    assert rep.percentages[TYPE_INDEX["Uber"]] == rep.percentages.max()


def test_snapshot_t_is_snapped_to_grid():
    dip = _make_diploid("Uber", "Human", seed=3)
    # Ask for t=3; nearest grid points are 2 and 4, so t_used is one of them.
    snap = PopsimSnapshot(t=3, f_T=np.ones(8) / 8)
    rep = decode_composition(dip, snap)
    assert rep.t_used in (2, 4)


def test_cached_emission_tensor_gives_same_answer():
    E = build_emission_tensor(FOUNDER_PANEL, "Baseline")
    dip = _make_diploid("Goon", "Community", seed=5)
    snap = uniform_snapshot(t=0)
    rep_fresh = decode_composition(dip, snap)
    rep_cached = decode_composition(dip, snap, emission_tensor=E)
    assert np.array_equal(rep_fresh.percentages, rep_cached.percentages)
    assert np.allclose(rep_fresh.composition, rep_cached.composition)


def test_biased_f_T_shifts_percentages():
    """A strong f_T bias toward Uber should raise the reported Uber percentage
    on an ambiguous diploid."""
    dip = _make_diploid("Uber", "Human", seed=6)
    uniform = uniform_snapshot(t=32)  # ambient t so emissions are smeared
    uber_biased = biased_snapshot({"Uber": 0.9, "Human": 0.1}, t=32)
    rep_uniform = decode_composition(dip, uniform)
    rep_biased = decode_composition(dip, uber_biased)
    assert (rep_biased.composition[TYPE_INDEX["Uber"]]
            > rep_uniform.composition[TYPE_INDEX["Uber"]])
