import numpy as np
import pytest

from hmm.compose import DEFAULT_MIN_THRESHOLD, compose_percentages


def test_sums_to_100():
    p = np.array([0.4, 0.35, 0.15, 0.1])
    out = compose_percentages(p)
    assert out.sum() == 100


def test_below_threshold_dropped():
    p = np.array([0.7, 0.25, 0.03, 0.02])
    out = compose_percentages(p, min_threshold=0.05)
    assert out[2] == 0
    assert out[3] == 0
    assert out.sum() == 100


def test_all_below_threshold_keeps_argmax():
    p = np.array([0.03, 0.04, 0.02, 0.01])  # everything sub-5% after renorm
    # After normalization, largest is 0.04/0.10 = 0.4, so it's above 0.05 —
    # try a case where every entry is truly small.
    p2 = np.array([0.02, 0.02, 0.02, 0.02])
    out = compose_percentages(p2, min_threshold=0.9)
    assert out.sum() == 100
    # Exactly one nonzero.
    assert (out > 0).sum() == 1


def test_uniform_composition():
    p = np.array([0.25, 0.25, 0.25, 0.25])
    out = compose_percentages(p, min_threshold=0.05)
    assert out.sum() == 100
    # Should be 25/25/25/25.
    assert out.tolist() == [25, 25, 25, 25]


def test_thirds_round_correctly():
    p = np.array([1 / 3, 1 / 3, 1 / 3])
    out = compose_percentages(p, min_threshold=0.05)
    assert out.sum() == 100
    assert sorted(out.tolist()) == [33, 33, 34]


def test_default_threshold():
    p = np.array([0.94, 0.04, 0.02])
    out = compose_percentages(p)  # default 0.05
    assert out.sum() == 100
    assert out[1] == 0 and out[2] == 0
    assert out[0] == 100


def test_bad_input_rejected():
    with pytest.raises(ValueError):
        compose_percentages(np.array([[0.5, 0.5]]))
    with pytest.raises(ValueError):
        compose_percentages(np.array([-0.1, 1.1]))
    with pytest.raises(ValueError):
        compose_percentages(np.zeros(3))


def test_unnormalized_input_normalized_first():
    a = compose_percentages(np.array([2.0, 3.0, 5.0]))
    b = compose_percentages(np.array([0.2, 0.3, 0.5]))
    assert a.tolist() == b.tolist()


def test_deterministic_tiebreak():
    """Two equal-remainder entries should give the extra point to the first
    (stable sort) — the output must be identical across repeated calls."""
    p = np.array([0.335, 0.335, 0.33])
    out_a = compose_percentages(p)
    out_b = compose_percentages(p)
    assert out_a.tolist() == out_b.tolist()
