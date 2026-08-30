"""Analytic founder priors π_{T, ℓ, a} from §2 of synthetic-population-genetics.md.

Produces an 8 × 9 × 256 stochastic tensor whose (T, ℓ, ·) slice is the per-locus
allele distribution induced by `Haploid::from_type(T)` in `dc34-api/src/lib.rs`.

No simulation required — every entry is closed-form uniform-on-a-range with two
Goon/Uber delta-spike overrides.
"""

from __future__ import annotations

import numpy as np

N_LOCI = 9
N_ALLELES = 256

LOCUS = dict(
    cd_period=0, cd_rate=1, cd_dir=2, sat=3,
    hue_ratedir=4, hue_base=5, hue_bound=6,
    chaser=7, nonlin=8,
)

TYPES = ("Goon", "Community", "Village", "Human",
         "Other", "CtfContest", "Uber", "None")
N_TYPES = len(TYPES)
TYPE_INDEX = {t: i for i, t in enumerate(TYPES)}

# The `None` type is a theoretical category in `dc34-api`; zero real DC34
# badges were `None`-typed. The ancestry HMM excludes it from every stage so
# no user's composition ever shows non-zero None mass, and Human/None ties
# (their founder ranges differ only in cd_period_max) always resolve to Human.
INACTIVE_TYPES = ("None",)
ACTIVE_TYPE_MASK = np.array(
    [T not in INACTIVE_TYPES for T in TYPES], dtype=bool
)
ACTIVE_TYPES = tuple(T for T in TYPES if T not in INACTIVE_TYPES)
N_ACTIVE_TYPES = len(ACTIVE_TYPES)
ACTIVE_TYPE_INDICES = np.array(
    [i for i, T in enumerate(TYPES) if T not in INACTIVE_TYPES], dtype=np.int64
)

# Per-type ranges from §2 of the pop-gen doc.
# hue = hue_base range (drawn uniformly on this interval, then hue_bound is
# drawn uniformly from [hue_base, hue_range_end]).
_TYPE_TABLE = {
    "Goon":       dict(hue=(0, 20),    sat=(160, 255), chaser=(90, 255),
                        nonlin=(0, 255), cd_dir=(0, 255), cd_period_max=4),
    "Community":  dict(hue=(32, 80),   sat=(32, 160),  chaser=(90, 255),
                        nonlin=(0, 255), cd_dir=(0, 255), cd_period_max=2),
    "Village":    dict(hue=(80, 128),  sat=(32, 160),  chaser=(90, 255),
                        nonlin=(0, 255), cd_dir=(0, 45),  cd_period_max=4),
    "Human":      dict(hue=(128, 160), sat=(32, 255),  chaser=(90, 255),
                        nonlin=(0, 255), cd_dir=(0, 255), cd_period_max=5),
    "Other":      dict(hue=(160, 192), sat=(16, 255),  chaser=(0, 255),
                        nonlin=(0, 90),  cd_dir=(0, 255), cd_period_max=6),
    "CtfContest": dict(hue=(192, 220), sat=(16, 255),  chaser=(90, 255),
                        nonlin=(0, 90),  cd_dir=(0, 255), cd_period_max=6),
    "Uber":       dict(hue=(220, 255), sat=(130, 255), chaser=(0, 45),
                        nonlin=(0, 44),  cd_dir=(0, 45),  cd_period_max=3),
    "None":       dict(hue=(128, 160), sat=(32, 255),  chaser=(90, 255),
                        nonlin=(0, 255), cd_dir=(0, 255), cd_period_max=4),
}


def _uniform_range(lo: int, hi: int) -> np.ndarray:
    """PMF uniform on integer range [lo, hi] inclusive, over 0..=255."""
    out = np.zeros(N_ALLELES, dtype=np.float64)
    n = hi - lo + 1
    out[lo:hi + 1] = 1.0 / n
    return out


def _hue_base_pmf(hue_lo: int, hue_hi: int) -> np.ndarray:
    """P(hue_base = a) marginal after Haploid::from_type."""
    return _uniform_range(hue_lo, hue_hi)


def _hue_bound_pmf(hue_lo: int, hue_hi: int) -> np.ndarray:
    """P(hue_bound = a) after marginalizing over hue_base.

    Sampling: hue_base ~ Uniform[hue_lo, hue_hi]; hue_bound ~ Uniform[hue_base, hue_hi].
    So P(bound = b) = sum_{a <= b, a in [hue_lo, hue_hi]} P(base=a) * 1/(hue_hi - a + 1),
    for b in [hue_lo, hue_hi]. Outside that support the pmf is 0.
    """
    out = np.zeros(N_ALLELES, dtype=np.float64)
    n_base = hue_hi - hue_lo + 1
    for a in range(hue_lo, hue_hi + 1):
        span = hue_hi - a + 1
        for b in range(a, hue_hi + 1):
            out[b] += (1.0 / n_base) * (1.0 / span)
    return out


def build_founder_panel() -> np.ndarray:
    """Return the 8 × 9 × 256 founder prior tensor π_{T, ℓ, a}."""
    pi = np.zeros((N_TYPES, N_LOCI, N_ALLELES), dtype=np.float64)

    for T, row in _TYPE_TABLE.items():
        ti = TYPE_INDEX[T]

        # cd_period: uniform on 0..=cd_period_max
        pi[ti, LOCUS["cd_period"]] = _uniform_range(0, row["cd_period_max"])

        # cd_rate: unconstrained u8 for all types
        pi[ti, LOCUS["cd_rate"]] = _uniform_range(0, 255)

        # cd_dir: type-specific range
        pi[ti, LOCUS["cd_dir"]] = _uniform_range(*row["cd_dir"])

        # sat: type-specific range
        pi[ti, LOCUS["sat"]] = _uniform_range(*row["sat"])

        # hue_ratedir: unconstrained u8 for all types
        pi[ti, LOCUS["hue_ratedir"]] = _uniform_range(0, 255)

        # hue_base: uniform on the type's hue range, with Goon override.
        hue_lo, hue_hi = row["hue"]
        if T == "Goon":
            # `hue_base` is forced to 0 regardless of the range roll.
            spike = np.zeros(N_ALLELES, dtype=np.float64)
            spike[0] = 1.0
            pi[ti, LOCUS["hue_base"]] = spike
        else:
            pi[ti, LOCUS["hue_base"]] = _hue_base_pmf(hue_lo, hue_hi)

        # hue_bound: marginal over hue_base, with Uber override.
        if T == "Uber":
            # `hue_bound` is forced to 255 regardless of the range roll.
            spike = np.zeros(N_ALLELES, dtype=np.float64)
            spike[255] = 1.0
            pi[ti, LOCUS["hue_bound"]] = spike
        else:
            pi[ti, LOCUS["hue_bound"]] = _hue_bound_pmf(hue_lo, hue_hi)

        # chaser: type-specific range
        pi[ti, LOCUS["chaser"]] = _uniform_range(*row["chaser"])

        # nonlin: type-specific range
        pi[ti, LOCUS["nonlin"]] = _uniform_range(*row["nonlin"])

    return pi


FOUNDER_PANEL = build_founder_panel()
