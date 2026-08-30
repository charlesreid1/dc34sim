"""Mutation kernels K_rate for the Gray-code mutation model of §5, pop-gen doc.

For each rate level, K_rate[a, a'] = P(new_allele = a' | old_allele = a) under
one mutation pass. Each row sums to 1.

`cd_period` uses a specialized kernel that folds the shipped `% 7` post-processing
into the transition (§5.2, "cd_period wraps mod 7 post-mutation"). Every other
locus uses the general Gray-code kernel.

K^t for the fixed grid t ∈ {0, 1, 2, 4, 8, 16, 32, 64, 128} is precomputed via
repeated squaring per rate.
"""

from __future__ import annotations

import numpy as np

N_ALLELES = 256

# Rate table from §5.1.
RATE_VALUES = {"None": 0, "Baseline": 64, "Elevated": 100,
               "Radioactive": 140, "Apocalyptic": 240}
RATE_BITS = {"None": 0, "Baseline": 0x01, "Elevated": 0x03,
             "Radioactive": 0x07, "Apocalyptic": 0x1F}

RATES = ("None", "Baseline", "Elevated", "Radioactive", "Apocalyptic")

# Fixed grid of t values (§4.3).
T_GRID = (0, 1, 2, 4, 8, 16, 32, 64, 128)

# Gray-code LUTs.
_g = np.arange(N_ALLELES, dtype=np.uint8)
GRAY_ENCODE = (_g ^ (_g >> 1)).astype(np.uint8)
GRAY_DECODE = np.empty(N_ALLELES, dtype=np.uint8)
GRAY_DECODE[GRAY_ENCODE] = _g


def build_general_kernel(rate_name: str) -> np.ndarray:
    """Return a 256 × 256 stochastic K for one mutation pass at `rate_name`.

    Semantics (§5.1, §5.2): with per-locus probability p = value/256, apply
    `gray_decode(gray_encode(a) ^ (bits << shift))` where shift ~ Uniform{0..7}.
    Otherwise (prob 1 - p) keep `a`.
    """
    val = RATE_VALUES[rate_name]
    bits = RATE_BITS[rate_name]
    p = val / 256.0
    K = np.zeros((N_ALLELES, N_ALLELES), dtype=np.float64)

    # No-op branch.
    K += (1.0 - p) * np.eye(N_ALLELES, dtype=np.float64)

    if val == 0:
        return K

    # Fire branch: uniform over 8 shift choices.
    enc = GRAY_ENCODE.astype(np.int64)          # (256,)
    weight = p / 8.0
    for shift in range(8):
        mask = ((bits << shift) & 0xFF)
        flipped_enc = np.bitwise_xor(enc, mask).astype(np.uint8)
        new_alleles = GRAY_DECODE[flipped_enc]  # (256,) -- new allele for each old
        # Add weight to K[old, new] for every old.
        np.add.at(K, (np.arange(N_ALLELES), new_alleles), weight)

    return K


def build_cd_period_kernel(rate_name: str) -> np.ndarray:
    """Kernel for `cd_period`, where the shipped `mutated %= 7` is applied
    post-mutation on fire branches only. Support stays inside 0..255 but is
    concentrated on 0..=6; entries outside 0..=6 in the old-axis are still
    valid (mutation could in principle land there transiently), but the fire
    branch always folds mod 7.

    We build the full 256×256 kernel so it composes with itself under matmul
    the same way the general kernel does.
    """
    val = RATE_VALUES[rate_name]
    bits = RATE_BITS[rate_name]
    p = val / 256.0
    K = np.zeros((N_ALLELES, N_ALLELES), dtype=np.float64)

    K += (1.0 - p) * np.eye(N_ALLELES, dtype=np.float64)

    if val == 0:
        return K

    enc = GRAY_ENCODE.astype(np.int64)
    weight = p / 8.0
    for shift in range(8):
        mask = ((bits << shift) & 0xFF)
        flipped_enc = np.bitwise_xor(enc, mask).astype(np.uint8)
        new_alleles = GRAY_DECODE[flipped_enc] % 7   # <-- the mod 7 fold
        np.add.at(K, (np.arange(N_ALLELES), new_alleles), weight)

    return K


def matrix_power_grid(K: np.ndarray, grid=T_GRID) -> dict:
    """Return {t: K^t} for every t in `grid` via repeated squaring.

    Grid values need not be powers of two; we compute K^t individually with
    numpy's optimized `matrix_power` which uses binary exponentiation.
    """
    return {t: np.linalg.matrix_power(K, t) if t > 0 else np.eye(N_ALLELES)
            for t in grid}


def build_emission_tensor(pi_TLA: np.ndarray, rate_name: str,
                          grid=T_GRID) -> np.ndarray:
    """Return E[t_index, T, ℓ, a] = (π_{T,ℓ} · K^t)(a) for every grid t.

    Uses the cd_period-specific kernel for locus 0 and the general kernel for
    every other locus. Shape: (len(grid), n_types, n_loci, n_alleles).
    """
    from hmm.founder_panel import LOCUS
    cd_period_locus = LOCUS["cd_period"]

    K_gen = build_general_kernel(rate_name)
    K_cd = build_cd_period_kernel(rate_name)

    Kpow_gen = matrix_power_grid(K_gen, grid)
    Kpow_cd = matrix_power_grid(K_cd, grid)

    n_types, n_loci, n_alleles = pi_TLA.shape
    E = np.empty((len(grid), n_types, n_loci, n_alleles), dtype=np.float64)
    for ti, t in enumerate(grid):
        for ell in range(n_loci):
            K = Kpow_cd[t] if ell == cd_period_locus else Kpow_gen[t]
            E[ti, :, ell, :] = pi_TLA[:, ell, :] @ K

    return E


def snap_t_to_grid(t: float, grid=T_GRID) -> int:
    """Snap a continuous t to the nearest grid index (§4.3)."""
    grid_arr = np.asarray(grid)
    return int(np.argmin(np.abs(grid_arr - t)))
