"""Hopf-curve utilities for the Kim-Forger SNF model."""

import numpy as np
from scipy.optimize import brentq

A_GENERALIZED_HOPF = 13 / 1616
A_SINGULAR = 1 / 8


def hopf_curve(A: float | np.ndarray) -> float | np.ndarray:
    """Return `K_d` on the Hopf curve as a function of `A`.

    This is the formula presented in the manuscript.
    """

    A = np.asarray(A)
    return (128 * A**2 + 240 * A + 49 - (7 + 16 * A) * np.sqrt(256 * A + 49)) / (16 * (1 - 8 * A))


def find_hopf_points(
    K_d: float,
    *,
    A_min: float = 1.0e-8,
    A_max: float = A_SINGULAR - 1.0e-8,
    n_brackets: int = 4000,
    require_two: bool = True,
) -> tuple[float, ...]:
    """Find values of `A` where the Hopf curve has the requested `K_d`."""

    A_grid = np.linspace(A_min, A_max, n_brackets + 1)
    values = hopf_curve(A_grid) - K_d
    roots: list[float] = []

    for left, right, f_left, f_right in zip(A_grid[:-1], A_grid[1:], values[:-1], values[1:]):
        if not (np.isfinite(f_left) and np.isfinite(f_right)):
            continue

        if f_left == 0:
            roots.append(float(left))
        elif f_left * f_right < 0:
            roots.append(float(brentq(lambda A: hopf_curve(A) - K_d, left, right)))

    if values[-1] == 0:
        roots.append(float(A_grid[-1]))

    unique_roots = _deduplicate_roots(roots)
    if require_two and len(unique_roots) != 2:
        raise ValueError(f"Expected two Hopf points for K_d={K_d}, found {len(unique_roots)}.")

    return tuple(unique_roots)


def periodic_interval(
    K_d: float,
    *,
    margin_fraction: float = 0.005,
    margin_absolute: float = 1.0e-5,
) -> tuple[float, float]:
    """Return an inward-shifted `A` interval inside the periodic regime."""

    left, right = find_hopf_points(K_d)
    width = right - left
    margin = max(margin_absolute, margin_fraction * width)
    if 2 * margin >= width:
        raise ValueError(f"Requested margin is too large for K_d={K_d}.")
    return left + margin, right - margin


def _deduplicate_roots(roots: list[float], *, tol: float = 1.0e-9) -> list[float]:
    """Return sorted roots, removing adjacent values within `tol`."""
    roots = sorted(roots)
    unique: list[float] = []
    for root in roots:
        if not unique or abs(root - unique[-1]) > tol:
            unique.append(root)
    return unique
