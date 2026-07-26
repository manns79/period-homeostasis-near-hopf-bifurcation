"""Utilities for locating infinitesimal period homeostasis points."""

from dataclasses import dataclass

import numpy as np

from .hopf import periodic_interval
from .period import PeriodSettings, compute_period_curve, finite_difference_period_derivative


@dataclass(frozen=True)
class IphPoint:
    """A detected critical point of the period input-output function."""

    K_d: float
    A: float
    period: float
    derivative_left: float
    derivative_right: float


def find_iph_points_from_curve(K_d: float, A_grid: np.ndarray, periods: np.ndarray) -> list[IphPoint]:
    """Find approximate IPH points from a sampled period curve."""

    derivatives = finite_difference_period_derivative(A_grid, periods)
    points: list[IphPoint] = []

    for i in range(len(A_grid) - 1):
        left_derivative = derivatives[i]
        right_derivative = derivatives[i + 1]

        if not (np.isfinite(left_derivative) and np.isfinite(right_derivative)):
            continue

        if left_derivative == 0:
            A_iph = float(A_grid[i])
            period = float(periods[i])
        elif left_derivative * right_derivative < 0:
            weight = abs(left_derivative) / (abs(left_derivative) + abs(right_derivative))
            A_iph = float(A_grid[i] + weight * (A_grid[i + 1] - A_grid[i]))
            period = float(np.interp(A_iph, A_grid, periods))
        else:
            continue

        points.append(
            IphPoint(
                K_d=K_d,
                A=A_iph,
                period=period,
                derivative_left=float(left_derivative),
                derivative_right=float(right_derivative),
            )
        )

    return points


def find_iph_points(
    K_d: float,
    *,
    n_A: int = 200,
    settings: PeriodSettings | None = None,
    A_bounds: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[IphPoint]]:
    """Compute a period curve for one `K_d` slice and locate IPH points."""

    if A_bounds is None:
        A_bounds = periodic_interval(K_d)

    A_grid = np.linspace(A_bounds[0], A_bounds[1], n_A)
    A_grid, periods = compute_period_curve(K_d, A_grid, settings)
    points = find_iph_points_from_curve(K_d, A_grid, periods)
    return A_grid, periods, points
