"""Period input-output calculations."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .model import DEFAULT_INITIAL_CONDITION, rhs


class PeriodComputationError(RuntimeError):
    """Raised when the period cannot be computed from detected maxima."""


@dataclass(frozen=True)
class PeriodSettings:
    """Settings for period calculations."""

    initial_condition: tuple[float, float, float] = DEFAULT_INITIAL_CONDITION
    t_start: float = 0.0
    t_end: float = 80.0
    transient_time: float = 50.0
    method: str = "LSODA"
    rtol: float = 1.0e-9
    atol: float = 1.0e-11
    max_step: float = 0.01
    periods_to_average: int = 1


def compute_period(K_d: float, A: float, settings: PeriodSettings | None = None) -> float:
    """Compute the period by measuring times between successive maxima of `M(t)`."""

    settings = settings or PeriodSettings()
    maxima = detect_maxima_times(K_d, A, settings)

    if len(maxima) < settings.periods_to_average + 1:
        raise PeriodComputationError(
            f"Need {settings.periods_to_average + 1} maxima after transient; found {len(maxima)} "
            f"for K_d={K_d}, A={A}."
        )

    differences = np.diff(maxima)
    if settings.periods_to_average == 1:
        return float(differences[0])

    return float(np.mean(differences[-settings.periods_to_average :]))


def detect_maxima_times(K_d: float, A: float, settings: PeriodSettings | None = None) -> np.ndarray:
    """Integrate the model and return times where `M(t)` has a local maximum."""

    settings = settings or PeriodSettings()

    def dMdt_event(t: float, state: np.ndarray) -> float:
        return rhs(t, state, A, K_d)[0]

    dMdt_event.direction = -1

    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, A, K_d),
        t_span=(settings.t_start, settings.t_end),
        y0=settings.initial_condition,
        method=settings.method,
        rtol=settings.rtol,
        atol=settings.atol,
        max_step=settings.max_step,
        events=dMdt_event,
    )

    if not sol.success:
        raise PeriodComputationError(f"Integration failed for K_d={K_d}, A={A}: {sol.message}")

    maxima = sol.t_events[0]
    return maxima[maxima >= settings.transient_time]


def compute_period_curve(
    K_d: float,
    A_grid: np.ndarray,
    settings: PeriodSettings | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute `P(K_d, A)` over a grid of input values."""

    settings = settings or PeriodSettings()
    periods = np.array([compute_period(K_d, float(A), settings) for A in A_grid], dtype=float)
    return np.asarray(A_grid, dtype=float), periods


def finite_difference_period_derivative(A_grid: np.ndarray, periods: np.ndarray) -> np.ndarray:
    """Estimate dP/dA from a computed period curve."""

    A_grid = np.asarray(A_grid, dtype=float)
    periods = np.asarray(periods, dtype=float)
    edge_order = 2 if len(A_grid) >= 3 else 1
    return np.gradient(periods, A_grid, edge_order=edge_order)
