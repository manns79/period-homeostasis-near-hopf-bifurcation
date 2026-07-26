"""Dimensionless Kim-Forger SNF model."""

import numpy as np


DEFAULT_INITIAL_CONDITION = (0.01, 0.01, 0.006)


def free_bmal_fraction(P: float | np.ndarray, A: float, K_d: float) -> float | np.ndarray:
    """Return the free BMAL1:CLOCK fraction used in the SNF model."""

    return (1 - P / A - K_d / A + np.sqrt((1 - P / A - K_d / A) ** 2 + 4 * K_d / A)) / 2


def rhs(t: float, state: tuple[float, float, float] | np.ndarray, A: float, K_d: float) -> list[float]:
    """Right-hand side of the dimensionless Kim-Forger SNF model with the argument order matching `scipy.integrate.solve_ivp`.
    """

    M, P_c, P = state
    f = free_bmal_fraction(P, A, K_d)
    return [f - M, M - P_c, P_c - P]


