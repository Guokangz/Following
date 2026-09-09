"""Genuine order-by-order perturbative hierarchy for a driven 2LS cavity.

The hierarchy is obtained by inserting

``a=sum eta^n a_n``, ``s=sum eta^n s_n``, and ``p=sum eta^n p_n``

into the nonlinear mean-field equations.  It propagates the coefficients
directly and does not infer them by fitting multiple full solutions.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from .mean_field import MeanFieldTrajectory, TwoLevelMeanFieldParameters


ComplexEnvelope = Callable[[float], complex]


@dataclass(frozen=True)
class PerturbativeTrajectory:
    time: NDArray[np.float64]
    cavity_orders: NDArray[np.complex128]
    coherence_orders: NDArray[np.complex128]
    excited_population_orders: NDArray[np.float64]

    @property
    def max_order(self) -> int:
        return self.cavity_orders.shape[0] - 1

    def reconstruct(self, eta: float, *, through_order: int | None = None) -> MeanFieldTrajectory:
        """Reconstruct the physical trajectory through a selected order."""

        selected_order = self.max_order if through_order is None else through_order
        if selected_order < 0 or selected_order > self.max_order:
            raise ValueError("through_order lies outside the propagated hierarchy")
        powers = eta ** np.arange(selected_order + 1)
        cavity = np.tensordot(powers, self.cavity_orders[: selected_order + 1], axes=(0, 0))
        coherence = np.tensordot(
            powers,
            self.coherence_orders[: selected_order + 1],
            axes=(0, 0),
        )
        population = np.tensordot(
            powers,
            self.excited_population_orders[: selected_order + 1],
            axes=(0, 0),
        )
        return MeanFieldTrajectory(
            time=self.time,
            cavity_amplitude=np.asarray(cavity, dtype=complex),
            molecular_coherence=np.asarray(coherence, dtype=complex),
            excited_population=np.asarray(population, dtype=float),
        )


def _hierarchy_rhs(
    time: float,
    flat_state: NDArray[np.float64],
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    envelope: ComplexEnvelope,
    max_order: int,
) -> NDArray[np.float64]:
    state = flat_state.reshape(max_order + 1, 5)
    cavity = state[:, 0] + 1j * state[:, 1]
    coherence = state[:, 2] + 1j * state[:, 3]
    population = state[:, 4]

    derivative = np.zeros_like(state)
    delta_c = parameters.omega_c - drive_frequency
    delta_m = parameters.omega_0 - drive_frequency
    coupling = parameters.collective_coupling

    for order in range(1, max_order + 1):
        drive_term = envelope(time) if order == 1 else 0.0j
        cavity_derivative = -(
            0.5 * parameters.kappa + 1j * delta_c
        ) * cavity[order] - 1j * coupling * coherence[order] - drive_term

        cavity_population_convolution = sum(
            cavity[order - lower_order] * population[lower_order]
            for lower_order in range(1, order)
        )
        coherence_derivative = -(
            parameters.coherence_decay_rate + 1j * delta_m
        ) * coherence[order] - 1j * coupling * cavity[order]
        coherence_derivative += 2j * coupling * cavity_population_convolution

        field_coherence_convolution = sum(
            np.conjugate(cavity[lower_order]) * coherence[order - lower_order]
            for lower_order in range(1, order)
        )
        population_derivative = (
            -parameters.gamma * population[order]
            - 2.0 * coupling * np.imag(field_coherence_convolution)
        )

        derivative[order, 0] = cavity_derivative.real
        derivative[order, 1] = cavity_derivative.imag
        derivative[order, 2] = coherence_derivative.real
        derivative[order, 3] = coherence_derivative.imag
        derivative[order, 4] = population_derivative

    return derivative.ravel()


def propagate_perturbative_hierarchy(
    time: ArrayLike,
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    envelope: ComplexEnvelope,
    max_order: int = 3,
    rtol: float = 1e-10,
    atol: float = 1e-12,
    method: str = "DOP853",
) -> PerturbativeTrajectory:
    """Propagate coefficients from order zero through ``max_order``."""

    if max_order < 1:
        raise ValueError("max_order must be at least one")
    time_array = np.asarray(time, dtype=float)
    if time_array.ndim != 1 or time_array.size < 2 or np.any(np.diff(time_array) <= 0):
        raise ValueError("time must be a strictly increasing 1D grid")

    initial_state = np.zeros((max_order + 1, 5), dtype=float)
    solution = solve_ivp(
        lambda current_time, state: _hierarchy_rhs(
            current_time,
            state,
            parameters=parameters,
            drive_frequency=drive_frequency,
            envelope=envelope,
            max_order=max_order,
        ),
        (float(time_array[0]), float(time_array[-1])),
        initial_state.ravel(),
        t_eval=time_array,
        rtol=rtol,
        atol=atol,
        method=method,
    )
    if not solution.success:
        raise RuntimeError(f"perturbative propagation failed: {solution.message}")

    state = solution.y.T.reshape(time_array.size, max_order + 1, 5).transpose(1, 2, 0)
    return PerturbativeTrajectory(
        time=time_array,
        cavity_orders=state[:, 0] + 1j * state[:, 1],
        coherence_orders=state[:, 2] + 1j * state[:, 3],
        excited_population_orders=state[:, 4],
    )

