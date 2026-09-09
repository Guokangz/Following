"""Two-pulse perturbative hierarchy for a semiclassical 2LS cavity.

This is the reduced-Bloch-variable form of SM Eqs. (S.7)-(S.10), with
direct propagation of every coefficient ``(n)(m)`` rather than amplitude
fitting.
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
class PumpProbeTrajectory:
    time: NDArray[np.float64]
    cavity_orders: NDArray[np.complex128]
    coherence_orders: NDArray[np.complex128]
    excited_population_orders: NDArray[np.float64]

    @property
    def max_pump_order(self) -> int:
        return self.cavity_orders.shape[0] - 1

    @property
    def max_probe_order(self) -> int:
        return self.cavity_orders.shape[1] - 1

    def reconstruct(self, eta_pump: float, eta_probe: float) -> MeanFieldTrajectory:
        pump_powers = eta_pump ** np.arange(self.max_pump_order + 1)
        probe_powers = eta_probe ** np.arange(self.max_probe_order + 1)
        weights = pump_powers[:, None] * probe_powers[None, :]
        cavity = np.tensordot(weights, self.cavity_orders, axes=((0, 1), (0, 1)))
        coherence = np.tensordot(weights, self.coherence_orders, axes=((0, 1), (0, 1)))
        population = np.tensordot(
            weights,
            self.excited_population_orders,
            axes=((0, 1), (0, 1)),
        )
        return MeanFieldTrajectory(
            time=self.time,
            cavity_amplitude=np.asarray(cavity, dtype=complex),
            molecular_coherence=np.asarray(coherence, dtype=complex),
            excited_population=np.asarray(population, dtype=float),
        )


def _pump_probe_rhs(
    time: float,
    flat_state: NDArray[np.float64],
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    pump_envelope: ComplexEnvelope,
    probe_envelope: ComplexEnvelope,
    max_pump_order: int,
    max_probe_order: int,
) -> NDArray[np.float64]:
    shape = (max_pump_order + 1, max_probe_order + 1, 5)
    state = flat_state.reshape(shape)
    cavity = state[..., 0] + 1j * state[..., 1]
    coherence = state[..., 2] + 1j * state[..., 3]
    population = state[..., 4]
    derivative = np.zeros_like(state)

    delta_c = parameters.omega_c - drive_frequency
    delta_m = parameters.omega_0 - drive_frequency
    coupling = parameters.collective_coupling

    for pump_order in range(max_pump_order + 1):
        for probe_order in range(max_probe_order + 1):
            if pump_order == 0 and probe_order == 0:
                continue
            drive_term = 0.0j
            if pump_order == 1 and probe_order == 0:
                drive_term += pump_envelope(time)
            if pump_order == 0 and probe_order == 1:
                drive_term += probe_envelope(time)

            cavity_derivative = -(
                0.5 * parameters.kappa + 1j * delta_c
            ) * cavity[pump_order, probe_order]
            cavity_derivative -= 1j * coupling * coherence[pump_order, probe_order]
            cavity_derivative -= drive_term

            cavity_population_convolution = 0.0j
            field_coherence_convolution = 0.0j
            for pump_lower in range(pump_order + 1):
                for probe_lower in range(probe_order + 1):
                    pump_upper = pump_order - pump_lower
                    probe_upper = probe_order - probe_lower
                    cavity_population_convolution += (
                        cavity[pump_upper, probe_upper]
                        * population[pump_lower, probe_lower]
                    )
                    field_coherence_convolution += (
                        np.conjugate(cavity[pump_lower, probe_lower])
                        * coherence[pump_upper, probe_upper]
                    )

            coherence_derivative = -(
                parameters.coherence_decay_rate + 1j * delta_m
            ) * coherence[pump_order, probe_order]
            coherence_derivative -= 1j * coupling * cavity[pump_order, probe_order]
            coherence_derivative += 2j * coupling * cavity_population_convolution

            population_derivative = (
                -parameters.gamma * population[pump_order, probe_order]
                - 2.0 * coupling * np.imag(field_coherence_convolution)
            )

            derivative[pump_order, probe_order, 0] = cavity_derivative.real
            derivative[pump_order, probe_order, 1] = cavity_derivative.imag
            derivative[pump_order, probe_order, 2] = coherence_derivative.real
            derivative[pump_order, probe_order, 3] = coherence_derivative.imag
            derivative[pump_order, probe_order, 4] = population_derivative
    return derivative.ravel()


def propagate_pump_probe_hierarchy(
    time: ArrayLike,
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    pump_envelope: ComplexEnvelope,
    probe_envelope: ComplexEnvelope,
    max_pump_order: int = 2,
    max_probe_order: int = 1,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    method: str = "DOP853",
) -> PumpProbeTrajectory:
    """Propagate all two-pulse coefficients up to ``(max_pump,max_probe)``."""

    if max_pump_order < 0 or max_probe_order < 0 or max_pump_order + max_probe_order < 1:
        raise ValueError("at least one positive pulse order is required")
    time_array = np.asarray(time, dtype=float)
    if time_array.ndim != 1 or time_array.size < 2 or np.any(np.diff(time_array) <= 0):
        raise ValueError("time must be a strictly increasing 1D grid")

    initial = np.zeros((max_pump_order + 1, max_probe_order + 1, 5), dtype=float)
    solution = solve_ivp(
        lambda current_time, state: _pump_probe_rhs(
            current_time,
            state,
            parameters=parameters,
            drive_frequency=drive_frequency,
            pump_envelope=pump_envelope,
            probe_envelope=probe_envelope,
            max_pump_order=max_pump_order,
            max_probe_order=max_probe_order,
        ),
        (float(time_array[0]), float(time_array[-1])),
        initial.ravel(),
        t_eval=time_array,
        rtol=rtol,
        atol=atol,
        method=method,
    )
    if not solution.success:
        raise RuntimeError(f"pump-probe hierarchy propagation failed: {solution.message}")

    state = solution.y.T.reshape(time_array.size, *initial.shape).transpose(1, 2, 3, 0)
    return PumpProbeTrajectory(
        time=time_array,
        cavity_orders=state[..., 0, :] + 1j * state[..., 1, :],
        coherence_orders=state[..., 2, :] + 1j * state[..., 3, :],
        excited_population_orders=state[..., 4, :],
    )

