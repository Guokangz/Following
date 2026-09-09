"""Explicit phase-resolved 2LS pump-probe hierarchy.

This module specializes SM Eqs. (S.39)-(S.43) to a two-level system and
propagates the two third-order positive-frequency cavity components
``(0,1)`` and ``(2,-1)``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from .mean_field import TwoLevelMeanFieldParameters


ComplexEnvelope = Callable[[float], complex]


@dataclass(frozen=True)
class PhaseResolvedTrajectory:
    time: NDArray[np.float64]
    alpha_pump: NDArray[np.complex128]
    coherence_pump: NDArray[np.complex128]
    alpha_probe: NDArray[np.complex128]
    coherence_probe: NDArray[np.complex128]
    pump_population_00: NDArray[np.float64]
    mixed_population_1m1: NDArray[np.complex128]
    alpha_3_01: NDArray[np.complex128]
    coherence_3_01: NDArray[np.complex128]
    alpha_3_2m1: NDArray[np.complex128]
    coherence_3_2m1: NDArray[np.complex128]

    @property
    def alpha_3_total_zero_phases(self) -> NDArray[np.complex128]:
        return self.alpha_3_01 + self.alpha_3_2m1


_VARIABLE_COUNT = 10


def _pack_complex(values: NDArray[np.complex128]) -> NDArray[np.float64]:
    return np.concatenate((values.real, values.imag))


def _unpack_complex(values: NDArray[np.float64]) -> NDArray[np.complex128]:
    return values[:_VARIABLE_COUNT] + 1j * values[_VARIABLE_COUNT:]


def _phase_rhs(
    time: float,
    real_state: NDArray[np.float64],
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    pump_envelope: ComplexEnvelope,
    probe_envelope: ComplexEnvelope,
) -> NDArray[np.float64]:
    (
        alpha_p,
        coherence_p,
        alpha_q,
        coherence_q,
        population_20,
        population_11_1m1,
        alpha_3_01,
        coherence_3_01,
        alpha_3_2m1,
        coherence_3_2m1,
    ) = _unpack_complex(real_state)

    delta_c = parameters.omega_c - drive_frequency
    delta_m = parameters.omega_0 - drive_frequency
    cavity_rate = 0.5 * parameters.kappa + 1j * delta_c
    molecular_rate = parameters.coherence_decay_rate + 1j * delta_m
    coupling = parameters.collective_coupling

    alpha_p_dot = -cavity_rate * alpha_p - 1j * coupling * coherence_p - pump_envelope(time)
    coherence_p_dot = -molecular_rate * coherence_p - 1j * coupling * alpha_p
    alpha_q_dot = -cavity_rate * alpha_q - 1j * coupling * coherence_q - probe_envelope(time)
    coherence_q_dot = -molecular_rate * coherence_q - 1j * coupling * alpha_q

    population_20_dot = (
        -parameters.gamma * population_20
        - 2.0 * coupling * np.imag(np.conjugate(alpha_p) * coherence_p)
    )
    population_11_1m1_dot = (
        -parameters.gamma * population_11_1m1
        - 1j * coupling * alpha_p * np.conjugate(coherence_q)
        + 1j * coupling * np.conjugate(alpha_q) * coherence_p
    )

    alpha_3_01_dot = -cavity_rate * alpha_3_01 - 1j * coupling * coherence_3_01
    coherence_3_01_dot = -molecular_rate * coherence_3_01 - 1j * coupling * alpha_3_01
    coherence_3_01_dot += 2j * coupling * (
        alpha_q * population_20
        + alpha_p * np.conjugate(population_11_1m1)
    )

    alpha_3_2m1_dot = -cavity_rate * alpha_3_2m1 - 1j * coupling * coherence_3_2m1
    coherence_3_2m1_dot = -molecular_rate * coherence_3_2m1 - 1j * coupling * alpha_3_2m1
    coherence_3_2m1_dot += 2j * coupling * alpha_p * population_11_1m1

    derivative = np.array(
        [
            alpha_p_dot,
            coherence_p_dot,
            alpha_q_dot,
            coherence_q_dot,
            population_20_dot,
            population_11_1m1_dot,
            alpha_3_01_dot,
            coherence_3_01_dot,
            alpha_3_2m1_dot,
            coherence_3_2m1_dot,
        ],
        dtype=complex,
    )
    return _pack_complex(derivative)


def propagate_phase_resolved_hierarchy(
    time: ArrayLike,
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    pump_envelope: ComplexEnvelope,
    probe_envelope: ComplexEnvelope,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    method: str = "DOP853",
) -> PhaseResolvedTrajectory:
    """Propagate explicit ``(0,1)`` and ``(2,-1)`` third-order components."""

    time_array = np.asarray(time, dtype=float)
    if time_array.ndim != 1 or time_array.size < 2 or np.any(np.diff(time_array) <= 0):
        raise ValueError("time must be a strictly increasing 1D grid")
    initial = np.zeros(2 * _VARIABLE_COUNT, dtype=float)
    solution = solve_ivp(
        lambda current_time, state: _phase_rhs(
            current_time,
            state,
            parameters=parameters,
            drive_frequency=drive_frequency,
            pump_envelope=pump_envelope,
            probe_envelope=probe_envelope,
        ),
        (float(time_array[0]), float(time_array[-1])),
        initial,
        t_eval=time_array,
        rtol=rtol,
        atol=atol,
        method=method,
    )
    if not solution.success:
        raise RuntimeError(f"phase-resolved propagation failed: {solution.message}")
    complex_state = solution.y[:_VARIABLE_COUNT] + 1j * solution.y[_VARIABLE_COUNT:]
    return PhaseResolvedTrajectory(
        time=time_array,
        alpha_pump=complex_state[0],
        coherence_pump=complex_state[1],
        alpha_probe=complex_state[2],
        coherence_probe=complex_state[3],
        pump_population_00=np.real(complex_state[4]),
        mixed_population_1m1=complex_state[5],
        alpha_3_01=complex_state[6],
        coherence_3_01=complex_state[7],
        alpha_3_2m1=complex_state[8],
        coherence_3_2m1=complex_state[9],
    )

