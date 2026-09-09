"""Transparent two-level semiclassical mean-field equations.

The equations are the rotating-wave, rotating-frame form of PRL Eqs. (2)
and (3), with the dissipative convention made explicit in SM Eqs. (S.14)
and (S.39).  The cavity amplitude is scaled by ``sqrt(N)`` so the equations
depend on the collective coupling ``G=g sqrt(N)`` only:

    da/dt = -(kappa/2 + i delta_c) a - i G s - u(t)
    ds/dt = -(Gamma_2 + i delta_m) s - i G a (1-2p)
    dp/dt = -gamma p - 2 G Im(a* s)

Here ``s=rho_eg``, ``p=rho_ee``, and
``Gamma_2=(gamma+gamma_phi)/2`` under the Supplemental Material's rate
convention.  The physical ratios ``alpha/input`` and scaled ``a/u`` are
identical, so input-output spectra are unaffected by the scaling.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp


ComplexDrive = Callable[[float], complex]


def two_level_operators(
    omega_0: float,
    *,
    transition_dipole: float = 1.0,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Construct ``H0`` and ``mu`` in the ordered basis ``(|g>,|e>)``."""

    hamiltonian = np.array([[0.0, 0.0], [0.0, omega_0]], dtype=complex)
    dipole = transition_dipole * np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    return hamiltonian, dipole


@dataclass(frozen=True)
class TwoLevelMeanFieldParameters:
    omega_c: float
    omega_0: float
    kappa: float
    collective_coupling: float
    gamma: float = 0.0
    gamma_phi: float = 0.0

    def __post_init__(self) -> None:
        if self.kappa <= 0:
            raise ValueError("kappa must be positive")
        if self.collective_coupling < 0:
            raise ValueError("collective_coupling must be non-negative")
        if self.gamma < 0 or self.gamma_phi < 0:
            raise ValueError("molecular rates must be non-negative")

    @property
    def coherence_decay_rate(self) -> float:
        """Coherence HWHM ``(gamma+gamma_phi)/2`` in the SM convention."""

        return 0.5 * (self.gamma + self.gamma_phi)

    @property
    def linear_full_linewidth(self) -> float:
        """The ``gamma`` entering linear-paper denominator ``+i gamma/2``."""

        return self.gamma + self.gamma_phi


@dataclass(frozen=True)
class MeanFieldTrajectory:
    time: NDArray[np.float64]
    cavity_amplitude: NDArray[np.complex128]
    molecular_coherence: NDArray[np.complex128]
    excited_population: NDArray[np.float64]

    @property
    def ground_population(self) -> NDArray[np.float64]:
        return 1.0 - self.excited_population

    @property
    def trace(self) -> NDArray[np.float64]:
        return self.ground_population + self.excited_population

    @property
    def minimum_density_eigenvalue(self) -> NDArray[np.float64]:
        half_difference = 0.5 - self.excited_population
        radius = np.sqrt(half_difference**2 + np.abs(self.molecular_coherence) ** 2)
        return 0.5 - radius


def _pack_state(cavity: complex, coherence: complex, population: float) -> NDArray[np.float64]:
    return np.array(
        [cavity.real, cavity.imag, coherence.real, coherence.imag, population],
        dtype=float,
    )


def _unpack_state(state: ArrayLike) -> tuple[complex, complex, float]:
    values = np.asarray(state, dtype=float)
    return values[0] + 1j * values[1], values[2] + 1j * values[3], float(values[4])


def mean_field_rhs(
    time: float,
    state: NDArray[np.float64],
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    drive: ComplexDrive,
) -> NDArray[np.float64]:
    """Right-hand side of the rotating-frame nonlinear mean-field ODE."""

    cavity, coherence, population = _unpack_state(state)
    delta_c = parameters.omega_c - drive_frequency
    delta_m = parameters.omega_0 - drive_frequency
    coupling = parameters.collective_coupling

    cavity_derivative = -(
        0.5 * parameters.kappa + 1j * delta_c
    ) * cavity - 1j * coupling * coherence - drive(time)
    coherence_derivative = -(
        parameters.coherence_decay_rate + 1j * delta_m
    ) * coherence - 1j * coupling * cavity * (1.0 - 2.0 * population)
    population_derivative = (
        -parameters.gamma * population
        - 2.0 * coupling * np.imag(np.conjugate(cavity) * coherence)
    )
    return _pack_state(cavity_derivative, coherence_derivative, population_derivative)


def propagate_mean_field(
    time: ArrayLike,
    *,
    parameters: TwoLevelMeanFieldParameters,
    drive_frequency: float,
    drive: ComplexDrive,
    initial_cavity: complex = 0.0j,
    initial_coherence: complex = 0.0j,
    initial_excited_population: float = 0.0,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    method: str = "DOP853",
) -> MeanFieldTrajectory:
    """Propagate the full nonlinear mean-field equations with ``solve_ivp``."""

    time_array = np.asarray(time, dtype=float)
    if time_array.ndim != 1 or time_array.size < 2 or np.any(np.diff(time_array) <= 0):
        raise ValueError("time must be a strictly increasing 1D grid")
    if not 0.0 <= initial_excited_population <= 1.0:
        raise ValueError("initial_excited_population must lie in [0,1]")

    initial_state = _pack_state(
        initial_cavity,
        initial_coherence,
        initial_excited_population,
    )
    solution = solve_ivp(
        lambda current_time, state: mean_field_rhs(
            current_time,
            state,
            parameters=parameters,
            drive_frequency=drive_frequency,
            drive=drive,
        ),
        (float(time_array[0]), float(time_array[-1])),
        initial_state,
        t_eval=time_array,
        rtol=rtol,
        atol=atol,
        method=method,
    )
    if not solution.success:
        raise RuntimeError(f"mean-field propagation failed: {solution.message}")

    cavity = solution.y[0] + 1j * solution.y[1]
    coherence = solution.y[2] + 1j * solution.y[3]
    population = solution.y[4]
    return MeanFieldTrajectory(
        time=time_array,
        cavity_amplitude=cavity,
        molecular_coherence=coherence,
        excited_population=population,
    )


def linear_cw_intracavity_response(
    drive_frequency: ArrayLike,
    *,
    parameters: TwoLevelMeanFieldParameters,
) -> NDArray[np.complex128]:
    """Return the first-order steady response ``a/u`` to a CW drive.

    The cavity equation uses the PRL sign ``... - input``.  Consequently,
    ``a/u = -i D^R`` and ``i a/u`` equals the linear-paper photon Green
    function exactly when linewidth conventions are matched.
    """

    frequency = np.asarray(drive_frequency, dtype=float)
    delta_c = parameters.omega_c - frequency
    delta_m = parameters.omega_0 - frequency
    molecular_denominator = parameters.coherence_decay_rate + 1j * delta_m
    return -1.0 / (
        0.5 * parameters.kappa
        + 1j * delta_c
        + parameters.collective_coupling**2 / molecular_denominator
    )
