"""Photon Green function and input-output observables.

This module implements Eqs. (27), (32)-(34) of Yuen-Zhou and Koner,
J. Chem. Phys. 160, 154107 (2024).  Reflection and absorption are
evaluated independently from Eq. (27), so flux conservation is a real
numerical check rather than a consequence of defining ``R=1-T-A``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class LinearSpectra:
    omega: NDArray[np.float64]
    susceptibility: NDArray[np.complex128]
    green: NDArray[np.complex128]
    transmission: NDArray[np.float64]
    reflection: NDArray[np.float64]
    absorption: NDArray[np.float64]

    @property
    def flux_sum(self) -> NDArray[np.float64]:
        return self.transmission + self.reflection + self.absorption


def _validate_decay_rates(kappa_left: float, kappa_right: float) -> float:
    if kappa_left < 0 or kappa_right < 0:
        raise ValueError("mirror escape rates must be non-negative")
    kappa = kappa_left + kappa_right
    if kappa <= 0:
        raise ValueError("total cavity escape rate must be positive")
    return kappa


def photon_green_function(
    omega: ArrayLike,
    susceptibility: ArrayLike,
    *,
    omega_c: float,
    kappa: float,
) -> NDArray[np.complex128]:
    """Return ``D^R = 1/(omega-omega_c+i kappa/2+chi)``, Eq. (32)."""

    if kappa <= 0:
        raise ValueError("kappa must be positive")
    omega_array = np.asarray(omega, dtype=float)
    chi_array = np.asarray(susceptibility, dtype=complex)
    if omega_array.shape != chi_array.shape:
        raise ValueError("omega and susceptibility must have the same shape")
    return 1.0 / (omega_array - omega_c + 0.5j * kappa + chi_array)


def spectra_from_susceptibility(
    omega: ArrayLike,
    susceptibility: ArrayLike,
    *,
    omega_c: float,
    kappa_left: float,
    kappa_right: float,
) -> LinearSpectra:
    """Compute T, R, and A from a molecular susceptibility.

    The formulas correspond to incidence through the left mirror.  The
    total photon escape rate is ``kappa_left + kappa_right``.
    """

    kappa = _validate_decay_rates(kappa_left, kappa_right)
    omega_array = np.asarray(omega, dtype=float)
    chi_array = np.asarray(susceptibility, dtype=complex)
    green = photon_green_function(
        omega_array,
        chi_array,
        omega_c=omega_c,
        kappa=kappa,
    )

    green_abs_sq = np.abs(green) ** 2
    transmission = kappa_left * kappa_right * green_abs_sq
    reflection = 1.0 + 2.0 * kappa_left * np.imag(green) + kappa_left**2 * green_abs_sq
    absorption = -kappa_left * (kappa * green_abs_sq + 2.0 * np.imag(green))

    return LinearSpectra(
        omega=omega_array,
        susceptibility=chi_array,
        green=green,
        transmission=np.asarray(transmission, dtype=float),
        reflection=np.asarray(reflection, dtype=float),
        absorption=np.asarray(absorption, dtype=float),
    )

