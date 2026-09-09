"""Molecular susceptibilities used by the linear-response theory.

The sign and linewidth conventions follow Yuen-Zhou and Koner,
J. Chem. Phys. 160, 154107 (2024), especially Eqs. (21), (37), and
(38): a passive transition has a positive imaginary susceptibility.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import simpson
from scipy.special import gammaln, wofz


def thermal_population_factor(delta_over_kbt: ArrayLike) -> NDArray[np.float64]:
    """Return the 2LS population difference ``tanh(DeltaE / 2 k_B T)``.

    ``delta_over_kbt`` is dimensionless.  Passing ``np.inf`` gives the
    zero-temperature value one without assigning any SI units.
    """

    ratio = np.asarray(delta_over_kbt, dtype=float)
    return np.tanh(0.5 * ratio)


def two_level_susceptibility(
    omega: ArrayLike,
    *,
    omega_exc: float,
    gamma: float,
    collective_coupling: float,
    population_factor: float = 1.0,
) -> NDArray[np.complex128]:
    """Susceptibility of identical two-level emitters, Eq. (37).

    Parameters
    ----------
    omega:
        Angular-frequency grid in arbitrary but consistent units.
    omega_exc:
        Molecular transition frequency.
    gamma:
        Full homogeneous molecular linewidth.  The retarded denominator
        contains ``+ i gamma/2``.
    collective_coupling:
        Zero-temperature collective coupling ``G = sqrt(N) |lambda mu|``.
    population_factor:
        Stationary population difference ``p_g-p_e``.  For a thermal 2LS
        this is ``tanh(DeltaE / 2 k_B T)``.

    Returns
    -------
    numpy.ndarray
        ``chi(omega) = -G^2 population_factor /
        (omega-omega_exc+i gamma/2)``.
    """

    if gamma < 0:
        raise ValueError("gamma must be non-negative")
    if collective_coupling < 0:
        raise ValueError("collective_coupling must be non-negative")

    omega_array = np.asarray(omega, dtype=float)
    denominator = omega_array - omega_exc + 0.5j * gamma
    return -(collective_coupling**2) * population_factor / denominator


def two_level_poles(
    *,
    omega_c: float,
    omega_exc: float,
    kappa: float,
    gamma: float,
    collective_coupling: float,
    population_factor: float = 1.0,
) -> NDArray[np.complex128]:
    """Return the two complex poles of the damped 2LS-cavity denominator."""

    effective_coupling_sq = collective_coupling**2 * population_factor
    cavity_pole = omega_c - 0.5j * kappa
    molecule_pole = omega_exc - 0.5j * gamma
    coefficients = [
        1.0,
        -(cavity_pole + molecule_pole),
        cavity_pole * molecule_pole - effective_coupling_sq,
    ]
    return np.sort_complex(np.roots(coefficients))


def gaussian_density(
    excitation_frequency: ArrayLike,
    *,
    mean: float,
    sigma: float,
) -> NDArray[np.float64]:
    """Normalized Gaussian distribution in linear-paper Eq. (42)."""

    if sigma <= 0:
        raise ValueError("sigma must be positive")
    frequency = np.asarray(excitation_frequency, dtype=float)
    return np.exp(-0.5 * ((frequency - mean) / sigma) ** 2) / (
        np.sqrt(2.0 * np.pi) * sigma
    )


def lorentzian_density(
    excitation_frequency: ArrayLike,
    *,
    mean: float,
    sigma: float,
    normalized: bool = True,
) -> NDArray[np.float64]:
    """Lorentzian disorder density associated with Eq. (43).

    The printed accepted manuscript has ``(sigma/2)^2`` in the numerator,
    which integrates to ``sigma/2`` and therefore is not a probability
    density.  The default implements the normalized Cauchy density, with
    HWHM ``sigma/2``.  ``normalized=False`` reproduces the printed typo for
    explicit diagnostics only.
    """

    if sigma <= 0:
        raise ValueError("sigma must be positive")
    frequency = np.asarray(excitation_frequency, dtype=float)
    half_width = sigma / 2.0
    numerator = half_width if normalized else half_width**2
    return numerator / (np.pi * ((frequency - mean) ** 2 + half_width**2))


def inhomogeneous_susceptibility_grid(
    omega: ArrayLike,
    excitation_grid: ArrayLike,
    density: ArrayLike,
    *,
    gamma: float,
    collective_coupling: float,
    chunk_size: int = 256,
) -> NDArray[np.complex128]:
    """Evaluate Eq. (41) by deterministic Simpson integration.

    The density is normalized numerically before integration.  This makes
    finite-range convergence transparent while preventing an integration
    range from silently changing the collective oscillator strength.
    """

    if gamma < 0:
        raise ValueError("gamma must be non-negative")
    if collective_coupling < 0:
        raise ValueError("collective_coupling must be non-negative")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    omega_array = np.asarray(omega, dtype=float)
    excitation = np.asarray(excitation_grid, dtype=float)
    probability = np.asarray(density, dtype=float)
    if excitation.ndim != 1 or probability.ndim != 1 or excitation.shape != probability.shape:
        raise ValueError("excitation_grid and density must be matching 1D arrays")
    if np.any(np.diff(excitation) <= 0):
        raise ValueError("excitation_grid must be strictly increasing")
    normalization = float(simpson(probability, x=excitation))
    if not np.isfinite(normalization) or normalization <= 0:
        raise ValueError("density has invalid numerical normalization")
    probability = probability / normalization

    flat_omega = omega_array.ravel()
    result = np.empty(flat_omega.size, dtype=complex)
    for start in range(0, flat_omega.size, chunk_size):
        stop = min(start + chunk_size, flat_omega.size)
        denominator = (
            flat_omega[start:stop, None]
            - excitation[None, :]
            + 0.5j * gamma
        )
        result[start:stop] = -(collective_coupling**2) * simpson(
            probability[None, :] / denominator,
            x=excitation,
            axis=1,
        )
    return result.reshape(omega_array.shape)


def gaussian_disorder_susceptibility(
    omega: ArrayLike,
    *,
    mean: float,
    sigma: float,
    gamma: float,
    collective_coupling: float,
) -> NDArray[np.complex128]:
    """Analytic Voigt evaluation of Eqs. (41)-(42).

    The result is benchmarked against ``inhomogeneous_susceptibility_grid``
    in the reproduction script and tests.
    """

    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if gamma <= 0:
        raise ValueError("gamma must be positive for the retarded Voigt expression")
    omega_array = np.asarray(omega, dtype=float)
    zeta = (omega_array - mean + 0.5j * gamma) / (np.sqrt(2.0) * sigma)
    integral = -1j * np.sqrt(np.pi) / (np.sqrt(2.0) * sigma) * wofz(zeta)
    return -(collective_coupling**2) * integral


def lorentzian_disorder_susceptibility(
    omega: ArrayLike,
    *,
    mean: float,
    sigma: float,
    gamma: float,
    collective_coupling: float,
) -> NDArray[np.complex128]:
    """Analytic convolution for normalized Lorentzian disorder.

    With Eq. (43)'s width convention, disorder adds ``sigma`` to the full
    homogeneous linewidth ``gamma``.
    """

    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if gamma < 0:
        raise ValueError("gamma must be non-negative")
    omega_array = np.asarray(omega, dtype=float)
    return -(collective_coupling**2) / (
        omega_array - mean + 0.5j * (gamma + sigma)
    )


def franck_condon_weights(
    huang_rhys: float,
    *,
    tolerance: float = 1e-12,
    m_max: int | None = None,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Return Poisson Franck-Condon weights ``exp(-S) S^m/m!``.

    If ``m_max`` is omitted, it is increased until the omitted probability
    is below ``tolerance``.
    """

    if huang_rhys < 0:
        raise ValueError("huang_rhys must be non-negative")
    if not 0 < tolerance < 1:
        raise ValueError("tolerance must lie between zero and one")
    if m_max is not None and m_max < 0:
        raise ValueError("m_max must be non-negative")

    if huang_rhys == 0:
        return np.array([0], dtype=int), np.array([1.0])

    chosen_m_max = 0 if m_max is None else int(m_max)
    while True:
        m = np.arange(chosen_m_max + 1, dtype=int)
        weights = np.exp(-huang_rhys + m * np.log(huang_rhys) - gammaln(m + 1.0))
        if m_max is not None or 1.0 - float(np.sum(weights)) <= tolerance:
            return m, weights
        chosen_m_max += 1
        if chosen_m_max > 100_000:
            raise RuntimeError("Franck-Condon truncation did not converge")


def vibronic_susceptibility(
    omega: ArrayLike,
    *,
    omega_exc: float,
    omega_v: float,
    huang_rhys: float,
    gamma: float,
    collective_coupling: float,
    m_max: int,
) -> NDArray[np.complex128]:
    """Vibronic 2LS susceptibility from linear-paper Eq. (44)."""

    if omega_v <= 0:
        raise ValueError("omega_v must be positive")
    m, weights = franck_condon_weights(huang_rhys, m_max=m_max)
    transitions = omega_exc - huang_rhys * omega_v + m * omega_v
    omega_array = np.asarray(omega, dtype=float)
    denominator = omega_array[..., None] - transitions + 0.5j * gamma
    return -(collective_coupling**2) * np.sum(weights / denominator, axis=-1)


def three_level_susceptibility(
    omega: ArrayLike,
    *,
    level_energies: ArrayLike,
    populations: ArrayLike,
    gamma: float,
    collective_coupling: float,
    relative_dipoles: ArrayLike | None = None,
) -> NDArray[np.complex128]:
    """Positive-frequency 3LS susceptibility from Eq. (46).

    ``collective_coupling`` is the coupling for a unit relative transition
    dipole.  Only uphill transitions ``i < j`` are retained, consistently
    with the RWA used to derive the spectra.
    """

    energies = np.asarray(level_energies, dtype=float)
    probability = np.asarray(populations, dtype=float)
    if energies.shape != (3,) or probability.shape != (3,):
        raise ValueError("level_energies and populations must each contain three values")
    if np.any(np.diff(energies) <= 0):
        raise ValueError("level energies must be strictly increasing")
    if np.any(probability < 0) or not np.isclose(np.sum(probability), 1.0, atol=1e-12):
        raise ValueError("populations must be non-negative and sum to one")
    if gamma < 0 or collective_coupling < 0:
        raise ValueError("gamma and collective_coupling must be non-negative")

    if relative_dipoles is None:
        dipoles = np.ones((3, 3), dtype=float)
    else:
        dipoles = np.asarray(relative_dipoles, dtype=float)
        if dipoles.shape != (3, 3):
            raise ValueError("relative_dipoles must be a 3x3 array")

    omega_array = np.asarray(omega, dtype=float)
    chi = np.zeros_like(omega_array, dtype=complex)
    for lower in range(3):
        for upper in range(lower + 1, 3):
            transition = energies[upper] - energies[lower]
            population_difference = probability[lower] - probability[upper]
            coupling_sq = collective_coupling**2 * abs(dipoles[upper, lower]) ** 2
            chi -= coupling_sq * population_difference / (
                omega_array - transition + 0.5j * gamma
            )
    return chi
