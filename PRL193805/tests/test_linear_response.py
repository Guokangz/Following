from __future__ import annotations

import unittest

import numpy as np

from src.linear.input_output import spectra_from_susceptibility
from scipy.integrate import quad

from src.linear.susceptibility import (
    franck_condon_weights,
    gaussian_density,
    gaussian_disorder_susceptibility,
    inhomogeneous_susceptibility_grid,
    lorentzian_density,
    lorentzian_disorder_susceptibility,
    thermal_population_factor,
    three_level_susceptibility,
    two_level_susceptibility,
)


class TestLinearInputOutput(unittest.TestCase):
    def setUp(self) -> None:
        self.omega = np.linspace(-5.0, 5.0, 20_001)
        self.kappa = 0.1
        self.kappa_left = self.kappa_right = self.kappa / 2.0

    def test_empty_balanced_cavity(self) -> None:
        spectra = spectra_from_susceptibility(
            self.omega,
            np.zeros_like(self.omega, dtype=complex),
            omega_c=0.0,
            kappa_left=self.kappa_left,
            kappa_right=self.kappa_right,
        )
        expected_t = (self.kappa / 2.0) ** 2 / (
            self.omega**2 + (self.kappa / 2.0) ** 2
        )
        self.assertLess(np.max(np.abs(spectra.transmission - expected_t)), 2e-15)
        self.assertLess(np.max(np.abs(spectra.absorption)), 2e-15)
        self.assertLess(np.max(np.abs(spectra.flux_sum - 1.0)), 2e-15)
        center = self.omega.size // 2
        self.assertAlmostEqual(spectra.transmission[center], 1.0, places=14)
        self.assertAlmostEqual(spectra.reflection[center], 0.0, places=14)

    def test_fig3a_flux_conservation_and_passivity(self) -> None:
        chi = two_level_susceptibility(
            self.omega,
            omega_exc=0.0,
            gamma=0.3,
            collective_coupling=2.0,
        )
        spectra = spectra_from_susceptibility(
            self.omega,
            chi,
            omega_c=0.0,
            kappa_left=self.kappa_left,
            kappa_right=self.kappa_right,
        )
        self.assertLess(np.max(np.abs(spectra.flux_sum - 1.0)), 5e-14)
        self.assertGreaterEqual(np.min(spectra.transmission), -2e-14)
        self.assertGreaterEqual(np.min(spectra.reflection), -2e-14)
        self.assertGreaterEqual(np.min(spectra.absorption), -2e-14)

    def test_eq38_matches_general_eq34_interface(self) -> None:
        gamma = 0.3
        collective_coupling = 2.0
        chi = two_level_susceptibility(
            self.omega,
            omega_exc=0.0,
            gamma=gamma,
            collective_coupling=collective_coupling,
        )
        spectra = spectra_from_susceptibility(
            self.omega,
            chi,
            omega_c=0.0,
            kappa_left=self.kappa_left,
            kappa_right=self.kappa_right,
        )
        denominator = (
            (self.omega + 0.5j * self.kappa)
            * (self.omega + 0.5j * gamma)
            - collective_coupling**2
        )
        expected_t = (
            self.kappa_left
            * self.kappa_right
            * (self.omega**2 + (gamma / 2.0) ** 2)
            / np.abs(denominator) ** 2
        )
        expected_a = (
            self.kappa_left
            * gamma
            * collective_coupling**2
            / np.abs(denominator) ** 2
        )
        self.assertLess(np.max(np.abs(spectra.transmission - expected_t)), 2e-14)
        self.assertLess(np.max(np.abs(spectra.absorption - expected_a)), 2e-14)

    def test_no_molecular_damping_has_no_absorption(self) -> None:
        # Avoid exactly omega=omega_exc, where the zero-width susceptibility
        # is distribution-valued rather than an ordinary finite function.
        omega = self.omega + 1.0e-7
        chi = two_level_susceptibility(
            omega,
            omega_exc=0.0,
            gamma=0.0,
            collective_coupling=2.0,
        )
        spectra = spectra_from_susceptibility(
            omega,
            chi,
            omega_c=0.0,
            kappa_left=self.kappa_left,
            kappa_right=self.kappa_right,
        )
        self.assertLess(np.max(np.abs(spectra.absorption)), 5e-12)
        self.assertLess(np.max(np.abs(spectra.flux_sum - 1.0)), 5e-12)

    def test_thermal_population_factor(self) -> None:
        ratios = np.array([np.inf, 1.1, 0.1, 0.0])
        factors = thermal_population_factor(ratios)
        np.testing.assert_allclose(
            factors,
            [1.0, np.tanh(0.55), np.tanh(0.05), 0.0],
            rtol=0.0,
            atol=1e-15,
        )

    def test_disorder_distributions_are_normalized(self) -> None:
        gaussian_integral = quad(
            lambda x: float(gaussian_density(x, mean=0.0, sigma=1.0)),
            -np.inf,
            np.inf,
        )[0]
        lorentzian_integral = quad(
            lambda x: float(lorentzian_density(x, mean=0.0, sigma=1.0)),
            -np.inf,
            np.inf,
        )[0]
        printed_eq43_integral = quad(
            lambda x: float(
                lorentzian_density(x, mean=0.0, sigma=1.0, normalized=False)
            ),
            -np.inf,
            np.inf,
        )[0]
        self.assertAlmostEqual(gaussian_integral, 1.0, places=12)
        self.assertAlmostEqual(lorentzian_integral, 1.0, places=12)
        self.assertAlmostEqual(printed_eq43_integral, 0.5, places=12)

    def test_gaussian_grid_matches_voigt_expression(self) -> None:
        omega = np.linspace(-4.0, 4.0, 101)
        excitation = np.linspace(-8.0, 8.0, 16_001)
        density = gaussian_density(excitation, mean=0.0, sigma=1.0)
        grid_chi = inhomogeneous_susceptibility_grid(
            omega,
            excitation,
            density,
            gamma=0.1,
            collective_coupling=1.5,
        )
        analytic_chi = gaussian_disorder_susceptibility(
            omega,
            mean=0.0,
            sigma=1.0,
            gamma=0.1,
            collective_coupling=1.5,
        )
        self.assertLess(np.max(np.abs(grid_chi - analytic_chi)), 2e-11)

    def test_lorentzian_disorder_adds_linewidth(self) -> None:
        chi = lorentzian_disorder_susceptibility(
            self.omega,
            mean=0.0,
            sigma=1.0,
            gamma=0.1,
            collective_coupling=1.5,
        )
        expected = two_level_susceptibility(
            self.omega,
            omega_exc=0.0,
            gamma=1.1,
            collective_coupling=1.5,
        )
        np.testing.assert_allclose(chi, expected, rtol=2e-15, atol=2e-15)

    def test_franck_condon_truncation_and_mean(self) -> None:
        m, weights = franck_condon_weights(3.0, tolerance=1e-12)
        self.assertLess(1.0 - np.sum(weights), 1e-12)
        self.assertAlmostEqual(float(np.sum(m * weights)), 3.0, places=9)

    def test_equal_three_level_populations_are_transparent(self) -> None:
        chi = three_level_susceptibility(
            self.omega,
            level_energies=[0.0, 1.0, 3.0],
            populations=[1 / 3, 1 / 3, 1 / 3],
            gamma=0.3,
            collective_coupling=1.0,
        )
        self.assertLess(np.max(np.abs(chi)), 1e-15)


if __name__ == "__main__":
    unittest.main()
