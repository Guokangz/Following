from __future__ import annotations

import unittest

import numpy as np
from scipy.signal import find_peaks

from src.nonlinear.mean_field import (
    TwoLevelMeanFieldParameters,
    linear_cw_intracavity_response,
    propagate_mean_field,
    two_level_operators,
)


class TestTwoLevelMeanField(unittest.TestCase):
    def setUp(self) -> None:
        self.parameters = TwoLevelMeanFieldParameters(
            omega_c=0.0,
            omega_0=0.0,
            kappa=0.1,
            collective_coupling=2.0,
            gamma=0.1,
            gamma_phi=0.2,
        )

    def test_two_level_operators_are_hermitian(self) -> None:
        hamiltonian, dipole = two_level_operators(2.3, transition_dipole=0.7)
        np.testing.assert_allclose(hamiltonian, hamiltonian.conj().T)
        np.testing.assert_allclose(dipole, dipole.conj().T)
        np.testing.assert_allclose(np.linalg.eigvalsh(hamiltonian), [0.0, 2.3])
        self.assertAlmostEqual(abs(dipole[1, 0]), 0.7)

    def test_zero_drive_ground_state_is_stationary(self) -> None:
        time = np.linspace(0.0, 100.0, 501)
        trajectory = propagate_mean_field(
            time,
            parameters=self.parameters,
            drive_frequency=0.0,
            drive=lambda _: 0.0j,
        )
        self.assertEqual(np.max(np.abs(trajectory.cavity_amplitude)), 0.0)
        self.assertEqual(np.max(np.abs(trajectory.molecular_coherence)), 0.0)
        self.assertEqual(np.max(np.abs(trajectory.excited_population)), 0.0)
        self.assertEqual(np.max(np.abs(trajectory.trace - 1.0)), 0.0)

    def test_linear_cw_splitting_is_twice_collective_coupling(self) -> None:
        omega = np.linspace(-4.0, 4.0, 20_001)
        response = linear_cw_intracavity_response(omega, parameters=self.parameters)
        transmission = (self.parameters.kappa / 2.0) ** 2 * np.abs(response) ** 2
        peaks, _ = find_peaks(transmission, prominence=0.01 * np.max(transmission))
        strongest = np.sort(peaks[np.argsort(transmission[peaks])[-2:]])
        splitting = omega[strongest[1]] - omega[strongest[0]]
        self.assertAlmostEqual(splitting, 4.0, delta=0.02)

    def test_weak_full_cw_solution_matches_linear_response(self) -> None:
        drive_amplitude = 1e-6
        time = np.linspace(0.0, 800.0, 4_001)
        for frequency in (-2.0, 0.0, 2.0):
            trajectory = propagate_mean_field(
                time,
                parameters=self.parameters,
                drive_frequency=frequency,
                drive=lambda _, amplitude=drive_amplitude: amplitude,
                rtol=1e-12,
                atol=1e-15,
            )
            numerical = trajectory.cavity_amplitude[-1] / drive_amplitude
            expected = linear_cw_intracavity_response(frequency, parameters=self.parameters)
            self.assertAlmostEqual(numerical.real, expected.real, delta=2e-8)
            self.assertAlmostEqual(numerical.imag, expected.imag, delta=2e-8)
            self.assertLess(abs(trajectory.trace[-1] - 1.0), 1e-14)
            self.assertGreaterEqual(trajectory.excited_population[-1], -1e-12)
            self.assertLessEqual(trajectory.excited_population[-1], 1.0 + 1e-12)
            self.assertGreaterEqual(trajectory.minimum_density_eigenvalue[-1], -1e-10)


if __name__ == "__main__":
    unittest.main()
