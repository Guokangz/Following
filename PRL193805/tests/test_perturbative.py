from __future__ import annotations

import unittest

import numpy as np

from src.nonlinear.mean_field import TwoLevelMeanFieldParameters, propagate_mean_field
from src.nonlinear.perturbative import propagate_perturbative_hierarchy


class TestPerturbativeHierarchy(unittest.TestCase):
    def setUp(self) -> None:
        self.parameters = TwoLevelMeanFieldParameters(
            omega_c=0.0,
            omega_0=0.0,
            kappa=1.0,
            collective_coupling=3.0,
            gamma=0.2,
            gamma_phi=0.1,
        )
        self.time = np.linspace(0.0, 8.0, 2_001)
        self.envelope = lambda t: np.exp(-0.5 * ((t - 2.0) / 0.2) ** 2)

    def test_order_selection_rules(self) -> None:
        hierarchy = propagate_perturbative_hierarchy(
            self.time,
            parameters=self.parameters,
            drive_frequency=0.0,
            envelope=self.envelope,
            max_order=3,
        )
        self.assertEqual(np.max(np.abs(hierarchy.cavity_orders[0])), 0.0)
        self.assertEqual(np.max(np.abs(hierarchy.coherence_orders[0])), 0.0)
        self.assertEqual(np.max(np.abs(hierarchy.excited_population_orders[0])), 0.0)
        self.assertLess(np.max(np.abs(hierarchy.cavity_orders[2])), 1e-14)
        self.assertLess(np.max(np.abs(hierarchy.coherence_orders[2])), 1e-14)
        self.assertLess(np.max(np.abs(hierarchy.excited_population_orders[1])), 1e-14)
        self.assertLess(np.max(np.abs(hierarchy.excited_population_orders[3])), 1e-14)

    def test_third_order_reconstruction_matches_small_full_solution(self) -> None:
        hierarchy = propagate_perturbative_hierarchy(
            self.time,
            parameters=self.parameters,
            drive_frequency=0.0,
            envelope=self.envelope,
            max_order=3,
            rtol=1e-11,
            atol=1e-13,
        )
        eta = 0.05
        reconstruction = hierarchy.reconstruct(eta)
        full = propagate_mean_field(
            self.time,
            parameters=self.parameters,
            drive_frequency=0.0,
            drive=lambda t: eta * self.envelope(t),
            rtol=1e-11,
            atol=1e-13,
        )
        total_error = max(
            np.max(np.abs(full.cavity_amplitude - reconstruction.cavity_amplitude)),
            np.max(np.abs(full.molecular_coherence - reconstruction.molecular_coherence)),
            np.max(np.abs(full.excited_population - reconstruction.excited_population)),
        )
        # The leading omitted population term is fourth order, so an O(eta^4)
        # error is expected even though the field is reconstructed to n=3.
        self.assertLess(total_error, 1e-7)


if __name__ == "__main__":
    unittest.main()
