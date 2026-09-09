from __future__ import annotations

import unittest

import numpy as np

from src.nonlinear.mean_field import TwoLevelMeanFieldParameters
from src.nonlinear.pump_probe import propagate_pump_probe_hierarchy


class TestPumpProbeHierarchy(unittest.TestCase):
    def test_required_orders_and_selection_rules(self) -> None:
        parameters = TwoLevelMeanFieldParameters(
            omega_c=0.0,
            omega_0=0.0,
            kappa=1.0,
            collective_coupling=3.0,
            gamma=0.2,
            gamma_phi=0.1,
        )
        time = np.linspace(0.0, 8.0, 2_001)
        pump = lambda t: np.exp(-0.5 * ((t - 2.0) / 0.2) ** 2)
        probe = lambda t: np.exp(-0.5 * ((t - 3.0) / 0.2) ** 2)
        hierarchy = propagate_pump_probe_hierarchy(
            time,
            parameters=parameters,
            drive_frequency=0.0,
            pump_envelope=pump,
            probe_envelope=probe,
        )
        self.assertGreater(np.max(np.abs(hierarchy.cavity_orders[0, 1])), 1e-3)
        self.assertGreater(np.max(np.abs(hierarchy.cavity_orders[2, 1])), 1e-5)
        self.assertLess(np.max(np.abs(hierarchy.cavity_orders[1, 1])), 1e-13)
        self.assertGreater(
            np.max(np.abs(hierarchy.excited_population_orders[2, 0])),
            1e-4,
        )


if __name__ == "__main__":
    unittest.main()

