from __future__ import annotations

import unittest

import numpy as np

from src.nonlinear.mean_field import TwoLevelMeanFieldParameters
from src.nonlinear.phase_components import propagate_phase_resolved_hierarchy
from src.nonlinear.pump_probe import propagate_pump_probe_hierarchy


class TestPhaseComponents(unittest.TestCase):
    def test_components_sum_to_unresolved_alpha_21(self) -> None:
        parameters = TwoLevelMeanFieldParameters(
            omega_c=0.0,
            omega_0=0.0,
            kappa=1.0,
            collective_coupling=3.0,
            gamma=0.0,
            gamma_phi=0.1,
        )
        time = np.linspace(0.0, 10.0, 2_501)
        pump = lambda t: np.exp(-0.5 * ((t - 2.0) / 0.1) ** 2)
        probe = lambda t: np.exp(-0.5 * ((t - 3.0) / 0.1) ** 2)
        explicit = propagate_phase_resolved_hierarchy(
            time,
            parameters=parameters,
            drive_frequency=0.0,
            pump_envelope=pump,
            probe_envelope=probe,
            rtol=1e-10,
            atol=1e-12,
        )
        unresolved = propagate_pump_probe_hierarchy(
            time,
            parameters=parameters,
            drive_frequency=0.0,
            pump_envelope=pump,
            probe_envelope=probe,
            rtol=1e-10,
            atol=1e-12,
        )
        self.assertLess(
            np.max(
                np.abs(
                    explicit.alpha_3_total_zero_phases
                    - unresolved.cavity_orders[2, 1]
                )
            ),
            2e-10,
        )


if __name__ == "__main__":
    unittest.main()

