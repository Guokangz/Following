#!/usr/bin/env python3
"""Benchmark explicit phase hierarchy against numerical phase cycling."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.nonlinear.mean_field import TwoLevelMeanFieldParameters
from src.nonlinear.phase_components import propagate_phase_resolved_hierarchy
from src.nonlinear.pump_probe import propagate_pump_probe_hierarchy


def main() -> None:
    figure_dir = PROJECT_ROOT / "figures" / "nonlinear"
    data_dir = PROJECT_ROOT / "data" / "nonlinear"
    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    parameters = TwoLevelMeanFieldParameters(
        omega_c=0.0,
        omega_0=0.0,
        kappa=1.0,
        collective_coupling=3.0,
        gamma=0.0,
        gamma_phi=0.1,
    )
    time = np.linspace(0.0, 10.0, 2_501)
    pump_center = 2.0
    probe_center = 3.0
    pulse_width = 0.1
    pump = lambda t: np.exp(-0.5 * ((t - pump_center) / pulse_width) ** 2)
    probe = lambda t: np.exp(-0.5 * ((t - probe_center) / pulse_width) ** 2)

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
    explicit_sum_error = float(
        np.max(
            np.abs(
                explicit.alpha_3_total_zero_phases
                - unresolved.cavity_orders[2, 1]
            )
        )
    )

    phase_count = 6
    phases = 2.0 * np.pi * np.arange(phase_count) / phase_count
    phased_alpha_21 = np.empty((phase_count, phase_count, time.size), dtype=complex)
    for pump_index, pump_phase in enumerate(phases):
        for probe_index, probe_phase in enumerate(phases):
            phased = propagate_pump_probe_hierarchy(
                time,
                parameters=parameters,
                drive_frequency=0.0,
                pump_envelope=lambda t, phase=pump_phase: np.exp(-1j * phase) * pump(t),
                probe_envelope=lambda t, phase=probe_phase: np.exp(-1j * phase) * probe(t),
                rtol=2e-9,
                atol=2e-11,
            )
            phased_alpha_21[pump_index, probe_index] = phased.cavity_orders[2, 1]

    def extract(v_pump: int, v_probe: int) -> np.ndarray:
        phase_weight = np.exp(
            1j
            * (
                v_pump * phases[:, None]
                + v_probe * phases[None, :]
            )
        )
        return np.mean(phase_weight[..., None] * phased_alpha_21, axis=(0, 1))

    cycled_01 = extract(0, 1)
    cycled_2m1 = extract(2, -1)
    phase_cycle_error_01 = float(np.max(np.abs(cycled_01 - explicit.alpha_3_01)))
    phase_cycle_error_2m1 = float(np.max(np.abs(cycled_2m1 - explicit.alpha_3_2m1)))

    figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), constrained_layout=True)
    axes[0].plot(time, np.real(explicit.alpha_3_01), color="black", lw=1.4, label="explicit")
    axes[0].plot(time, np.real(cycled_01), color="tab:red", ls="--", lw=1.0, label="phase cycling")
    axes[0].set(xlabel=r"Time $t\kappa$", ylabel=r"Re $\alpha^{(2)(1)}_{(0,1)}$", title="(0,1)")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(time, np.real(explicit.alpha_3_2m1), color="black", lw=1.4, label="explicit")
    axes[1].plot(time, np.real(cycled_2m1), color="tab:blue", ls="--", lw=1.0, label="phase cycling")
    axes[1].set(
        xlabel=r"Time $t\kappa$",
        ylabel=r"Re $\alpha^{(2)(1)}_{(2,-1)}$",
        title="(2,-1)",
    )
    axes[1].legend(frameon=False, fontsize=8)
    figure.savefig(figure_dir / "phase_component_validation.png", dpi=240)
    figure.savefig(figure_dir / "phase_component_validation.pdf")
    plt.close(figure)

    np.savez_compressed(
        data_dir / "phase_component_validation.npz",
        time=time,
        phases=phases,
        explicit_alpha_01=explicit.alpha_3_01,
        explicit_alpha_2m1=explicit.alpha_3_2m1,
        phase_cycled_alpha_01=cycled_01,
        phase_cycled_alpha_2m1=cycled_2m1,
        unresolved_alpha_21=unresolved.cavity_orders[2, 1],
    )
    metrics = {
        "explicit_components_sum_vs_unresolved_max_error": explicit_sum_error,
        "phase_cycle_grid_per_pulse": phase_count,
        "phase_cycling_01_max_error": phase_cycle_error_01,
        "phase_cycling_2m1_max_error": phase_cycle_error_2m1,
    }
    metadata = {
        "model": "explicit SM S.39-S.43 2LS phase hierarchy",
        "phase_convention": "input envelope multiplied by exp(-i Phi); extraction uses exp(+i v.Phi)",
        "metrics": metrics,
    }
    (data_dir / "phase_component_validation.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

