#!/usr/bin/env python3
"""Validate the two-pulse hierarchy and the SM Eq. (S.20) FFT convention."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.nonlinear.fourier import positive_frequency_fft
from src.nonlinear.mean_field import TwoLevelMeanFieldParameters, propagate_mean_field
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
        gamma=0.2,
        gamma_phi=0.1,
    )
    dt = 0.002
    time = np.arange(0.0, 14.0, dt)
    pump_center = 2.0
    delay = 1.0
    probe_center = pump_center + delay
    pulse_width = 0.2
    pump = lambda t: np.exp(-0.5 * ((t - pump_center) / pulse_width) ** 2)
    probe = lambda t: np.exp(-0.5 * ((t - probe_center) / pulse_width) ** 2)

    hierarchy = propagate_pump_probe_hierarchy(
        time,
        parameters=parameters,
        drive_frequency=0.0,
        pump_envelope=pump,
        probe_envelope=probe,
        max_pump_order=2,
        max_probe_order=1,
        rtol=1e-10,
        atol=1e-12,
    )

    frequency, alpha_01_frequency = positive_frequency_fft(
        time,
        hierarchy.cavity_orders[0, 1],
        zero_padding_factor=4,
    )
    _, alpha_21_frequency = positive_frequency_fft(
        time,
        hierarchy.cavity_orders[2, 1],
        zero_padding_factor=4,
    )
    _, probe_frequency_numeric = positive_frequency_fft(
        time,
        np.array([probe(value) for value in time]),
        zero_padding_factor=4,
    )
    probe_frequency_analytic = pulse_width * np.exp(
        -0.5 * (frequency * pulse_width) ** 2
    ) * np.exp(1j * frequency * probe_center)
    gaussian_mask = np.abs(frequency) < 15.0
    gaussian_fft_error = float(
        np.max(
            np.abs(
                probe_frequency_numeric[gaussian_mask]
                - probe_frequency_analytic[gaussian_mask]
            )
        )
    )

    # Independent full-solution benchmark: subtract pump-only dynamics so
    # pump-only alpha^(3,0) terms cancel, then compare the probe-linear piece.
    eta_pump = 0.05
    eta_probe = 0.001
    full_both = propagate_mean_field(
        time,
        parameters=parameters,
        drive_frequency=0.0,
        drive=lambda t: eta_pump * pump(t) + eta_probe * probe(t),
        rtol=1e-10,
        atol=1e-12,
    )
    full_pump = propagate_mean_field(
        time,
        parameters=parameters,
        drive_frequency=0.0,
        drive=lambda t: eta_pump * pump(t),
        rtol=1e-10,
        atol=1e-12,
    )
    full_probe_induced = full_both.cavity_amplitude - full_pump.cavity_amplitude
    hierarchy_probe_induced = eta_probe * (
        hierarchy.cavity_orders[0, 1]
        + eta_pump**2 * hierarchy.cavity_orders[2, 1]
    )
    full_vs_hierarchy_error = float(
        np.max(np.abs(full_probe_induced - hierarchy_probe_induced))
    )

    spectral_mask = np.abs(frequency) <= 8.0
    third_order_interference = 2.0 * np.real(
        np.conjugate(alpha_01_frequency) * alpha_21_frequency
    )

    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.55), constrained_layout=True)
    axes[0].plot(time, np.abs(hierarchy.cavity_orders[0, 1]), label=r"$|\alpha^{(0)(1)}|$")
    axes[0].plot(time, np.abs(hierarchy.cavity_orders[2, 1]), label=r"$|\alpha^{(2)(1)}|$")
    axes[0].axvline(pump_center, color="tab:red", ls="--", lw=0.9, label="pump")
    axes[0].axvline(probe_center, color="tab:blue", ls=":", lw=0.9, label="probe")
    axes[0].set(xlabel=r"Time $t\kappa$", ylabel="Cavity coefficient")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(
        frequency[spectral_mask],
        third_order_interference[spectral_mask],
        color="tab:purple",
        lw=1.2,
    )
    axes[1].axvline(-3.0, color="tab:green", ls="--", lw=0.8)
    axes[1].axvline(3.0, color="tab:green", ls="--", lw=0.8)
    axes[1].set(
        xlabel=r"Rotating-frame frequency $\omega/\kappa$",
        ylabel=r"$2\operatorname{Re}[\alpha_{01}^*\alpha_{21}]$",
        xlim=(-8, 8),
    )
    figure.savefig(figure_dir / "pump_probe_hierarchy_validation.png", dpi=240)
    figure.savefig(figure_dir / "pump_probe_hierarchy_validation.pdf")
    plt.close(figure)

    np.savez_compressed(
        data_dir / "pump_probe_hierarchy_validation.npz",
        time=time,
        frequency=frequency,
        alpha_01_time=hierarchy.cavity_orders[0, 1],
        alpha_21_time=hierarchy.cavity_orders[2, 1],
        alpha_01_frequency=alpha_01_frequency,
        alpha_21_frequency=alpha_21_frequency,
        probe_frequency_numeric=probe_frequency_numeric,
        probe_frequency_analytic=probe_frequency_analytic,
        full_probe_induced=full_probe_induced,
        hierarchy_probe_induced=hierarchy_probe_induced,
    )
    metrics = {
        "gaussian_fft_max_error_abs_omega_below_15": gaussian_fft_error,
        "full_two_pulse_probe_induced_vs_hierarchy_max_error": full_vs_hierarchy_error,
        "max_abs_alpha_01": float(np.max(np.abs(hierarchy.cavity_orders[0, 1]))),
        "max_abs_alpha_21": float(np.max(np.abs(hierarchy.cavity_orders[2, 1]))),
        "max_abs_forbidden_alpha_11": float(np.max(np.abs(hierarchy.cavity_orders[1, 1]))),
    }
    metadata = {
        "model": "2LS pump-probe hierarchy through (2)(1)",
        "parameters": {
            "kappa": parameters.kappa,
            "sqrtN_g": parameters.collective_coupling,
            "gamma": parameters.gamma,
            "gamma_phi": parameters.gamma_phi,
            "pump_center": pump_center,
            "probe_center": probe_center,
            "pulse_width_standard_deviation": pulse_width,
        },
        "fourier_convention": "1/sqrt(2pi) integral exp(+i omega t) f(t) dt",
        "metrics": metrics,
    }
    (data_dir / "pump_probe_hierarchy_validation.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

