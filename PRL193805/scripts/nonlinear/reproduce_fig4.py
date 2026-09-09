#!/usr/bin/env python3
"""Reproduce PRL Fig. 4(a,b): phase-resolved nonlinear response."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.nonlinear.fourier import positive_frequency_fft
from src.nonlinear.mean_field import TwoLevelMeanFieldParameters
from src.nonlinear.phase_components import propagate_phase_resolved_hierarchy


def _gaussian(time: float | np.ndarray, center: float, width: float):
    return np.exp(-0.5 * ((np.asarray(time) - center) / width) ** 2)


def _component_lines(
    *,
    delay: float,
    dt: float,
    target_frequency: np.ndarray,
    pump_center: float,
    pulse_width: float,
    end_time: float,
    zero_padding_factor: int,
) -> tuple[np.ndarray, np.ndarray]:
    parameters = TwoLevelMeanFieldParameters(
        omega_c=0.0,
        omega_0=0.0,
        kappa=1.0,
        collective_coupling=3.0,
        gamma=0.0,
        gamma_phi=0.1,
    )
    time = np.arange(0.0, end_time, dt)
    probe_center = pump_center + delay
    trajectory = propagate_phase_resolved_hierarchy(
        time,
        parameters=parameters,
        drive_frequency=0.0,
        pump_envelope=lambda value: _gaussian(value, pump_center, pulse_width),
        probe_envelope=lambda value: _gaussian(value, probe_center, pulse_width),
        rtol=2e-9,
        atol=2e-11,
    )
    frequency, alpha_probe = positive_frequency_fft(
        time,
        trajectory.alpha_probe,
        zero_padding_factor=zero_padding_factor,
    )
    _, alpha_3_01 = positive_frequency_fft(
        time,
        trajectory.alpha_3_01,
        zero_padding_factor=zero_padding_factor,
    )
    _, alpha_3_2m1 = positive_frequency_fft(
        time,
        trajectory.alpha_3_2m1,
        zero_padding_factor=zero_padding_factor,
    )
    useful = np.abs(frequency) <= max(25.0, 2.0 * np.max(np.abs(target_frequency)))
    useful_frequency = frequency[useful]
    probe_abs_sq = pulse_width**2 * np.exp(
        -(useful_frequency * pulse_width) ** 2
    )
    prefactor = 0.25 * 2.0 / probe_abs_sq
    component_01_native = prefactor * np.real(
        np.conjugate(alpha_probe[useful]) * alpha_3_01[useful]
    )
    component_2m1_native = prefactor * np.real(
        np.conjugate(alpha_probe[useful]) * alpha_3_2m1[useful]
    )
    component_01 = np.interp(target_frequency, useful_frequency, component_01_native)
    component_2m1 = np.interp(target_frequency, useful_frequency, component_2m1_native)
    return component_01, component_2m1


def _exponential_peak_fit(
    delay: np.ndarray,
    envelope: np.ndarray,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Fit the oscillatory DQC signal's local-maximum envelope.

    A pointwise log fit is not meaningful because the UP/LP beating produces
    deep minima.  The peak separation is about half a Rabi period, so the
    minimum index separation below is fixed by the sampled delay grid rather
    than tuned to a target lifetime.
    """

    delay_step = float(delay[1] - delay[0])
    minimum_peak_separation = max(1, int(round(0.4 / delay_step)))
    peak_indices, _ = find_peaks(
        envelope,
        distance=minimum_peak_separation,
        prominence=np.max(envelope) * 1e-4,
    )
    peak_indices = peak_indices[
        (delay[peak_indices] >= 0.3) & (delay[peak_indices] <= 6.0)
    ]
    slope, intercept = np.polyfit(
        delay[peak_indices], np.log(envelope[peak_indices]), 1
    )
    fitted = np.exp(intercept + slope * delay)
    residual = np.log(envelope[peak_indices]) - (
        intercept + slope * delay[peak_indices]
    )
    total = np.log(envelope[peak_indices]) - np.mean(np.log(envelope[peak_indices]))
    r_squared = 1.0 - float(np.sum(residual**2) / np.sum(total**2))
    return float(-1.0 / slope), r_squared, fitted, peak_indices


def main() -> None:
    figure_dir = PROJECT_ROOT / "figures" / "nonlinear"
    data_dir = PROJECT_ROOT / "data" / "nonlinear"
    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    pump_center = 2.0
    pulse_width = 0.1
    delays = np.linspace(0.0, 8.0, 81)
    target_frequency = np.linspace(-5.0, 5.0, 501)
    dt = 0.005
    end_time = 24.0
    zero_padding_factor = 4

    component_01 = np.empty((delays.size, target_frequency.size))
    component_2m1 = np.empty_like(component_01)
    for delay_index, delay in enumerate(delays):
        component_01[delay_index], component_2m1[delay_index] = _component_lines(
            delay=float(delay),
            dt=dt,
            target_frequency=target_frequency,
            pump_center=pump_center,
            pulse_width=pulse_width,
            end_time=end_time,
            zero_padding_factor=zero_padding_factor,
        )
    print("completed phase-resolved delay scan")

    convergence_delays = np.array([0.0, 1.0, 2.0, 4.0, 8.0])
    timestep_errors = np.empty((2, convergence_delays.size))
    for index, delay in enumerate(convergence_delays):
        coarse = _component_lines(
            delay=float(delay),
            dt=dt,
            target_frequency=target_frequency,
            pump_center=pump_center,
            pulse_width=pulse_width,
            end_time=end_time,
            zero_padding_factor=zero_padding_factor,
        )
        fine = _component_lines(
            delay=float(delay),
            dt=dt / 2.0,
            target_frequency=target_frequency,
            pump_center=pump_center,
            pulse_width=pulse_width,
            end_time=end_time,
            zero_padding_factor=zero_padding_factor,
        )
        timestep_errors[0, index] = np.max(np.abs(coarse[0] - fine[0]))
        timestep_errors[1, index] = np.max(np.abs(coarse[1] - fine[1]))

    polariton_mask = (
        (np.abs(target_frequency - 3.0) < 0.7)
        | (np.abs(target_frequency + 3.0) < 0.7)
    )
    envelope_01 = np.sqrt(np.mean(component_01[:, polariton_mask] ** 2, axis=1))
    envelope_2m1 = np.sqrt(np.mean(component_2m1[:, polariton_mask] ** 2, axis=1))
    fitted_decay_time, fit_r_squared, fitted_envelope, peak_indices = (
        _exponential_peak_fit(
        delays,
        envelope_2m1,
        )
    )
    paper_order_of_magnitude_time = 1.0 / 1.1
    # With the amplitude decay rates appearing explicitly in SM (S.39)-(S.43),
    # each resonant polariton amplitude decays at (kappa+gamma_phi)/4.  The DQC
    # pathway contains two residual first-order amplitudes, hence its envelope
    # decays at (kappa+gamma_phi)/2.
    hierarchy_envelope_time = 2.0 / 1.1

    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.75), constrained_layout=True)
    for ax, data, label, title in zip(
        axes,
        (component_01, component_2m1),
        ("a", "b"),
        ("(0,1): bright + dark", "(2,-1): bright only"),
    ):
        limit = float(np.max(np.abs(data)))
        mesh = ax.pcolormesh(
            target_frequency,
            delays,
            data,
            shading="auto",
            cmap="seismic",
            vmin=-limit,
            vmax=limit,
        )
        ax.axvline(-3.0, color="limegreen", ls="--", lw=1.0)
        ax.axvline(3.0, color="limegreen", ls="--", lw=1.0)
        ax.set(
            xlim=(-5, 5),
            ylim=(0, 8),
            xlabel=r"Rotating-frame frequency $\omega/\kappa$",
            ylabel=r"Pulse delay $\tau_\Delta\kappa$",
            title=title,
        )
        ax.text(-0.14, 1.03, label, transform=ax.transAxes, fontsize=14, fontweight="bold")
        figure.colorbar(mesh, ax=ax, label=r"Phase-resolved $\Delta T^{(3)}$")
    figure.savefig(figure_dir / "fig4_phase_components_reproduction.png", dpi=240)
    figure.savefig(figure_dir / "fig4_phase_components_reproduction.pdf")
    plt.close(figure)

    decay_figure, decay_ax = plt.subplots(figsize=(5.2, 3.5), constrained_layout=True)
    decay_ax.semilogy(
        delays,
        envelope_2m1,
        "o",
        ms=2.5,
        color="0.7",
        label="oscillatory polariton-window RMS",
    )
    decay_ax.semilogy(
        delays[peak_indices],
        envelope_2m1[peak_indices],
        "o",
        ms=4,
        color="tab:blue",
        label="local maxima used in fit",
    )
    decay_ax.semilogy(
        delays,
        fitted_envelope,
        color="tab:red",
        lw=1.3,
        label=rf"fit $\tau={fitted_decay_time:.3f}\,\kappa^{{-1}}$",
    )
    hierarchy_curve = fitted_envelope[0] * np.exp(-delays / hierarchy_envelope_time)
    decay_ax.semilogy(
        delays,
        hierarchy_curve,
        color="black",
        ls="--",
        lw=1.0,
        label=rf"SM hierarchy: $2/(\kappa+\gamma_\phi)={hierarchy_envelope_time:.3f}$",
    )
    decay_ax.set(
        xlabel=r"Pulse delay $\tau_\Delta\kappa$",
        ylabel="DQC polariton-window RMS",
        xlim=(0, 8),
    )
    decay_ax.legend(frameon=False, fontsize=8)
    decay_figure.savefig(figure_dir / "fig4_dqc_decay_fit.png", dpi=240)
    decay_figure.savefig(figure_dir / "fig4_dqc_decay_fit.pdf")
    plt.close(decay_figure)

    np.savez_compressed(
        data_dir / "fig4_phase_components_reproduction.npz",
        frequency=target_frequency,
        delays=delays,
        component_01=component_01,
        component_2m1=component_2m1,
        envelope_01=envelope_01,
        envelope_2m1=envelope_2m1,
        fitted_envelope_2m1=fitted_envelope,
        decay_fit_peak_indices=peak_indices,
        convergence_delays=convergence_delays,
        timestep_errors=timestep_errors,
    )
    metrics = {
        "component_01_raw_range": [float(np.min(component_01)), float(np.max(component_01))],
        "component_2m1_raw_range": [
            float(np.min(component_2m1)),
            float(np.max(component_2m1)),
        ],
        "timestep_convergence_delays": convergence_delays.tolist(),
        "dt_0p005_vs_0p0025_max_errors_01_2m1": timestep_errors.tolist(),
        "fitted_dqc_decay_time_kappa_inverse": fitted_decay_time,
        "fit_log_space_r_squared": fit_r_squared,
        "paper_order_of_magnitude_kappa_plus_gamma_phi_inverse": (
            paper_order_of_magnitude_time
        ),
        "sm_hierarchy_dqc_envelope_time_2_over_kappa_plus_gamma_phi": (
            hierarchy_envelope_time
        ),
        "relative_difference_from_sm_hierarchy_envelope_time": float(
            (fitted_decay_time - hierarchy_envelope_time) / hierarchy_envelope_time
        ),
        "envelope_01_at_delay8": float(envelope_01[-1]),
        "envelope_2m1_at_delay8": float(envelope_2m1[-1]),
    }
    metadata = {
        "paper": "Reitz, Koner, and Yuen-Zhou, PRL 134, 193803 (2025)",
        "figure": "4(a,b)",
        "equations": [9, 10, "SM S.29-S.35", "SM S.39-S.43 with g23=0"],
        "parameters": {
            "kappa": 1.0,
            "sqrtN_g": 3.0,
            "gamma": 0.0,
            "gamma_phi": 0.1,
            "tau_w_standard_deviation": pulse_width,
            "pump_arrival": pump_center,
            "all_detunings": 0.0,
        },
        "numerical": {
            "time_step": dt,
            "end_time": end_time,
            "zero_padding_factor": zero_padding_factor,
            "decay_peak_fit_window": [0.3, 6.0],
            "decay_definition": (
                "local maxima of the polariton-window RMS; the pointwise RMS "
                "oscillates due to UP/LP beating"
            ),
        },
        "metrics": metrics,
    }
    (data_dir / "fig4_phase_components_reproduction.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
