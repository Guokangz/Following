#!/usr/bin/env python3
"""Reproduce PRL Fig. 3 from the direct pump-probe hierarchy."""

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
from src.nonlinear.mean_field import TwoLevelMeanFieldParameters
from src.nonlinear.perturbative import propagate_perturbative_hierarchy
from src.nonlinear.pump_probe import propagate_pump_probe_hierarchy


def _gaussian(time: float | np.ndarray, center: float, width: float):
    return np.exp(-0.5 * ((np.asarray(time) - center) / width) ** 2)


def _interpolate_real(
    source_frequency: np.ndarray,
    source_values: np.ndarray,
    target_frequency: np.ndarray,
) -> np.ndarray:
    return np.interp(target_frequency, source_frequency, np.real(source_values))


def _dt_line(
    *,
    delay: float,
    gamma_phi: float,
    dt: float,
    target_frequency: np.ndarray,
    pump_center: float,
    pulse_width: float,
    end_time: float,
    zero_padding_factor: int,
) -> np.ndarray:
    parameters = TwoLevelMeanFieldParameters(
        omega_c=0.0,
        omega_0=0.0,
        kappa=1.0,
        collective_coupling=3.0,
        gamma=0.0,
        gamma_phi=gamma_phi,
    )
    time = np.arange(0.0, end_time, dt)
    probe_center = pump_center + delay
    hierarchy = propagate_pump_probe_hierarchy(
        time,
        parameters=parameters,
        drive_frequency=0.0,
        pump_envelope=lambda value: _gaussian(value, pump_center, pulse_width),
        probe_envelope=lambda value: _gaussian(value, probe_center, pulse_width),
        max_pump_order=2,
        max_probe_order=1,
        rtol=2e-9,
        atol=2e-11,
    )
    frequency, alpha_01 = positive_frequency_fft(
        time,
        hierarchy.cavity_orders[0, 1],
        zero_padding_factor=zero_padding_factor,
    )
    _, alpha_21 = positive_frequency_fft(
        time,
        hierarchy.cavity_orders[2, 1],
        zero_padding_factor=zero_padding_factor,
    )
    # Eq. (8) is only needed on the plotted finite-frequency window.  Avoid
    # dividing by the underflowed Gaussian spectrum near the FFT Nyquist edge.
    useful = np.abs(frequency) <= max(25.0, 2.0 * np.max(np.abs(target_frequency)))
    useful_frequency = frequency[useful]
    probe_spectrum_abs_sq = pulse_width**2 * np.exp(
        -(useful_frequency * pulse_width) ** 2
    )
    dt_native = (
        0.25
        * 2.0
        * np.real(np.conjugate(alpha_01[useful]) * alpha_21[useful])
        / probe_spectrum_abs_sq
    )
    return _interpolate_real(useful_frequency, dt_native, target_frequency)


def _pump_dynamics(
    *,
    gamma_phi: float,
    time: np.ndarray,
    pump_center: float,
    pulse_width: float,
) -> tuple[np.ndarray, np.ndarray]:
    parameters = TwoLevelMeanFieldParameters(
        omega_c=0.0,
        omega_0=0.0,
        kappa=1.0,
        collective_coupling=3.0,
        gamma=0.0,
        gamma_phi=gamma_phi,
    )
    hierarchy = propagate_perturbative_hierarchy(
        time,
        parameters=parameters,
        drive_frequency=0.0,
        envelope=lambda value: _gaussian(value, pump_center, pulse_width),
        max_order=2,
        rtol=2e-10,
        atol=2e-12,
    )
    population = hierarchy.excited_population_orders[2]
    polarization_abs_sq = np.abs(hierarchy.coherence_orders[1]) ** 2
    return population, polarization_abs_sq


def main() -> None:
    figure_dir = PROJECT_ROOT / "figures" / "nonlinear"
    data_dir = PROJECT_ROOT / "data" / "nonlinear"
    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    pump_center = 2.0
    pulse_width = 0.05
    delays = np.linspace(0.0, 6.0, 61)
    target_frequency = np.linspace(-5.0, 5.0, 501)
    dt = 0.005
    end_time = 20.0
    zero_padding_factor = 4
    dephasing_values = np.array([0.0, 1.0])

    heatmaps = np.empty((2, delays.size, target_frequency.size))
    for dephasing_index, gamma_phi in enumerate(dephasing_values):
        for delay_index, delay in enumerate(delays):
            heatmaps[dephasing_index, delay_index] = _dt_line(
                delay=float(delay),
                gamma_phi=float(gamma_phi),
                dt=dt,
                target_frequency=target_frequency,
                pump_center=pump_center,
                pulse_width=pulse_width,
                end_time=end_time,
                zero_padding_factor=zero_padding_factor,
            )
        print(f"completed gamma_phi/kappa={gamma_phi:g}")

    # Timestep convergence at representative delays, using the raw Eq. (8)
    # coefficient without changing plotting limits.
    convergence_delays = np.array([0.0, 2.0, 4.0, 6.0])
    timestep_errors = np.empty((2, convergence_delays.size))
    for dephasing_index, gamma_phi in enumerate(dephasing_values):
        for delay_index, delay in enumerate(convergence_delays):
            coarse = _dt_line(
                delay=float(delay),
                gamma_phi=float(gamma_phi),
                dt=dt,
                target_frequency=target_frequency,
                pump_center=pump_center,
                pulse_width=pulse_width,
                end_time=end_time,
                zero_padding_factor=zero_padding_factor,
            )
            fine = _dt_line(
                delay=float(delay),
                gamma_phi=float(gamma_phi),
                dt=dt / 2.0,
                target_frequency=target_frequency,
                pump_center=pump_center,
                pulse_width=pulse_width,
                end_time=end_time,
                zero_padding_factor=zero_padding_factor,
            )
            timestep_errors[dephasing_index, delay_index] = np.max(np.abs(coarse - fine))

    pump_time = np.arange(0.0, 10.0, 0.0025)
    pump_population = np.empty((2, pump_time.size))
    pump_polarization_abs_sq = np.empty((2, pump_time.size))
    for index, gamma_phi in enumerate(dephasing_values):
        pump_population[index], pump_polarization_abs_sq[index] = _pump_dynamics(
            gamma_phi=float(gamma_phi),
            time=pump_time,
            pump_center=pump_center,
            pulse_width=pulse_width,
        )

    polariton_mask = (
        (np.abs(target_frequency - 3.0) < 0.5)
        | (np.abs(target_frequency + 3.0) < 0.5)
    )
    delay_envelope = np.sqrt(np.mean(heatmaps[..., polariton_mask] ** 2, axis=-1))
    cross_section_index = int(np.argmin(np.abs(delays - 4.0)))
    cross_section = heatmaps[1, cross_section_index]

    common_limit = float(np.max(np.abs(heatmaps)))
    figure, axes = plt.subplots(2, 2, figsize=(9.0, 6.6), constrained_layout=True)
    heatmap_axes = axes[0]
    for index, (ax, label, title) in enumerate(
        zip(
            heatmap_axes,
            ("a", "b"),
            (r"no dephasing: $\gamma_\phi=0$", r"with dephasing: $\gamma_\phi=\kappa$"),
        )
    ):
        mesh = ax.pcolormesh(
            target_frequency,
            delays,
            heatmaps[index],
            shading="auto",
            cmap="seismic",
            vmin=-common_limit,
            vmax=common_limit,
        )
        ax.axvline(-3.0, color="limegreen", ls="--", lw=1.0)
        ax.axvline(3.0, color="limegreen", ls="--", lw=1.0)
        ax.set(
            xlim=(-5, 5),
            ylim=(0, 6),
            xlabel=r"Rotating-frame frequency $\omega/\kappa$",
            ylabel=r"Pulse delay $\tau_\Delta\kappa$",
            title=title,
        )
        ax.text(-0.14, 1.03, label, transform=ax.transAxes, fontsize=14, fontweight="bold")
        figure.colorbar(mesh, ax=ax, label=r"$\Delta T^{(3)}(\omega)$")
    if np.max(np.abs(cross_section)) > 0:
        overlay = 4.0 + 0.38 * cross_section / np.max(np.abs(cross_section))
        heatmap_axes[1].plot(target_frequency, overlay, color="black", lw=1.0)

    for index, (ax, label, gamma_phi) in enumerate(
        zip(axes[1], ("c", "d"), dephasing_values)
    ):
        population_scale = np.max(np.abs(pump_population[index]))
        polarization_scale = np.max(pump_polarization_abs_sq[index])
        ax.plot(
            pump_time,
            _gaussian(pump_time, pump_center, pulse_width),
            color="lightcoral",
            lw=1.0,
            label=r"$f_p(t)$",
        )
        for shown_delay in (1.0, 2.0):
            ax.plot(
                pump_time,
                _gaussian(pump_time, pump_center + shown_delay, pulse_width),
                color="lightcoral",
                ls="--",
                lw=0.8,
            )
        ax.plot(
            pump_time,
            pump_population[index] / population_scale,
            color="tab:green",
            lw=1.2,
            label="Population",
        )
        ax.plot(
            pump_time,
            pump_polarization_abs_sq[index] / polarization_scale,
            color="tab:blue",
            ls="--",
            lw=1.2,
            label=r"$|P(t)|^2$",
        )
        ax.set(
            xlim=(0, 10),
            ylim=(-0.03, 1.08),
            xlabel=r"Time $t\kappa$",
            ylabel="Normalized pump response",
            title=rf"$\gamma_\phi/\kappa={gamma_phi:g}$",
        )
        ax.text(-0.14, 1.03, label, transform=ax.transAxes, fontsize=14, fontweight="bold")
        ax.legend(frameon=False, fontsize=7)
    figure.savefig(figure_dir / "fig3_reproduction.png", dpi=240)
    figure.savefig(figure_dir / "fig3_reproduction.pdf")
    plt.close(figure)

    cross_figure, cross_ax = plt.subplots(figsize=(5.0, 3.4), constrained_layout=True)
    cross_ax.plot(target_frequency, cross_section, color="black", lw=1.4)
    cross_ax.axvline(-3.0, color="limegreen", ls="--", lw=0.9)
    cross_ax.axvline(3.0, color="limegreen", ls="--", lw=0.9)
    cross_ax.set(
        xlabel=r"Rotating-frame frequency $\omega/\kappa$",
        ylabel=r"$\Delta T^{(3)}(\omega)$",
        title=r"$\gamma_\phi=\kappa$, $\tau_\Delta=4\kappa^{-1}$",
        xlim=(-5, 5),
    )
    cross_figure.savefig(figure_dir / "fig3b_delay4_cross_section.png", dpi=240)
    cross_figure.savefig(figure_dir / "fig3b_delay4_cross_section.pdf")
    plt.close(cross_figure)

    np.savez_compressed(
        data_dir / "fig3_reproduction.npz",
        frequency=target_frequency,
        delays=delays,
        gamma_phi_values=dephasing_values,
        differential_transmission=heatmaps,
        delay_envelope=delay_envelope,
        cross_section_delay4=cross_section,
        pump_time=pump_time,
        pump_population_order2=pump_population,
        pump_polarization_order1_abs_sq=pump_polarization_abs_sq,
        convergence_delays=convergence_delays,
        timestep_errors=timestep_errors,
    )
    metrics = {
        "raw_dt_min_each_dephasing": np.min(heatmaps, axis=(1, 2)).tolist(),
        "raw_dt_max_each_dephasing": np.max(heatmaps, axis=(1, 2)).tolist(),
        "timestep_convergence_delays": convergence_delays.tolist(),
        "dt_0p005_vs_0p0025_max_errors": timestep_errors.tolist(),
        "late_delay_envelope_at_6": delay_envelope[:, -1].tolist(),
        "late_pump_population_order2_at_t10": pump_population[:, -1].tolist(),
        "late_pump_polarization_abs_sq_at_t10": pump_polarization_abs_sq[:, -1].tolist(),
    }
    metadata = {
        "paper": "Reitz, Koner, and Yuen-Zhou, PRL 134, 193803 (2025)",
        "figure": "3(a-d)",
        "equations": [2, 3, 6, 8, "SM S.7-S.20"],
        "parameters": {
            "kappa": 1.0,
            "sqrtN_g": 3.0,
            "gamma": 0.0,
            "gamma_phi": dephasing_values.tolist(),
            "tau_w_standard_deviation": pulse_width,
            "omega_c_minus_omega_0": 0.0,
            "omega_p_minus_omega_0": 0.0,
            "omega_probe_minus_omega_0": 0.0,
            "pump_arrival": pump_center,
        },
        "numerical": {
            "time_step": dt,
            "end_time": end_time,
            "zero_padding_factor": zero_padding_factor,
            "input_amplitude_factor_eta_p_squared": 1.0,
        },
        "metrics": metrics,
    }
    (data_dir / "fig3_reproduction.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
