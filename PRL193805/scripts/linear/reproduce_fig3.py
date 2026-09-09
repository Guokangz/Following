#!/usr/bin/env python3
"""Reproduce linear-paper Fig. 3 from Eqs. (34), (37), and (38)."""

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

from src.linear.input_output import spectra_from_susceptibility
from src.linear.susceptibility import (
    thermal_population_factor,
    two_level_poles,
    two_level_susceptibility,
)


def _peak_positions(omega: np.ndarray, signal: np.ndarray) -> np.ndarray:
    indices, _ = find_peaks(signal, prominence=0.01 * np.max(signal))
    if indices.size < 2:
        raise RuntimeError("fewer than two polariton peaks were resolved")
    strongest = indices[np.argsort(signal[indices])[-2:]]
    return np.sort(omega[strongest])


def main() -> None:
    output_dir = PROJECT_ROOT / "figures" / "linear"
    data_dir = PROJECT_ROOT / "data" / "linear"
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    omega_c = 0.0
    omega_exc = 0.0
    kappa = 0.1
    gamma = 0.3
    collective_coupling = 2.0
    kappa_left = kappa_right = kappa / 2.0
    omega = np.linspace(-5.0, 5.0, 20_001)

    chi_zero_t = two_level_susceptibility(
        omega,
        omega_exc=omega_exc,
        gamma=gamma,
        collective_coupling=collective_coupling,
    )
    fig3a = spectra_from_susceptibility(
        omega,
        chi_zero_t,
        omega_c=omega_c,
        kappa_left=kappa_left,
        kappa_right=kappa_right,
    )

    max_flux_error = float(np.max(np.abs(fig3a.flux_sum - 1.0)))
    if max_flux_error > 5e-12:
        raise RuntimeError(f"T+R+A conservation failed: {max_flux_error:.3e}")

    absorption_peaks = _peak_positions(omega, fig3a.absorption)
    poles = two_level_poles(
        omega_c=omega_c,
        omega_exc=omega_exc,
        kappa=kappa,
        gamma=gamma,
        collective_coupling=collective_coupling,
    )
    pole_positions = np.sort(poles.real)

    temperature_ratios = np.array([np.inf, 1.1, 0.1])
    population_factors = thermal_population_factor(temperature_ratios)
    thermal_absorption = np.empty((temperature_ratios.size, omega.size))
    for index, factor in enumerate(population_factors):
        chi = two_level_susceptibility(
            omega,
            omega_exc=omega_exc,
            gamma=gamma,
            collective_coupling=collective_coupling,
            population_factor=float(factor),
        )
        thermal_absorption[index] = spectra_from_susceptibility(
            omega,
            chi,
            omega_c=omega_c,
            kappa_left=kappa_left,
            kappa_right=kappa_right,
        ).absorption

    population_factor_grid = np.linspace(0.0, 1.0, 501)
    thermal_pole_splitting = np.empty_like(population_factor_grid)
    for index, factor in enumerate(population_factor_grid):
        factor_poles = two_level_poles(
            omega_c=omega_c,
            omega_exc=omega_exc,
            kappa=kappa,
            gamma=gamma,
            collective_coupling=collective_coupling,
            population_factor=float(factor),
        )
        thermal_pole_splitting[index] = np.ptp(factor_poles.real)

    figure, axes = plt.subplots(1, 2, figsize=(8.2, 3.45), constrained_layout=True)
    ax = axes[0]
    ax.plot(omega, fig3a.transmission, color="red", lw=1.5, label="T")
    ax.plot(omega, fig3a.absorption, color="black", lw=1.5, label="A")
    ax.plot(omega, fig3a.reflection, color="blue", lw=1.5, label="R")
    ax.set(xlim=(-5, 5), ylim=(-0.02, 1.05), xlabel=r"Energy $(\omega-\omega_{ph})$", ylabel="Intensity")
    ax.legend(frameon=True, fontsize=9)
    ax.text(-0.17, 1.02, "a", transform=ax.transAxes, fontsize=14, fontweight="bold")

    ax = axes[1]
    styles = [("black", "--"), ("0.35", ":"), ("0.75", "-")]
    labels = [r"$\infty$", "1.1", "0.1"]
    for signal, (color, linestyle), label in zip(thermal_absorption, styles, labels):
        ax.plot(omega, signal, color=color, ls=linestyle, lw=1.6, label=label)
    ax.set(xlim=(-5, 5), ylim=(-0.01, 0.42), xlabel=r"Energy $(\omega-\omega_{ph})$", ylabel="Absorption")
    ax.legend(title=r"$\Delta E/k_BT$", frameon=True, fontsize=9, title_fontsize=9)
    ax.text(-0.17, 1.02, "b", transform=ax.transAxes, fontsize=14, fontweight="bold")

    png_path = output_dir / "fig3_reproduction.png"
    pdf_path = output_dir / "fig3_reproduction.pdf"
    figure.savefig(png_path, dpi=240)
    figure.savefig(pdf_path)
    plt.close(figure)

    # Stable per-panel filenames requested by the reproduction specification.
    for panel_name, panel_index in (("fig3a_reproduction", 0), ("fig3b_reproduction", 1)):
        panel_figure, panel_ax = plt.subplots(figsize=(4.2, 3.45), constrained_layout=True)
        if panel_index == 0:
            panel_ax.plot(omega, fig3a.transmission, color="red", lw=1.5, label="T")
            panel_ax.plot(omega, fig3a.absorption, color="black", lw=1.5, label="A")
            panel_ax.plot(omega, fig3a.reflection, color="blue", lw=1.5, label="R")
            panel_ax.set(xlim=(-5, 5), ylim=(-0.02, 1.05), ylabel="Intensity")
            panel_ax.legend(frameon=True, fontsize=9)
        else:
            for signal, (color, linestyle), label in zip(thermal_absorption, styles, labels):
                panel_ax.plot(omega, signal, color=color, ls=linestyle, lw=1.6, label=label)
            panel_ax.set(xlim=(-5, 5), ylim=(-0.01, 0.42), ylabel="Absorption")
            panel_ax.legend(title=r"$\Delta E/k_BT$", frameon=True, fontsize=9, title_fontsize=9)
        panel_ax.set_xlabel(r"Energy $(\omega-\omega_{ph})$")
        panel_figure.savefig(output_dir / f"{panel_name}.png", dpi=240)
        panel_figure.savefig(output_dir / f"{panel_name}.pdf")
        plt.close(panel_figure)

    auxiliary_figure, auxiliary_ax = plt.subplots(figsize=(4.4, 3.5), constrained_layout=True)
    auxiliary_ax.plot(population_factor_grid, thermal_pole_splitting, color="tab:purple", lw=1.8)
    auxiliary_ax.scatter(
        population_factors,
        np.interp(population_factors, population_factor_grid, thermal_pole_splitting),
        color="black",
        s=24,
        zorder=3,
        label="Fig. 3(b) temperatures",
    )
    auxiliary_ax.set(
        xlabel=r"Thermal population factor $\tanh(\Delta E/2k_BT)$",
        ylabel="Real-part pole splitting",
        xlim=(0, 1),
        ylim=(0, 4.15),
    )
    auxiliary_ax.legend(frameon=False, fontsize=8)
    auxiliary_figure.savefig(output_dir / "fig3_thermal_splitting_sanity.png", dpi=240)
    auxiliary_figure.savefig(output_dir / "fig3_thermal_splitting_sanity.pdf")
    plt.close(auxiliary_figure)

    np.savez_compressed(
        data_dir / "fig3_reproduction.npz",
        omega=omega,
        susceptibility=fig3a.susceptibility,
        photon_green=fig3a.green,
        transmission=fig3a.transmission,
        reflection=fig3a.reflection,
        absorption=fig3a.absorption,
        flux_sum=fig3a.flux_sum,
        temperature_ratios=temperature_ratios,
        thermal_population_factors=population_factors,
        thermal_absorption=thermal_absorption,
        population_factor_grid=population_factor_grid,
        thermal_pole_splitting=thermal_pole_splitting,
        complex_poles=poles,
        absorption_peak_positions=absorption_peaks,
    )

    metadata = {
        "paper": "Yuen-Zhou and Koner, JCP 160, 154107 (2024)",
        "figure": "3(a,b)",
        "equations": [34, 37, 38],
        "frequency_units": "arbitrary",
        "frequency_axis": "omega - omega_ph",
        "parameters": {
            "omega_ph": omega_c,
            "omega_exc": omega_exc,
            "kappa": kappa,
            "gamma": gamma,
            "sqrtN_g_zero_temperature": collective_coupling,
            "kappa_left": kappa_left,
            "kappa_right": kappa_right,
            "balanced_cavity": True,
            "balanced_cavity_status": "inferred from paper figure convention",
            "temperature_ratios_DeltaE_over_kBT": ["infinity", 1.1, 0.1],
        },
        "conventions": {
            "cavity_denominator": "+ i kappa/2",
            "molecular_denominator": "+ i gamma/2",
            "thermal_G_squared_scaling": "tanh(DeltaE/(2 k_B T))",
        },
        "metrics": {
            "max_abs_T_plus_R_plus_A_minus_1": max_flux_error,
            "complex_poles": [[float(z.real), float(z.imag)] for z in poles],
            "absorption_peak_positions": absorption_peaks.tolist(),
            "absorption_peak_splitting": float(np.diff(absorption_peaks)[0]),
            "pole_real_part_splitting": float(np.diff(pole_positions)[0]),
            "left_right_peak_symmetry_error": float(abs(np.sum(absorption_peaks))),
        },
    }
    metadata_path = data_dir / "fig3_reproduction.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(json.dumps(metadata["metrics"], indent=2))
    print(f"wrote {png_path.relative_to(PROJECT_ROOT)}")
    print(f"wrote {pdf_path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
