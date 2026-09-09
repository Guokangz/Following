#!/usr/bin/env python3
"""Reproduce linear-paper Fig. 6: arbitrary stationary 3LS populations."""

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
from src.linear.susceptibility import three_level_susceptibility


def main() -> None:
    output_dir = PROJECT_ROOT / "figures" / "linear"
    data_dir = PROJECT_ROOT / "data" / "linear"
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    detuning = np.linspace(-5.0, 5.0, 20_001)
    omega_12 = 1.0  # Inferred arbitrary-unit scale; only ratios are stated.
    omega_c = omega_12
    omega = omega_c + detuning
    level_energies = np.array([0.0, omega_12, 3.0 * omega_12])
    kappa = 0.1
    gamma = 0.3
    collective_coupling = 1.0
    population_cases = np.array(
        [
            [0.7, 0.2, 0.1],
            [0.48, 0.48, 0.04],
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        ]
    )

    susceptibility = np.empty((3, omega.size), dtype=complex)
    transmission = np.empty((3, omega.size))
    reflection = np.empty((3, omega.size))
    absorption = np.empty((3, omega.size))
    flux_error = np.empty(3)
    absorption_peak_counts = []
    for index, populations in enumerate(population_cases):
        susceptibility[index] = three_level_susceptibility(
            omega,
            level_energies=level_energies,
            populations=populations,
            gamma=gamma,
            collective_coupling=collective_coupling,
        )
        spectra = spectra_from_susceptibility(
            omega,
            susceptibility[index],
            omega_c=omega_c,
            kappa_left=kappa / 2.0,
            kappa_right=kappa / 2.0,
        )
        transmission[index] = spectra.transmission
        reflection[index] = spectra.reflection
        absorption[index] = spectra.absorption
        flux_error[index] = np.max(np.abs(spectra.flux_sum - 1.0))
        if np.max(absorption[index]) > 1e-12:
            peaks, _ = find_peaks(
                absorption[index],
                prominence=0.015 * np.max(absorption[index]),
            )
            absorption_peak_counts.append(int(peaks.size))
        else:
            absorption_peak_counts.append(0)

    empty_spectra = spectra_from_susceptibility(
        omega,
        np.zeros_like(omega, dtype=complex),
        omega_c=omega_c,
        kappa_left=kappa / 2.0,
        kappa_right=kappa / 2.0,
    )
    saturated_empty_cavity_error = float(
        max(
            np.max(np.abs(transmission[2] - empty_spectra.transmission)),
            np.max(np.abs(reflection[2] - empty_spectra.reflection)),
            np.max(np.abs(absorption[2] - empty_spectra.absorption)),
        )
    )

    figure, axes = plt.subplots(1, 3, figsize=(10.8, 3.25), constrained_layout=True, sharey=True)
    for index, (ax, label) in enumerate(zip(axes, ("a", "b", "c"))):
        ax.plot(detuning, transmission[index], color="red", lw=1.3, label="T")
        ax.plot(detuning, absorption[index], color="black", lw=1.3, label="A")
        ax.plot(detuning, reflection[index], color="blue", lw=1.3, label="R")
        ax.set(xlim=(-5, 5), ylim=(-0.02, 1.05), xlabel=r"Energy $(\omega-\omega_{ph})$")
        ax.legend(frameon=True, fontsize=8)
        ax.text(-0.16, 1.02, label, transform=ax.transAxes, fontsize=14, fontweight="bold")
        populations = population_cases[index]
        ax.text(
            0.04,
            0.93,
            rf"$p=({populations[0]:.2g},{populations[1]:.2g},{populations[2]:.2g})$",
            transform=ax.transAxes,
            va="top",
            fontsize=8,
        )
    axes[0].set_ylabel("Intensity")
    figure.savefig(output_dir / "fig6_reproduction.png", dpi=240)
    figure.savefig(output_dir / "fig6_reproduction.pdf")
    plt.close(figure)

    chi_figure, chi_axes = plt.subplots(1, 3, figsize=(10.8, 3.15), constrained_layout=True, sharey=True)
    for index, (ax, label) in enumerate(zip(chi_axes, ("a", "b", "c"))):
        ax.plot(detuning, np.imag(susceptibility[index]), color="tab:green", lw=1.4)
        ax.axhline(0.0, color="0.65", lw=0.7)
        ax.set(xlim=(-5, 5), xlabel=r"Energy $(\omega-\omega_{ph})$")
        ax.text(-0.16, 1.02, label, transform=ax.transAxes, fontsize=14, fontweight="bold")
    chi_axes[0].set_ylabel(r"$\operatorname{Im}\chi(\omega)$")
    chi_figure.savefig(output_dir / "fig6_susceptibility_sanity.png", dpi=240)
    chi_figure.savefig(output_dir / "fig6_susceptibility_sanity.pdf")
    plt.close(chi_figure)

    np.savez_compressed(
        data_dir / "fig6_reproduction.npz",
        detuning=detuning,
        omega=omega,
        level_energies=level_energies,
        populations=population_cases,
        susceptibility=susceptibility,
        transmission=transmission,
        reflection=reflection,
        absorption=absorption,
    )
    metrics = {
        "max_flux_errors": flux_error.tolist(),
        "max_abs_chi_each_case": np.max(np.abs(susceptibility), axis=1).tolist(),
        "absorption_peak_counts": absorption_peak_counts,
        "fully_saturated_vs_empty_cavity_max_error": saturated_empty_cavity_error,
    }
    metadata = {
        "paper": "Yuen-Zhou and Koner, JCP 160, 154107 (2024)",
        "figure": "6(a-c)",
        "equations": [34, 45, 46],
        "parameters": {
            "omega_12": omega_12,
            "omega_12_status": "inferred arbitrary-unit choice",
            "omega_ph": omega_c,
            "omega_23": 2.0 * omega_12,
            "kappa": kappa,
            "gamma": gamma,
            "sqrtN_g_each_transition": collective_coupling,
            "relative_transition_dipoles": "all one; inferred",
            "populations": population_cases.tolist(),
        },
        "positive_frequency_RWA_transitions": ["1-2", "1-3", "2-3"],
        "metrics": metrics,
    }
    (data_dir / "fig6_reproduction.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

