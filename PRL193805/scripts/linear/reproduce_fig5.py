#!/usr/bin/env python3
"""Reproduce linear-paper Fig. 5: vibronic coupling."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.linear.input_output import spectra_from_susceptibility
from src.linear.susceptibility import franck_condon_weights, vibronic_susceptibility


def main() -> None:
    output_dir = PROJECT_ROOT / "figures" / "linear"
    data_dir = PROJECT_ROOT / "data" / "linear"
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    omega = np.linspace(-5.0, 5.0, 20_001)
    omega_c = omega_exc = 0.0
    kappa = gamma = 0.1
    collective_coupling = 1.0
    omega_v = 0.3
    huang_rhys = 3.0
    m, weights = franck_condon_weights(huang_rhys, tolerance=1e-12)
    m_max = int(m[-1])

    chi = vibronic_susceptibility(
        omega,
        omega_exc=omega_exc,
        omega_v=omega_v,
        huang_rhys=huang_rhys,
        gamma=gamma,
        collective_coupling=collective_coupling,
        m_max=m_max,
    )
    spectra = spectra_from_susceptibility(
        omega,
        chi,
        omega_c=omega_c,
        kappa_left=kappa / 2.0,
        kappa_right=kappa / 2.0,
    )

    trial_m_max = np.unique(np.array([6, 9, 12, 15, 18, m_max], dtype=int))
    trial_m_max = trial_m_max[trial_m_max <= m_max]
    truncation_errors = []
    for trial in trial_m_max:
        trial_chi = vibronic_susceptibility(
            omega,
            omega_exc=omega_exc,
            omega_v=omega_v,
            huang_rhys=huang_rhys,
            gamma=gamma,
            collective_coupling=collective_coupling,
            m_max=int(trial),
        )
        trial_absorption = spectra_from_susceptibility(
            omega,
            trial_chi,
            omega_c=omega_c,
            kappa_left=kappa / 2.0,
            kappa_right=kappa / 2.0,
        ).absorption
        truncation_errors.append(float(np.max(np.abs(trial_absorption - spectra.absorption))))

    figure, ax = plt.subplots(figsize=(5.0, 3.55), constrained_layout=True)
    ax.plot(omega, spectra.transmission, color="red", lw=1.35, label="T")
    ax.plot(omega, spectra.absorption, color="black", lw=1.35, label="A")
    ax.plot(omega, spectra.reflection, color="blue", lw=1.35, label="R")
    ax.set(xlim=(-5, 5), ylim=(-0.02, 1.05), xlabel=r"Energy $(\omega-\omega_{ph})$", ylabel="Intensity")
    ax.legend(frameon=True, fontsize=8)
    figure.savefig(output_dir / "fig5_reproduction.png", dpi=240)
    figure.savefig(output_dir / "fig5_reproduction.pdf")
    plt.close(figure)

    weight_figure, weight_ax = plt.subplots(figsize=(4.7, 3.35), constrained_layout=True)
    weight_ax.stem(m, weights, linefmt="tab:purple", markerfmt="o", basefmt="black")
    weight_ax.set(
        xlabel="Vibronic quantum number m",
        ylabel="Franck-Condon weight",
        xlim=(-0.5, m_max + 0.5),
    )
    weight_ax.text(
        0.98,
        0.95,
        rf"$\sum_{{m=0}}^{{{m_max}}}w_m={np.sum(weights):.12f}$",
        transform=weight_ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    weight_figure.savefig(output_dir / "fig5_franck_condon_weights.png", dpi=240)
    weight_figure.savefig(output_dir / "fig5_franck_condon_weights.pdf")
    plt.close(weight_figure)

    np.savez_compressed(
        data_dir / "fig5_reproduction.npz",
        omega=omega,
        susceptibility=chi,
        transmission=spectra.transmission,
        reflection=spectra.reflection,
        absorption=spectra.absorption,
        flux_sum=spectra.flux_sum,
        vibronic_quantum_number=m,
        franck_condon_weights=weights,
        trial_m_max=trial_m_max,
        truncation_max_absorption_error=np.asarray(truncation_errors),
    )
    metrics = {
        "selected_m_max": m_max,
        "franck_condon_sum": float(np.sum(weights)),
        "omitted_weight_bound": float(1.0 - np.sum(weights)),
        "weighted_mean_m": float(np.sum(m * weights)),
        "trial_m_max": trial_m_max.tolist(),
        "truncation_max_absorption_error": truncation_errors,
        "max_flux_error": float(np.max(np.abs(spectra.flux_sum - 1.0))),
    }
    metadata = {
        "paper": "Yuen-Zhou and Koner, JCP 160, 154107 (2024)",
        "figure": "5",
        "equations": [34, 44],
        "parameters": {
            "omega_ph": omega_c,
            "omega_exc": omega_exc,
            "kappa": kappa,
            "gamma": gamma,
            "sqrtN_g": collective_coupling,
            "omega_v": omega_v,
            "S": huang_rhys,
        },
        "metrics": metrics,
    }
    (data_dir / "fig5_reproduction.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

