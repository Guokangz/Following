#!/usr/bin/env python3
"""Cross-check linear Eq. (34) against first-order semiclassical response."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.linear.input_output import spectra_from_susceptibility
from src.linear.susceptibility import two_level_susceptibility
from src.nonlinear.mean_field import (
    TwoLevelMeanFieldParameters,
    linear_cw_intracavity_response,
    propagate_mean_field,
)


def main() -> None:
    figure_dir = PROJECT_ROOT / "figures"
    data_dir = PROJECT_ROOT / "data"
    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    omega_c = omega_0 = 0.0
    kappa = 0.1
    collective_coupling = 2.0
    gamma_population = 0.1
    gamma_phi = 0.2
    gamma_linear = gamma_population + gamma_phi
    omega = np.linspace(-5.0, 5.0, 20_001)

    # Method A: linear paper Eqs. (34)/(38).
    chi = two_level_susceptibility(
        omega,
        omega_exc=omega_0,
        gamma=gamma_linear,
        collective_coupling=collective_coupling,
    )
    linear_spectra = spectra_from_susceptibility(
        omega,
        chi,
        omega_c=omega_c,
        kappa_left=kappa / 2.0,
        kappa_right=kappa / 2.0,
    )

    # Method B: PRL mean-field equations linearized about |g><g|.
    parameters = TwoLevelMeanFieldParameters(
        omega_c=omega_c,
        omega_0=omega_0,
        kappa=kappa,
        collective_coupling=collective_coupling,
        gamma=gamma_population,
        gamma_phi=gamma_phi,
    )
    semiclassical_response = linear_cw_intracavity_response(
        omega,
        parameters=parameters,
    )
    semiclassical_green = 1j * semiclassical_response
    semiclassical_transmission = (kappa / 2.0) ** 2 * np.abs(
        semiclassical_response
    ) ** 2

    complex_green_error = np.max(np.abs(semiclassical_green - linear_spectra.green))
    transmission_error = np.max(
        np.abs(semiclassical_transmission - linear_spectra.transmission)
    )

    # Independent full nonlinear ODE check in the infinitesimal-CW limit.
    cw_check_frequencies = np.array([-2.0, 0.0, 2.0])
    cw_drive_amplitude = 1e-6
    cw_time = np.linspace(0.0, 800.0, 4_001)
    cw_full_response = np.empty(cw_check_frequencies.size, dtype=complex)
    cw_max_trace_error = 0.0
    cw_min_density_eigenvalue = 1.0
    for index, frequency in enumerate(cw_check_frequencies):
        trajectory = propagate_mean_field(
            cw_time,
            parameters=parameters,
            drive_frequency=float(frequency),
            drive=lambda _, amplitude=cw_drive_amplitude: amplitude,
            rtol=1e-12,
            atol=1e-15,
        )
        cw_full_response[index] = trajectory.cavity_amplitude[-1] / cw_drive_amplitude
        cw_max_trace_error = max(
            cw_max_trace_error,
            float(np.max(np.abs(trajectory.trace - 1.0))),
        )
        cw_min_density_eigenvalue = min(
            cw_min_density_eigenvalue,
            float(np.min(trajectory.minimum_density_eigenvalue)),
        )
    cw_linear_response = linear_cw_intracavity_response(
        cw_check_frequencies,
        parameters=parameters,
    )
    cw_full_vs_linear_error = float(np.max(np.abs(cw_full_response - cw_linear_response)))

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(6.3, 5.0),
        constrained_layout=True,
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    axes[0].plot(
        omega,
        linear_spectra.transmission,
        color="black",
        lw=1.8,
        label="Linear Eq. (34)/(38)",
    )
    axes[0].plot(
        omega,
        semiclassical_transmission,
        color="tab:red",
        ls="--",
        lw=1.2,
        label="Semiclassical first order",
    )
    axes[0].set(ylabel="Transmission", xlim=(-5, 5), ylim=(-0.002, None))
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(
        omega,
        semiclassical_transmission - linear_spectra.transmission,
        color="tab:blue",
        lw=1.0,
    )
    axes[1].set(
        xlabel=r"Drive frequency $(\omega-\omega_c)$",
        ylabel="Difference",
        xlim=(-5, 5),
    )
    figure.savefig(figure_dir / "linear_vs_semiclassical_crosscheck.png", dpi=240)
    figure.savefig(figure_dir / "linear_vs_semiclassical_crosscheck.pdf")
    plt.close(figure)

    np.savez_compressed(
        data_dir / "linear_vs_semiclassical_crosscheck.npz",
        omega=omega,
        linear_green=linear_spectra.green,
        semiclassical_green=semiclassical_green,
        linear_transmission=linear_spectra.transmission,
        semiclassical_transmission=semiclassical_transmission,
        cw_check_frequencies=cw_check_frequencies,
        cw_full_response=cw_full_response,
        cw_linear_response=cw_linear_response,
    )
    metrics = {
        "max_complex_green_difference": float(complex_green_error),
        "max_transmission_difference": float(transmission_error),
        "full_nonlinear_weak_cw_vs_linear_response_max_difference": cw_full_vs_linear_error,
        "weak_cw_max_trace_error": cw_max_trace_error,
        "weak_cw_min_density_eigenvalue": cw_min_density_eigenvalue,
    }
    metadata = {
        "methods": {
            "A": "linear paper Eqs. (34)/(38)",
            "B": "PRL mean-field equations at first order / infinitesimal CW",
        },
        "parameters": {
            "omega_c": omega_c,
            "omega_0": omega_0,
            "kappa": kappa,
            "sqrtN_g": collective_coupling,
            "mean_field_gamma_population": gamma_population,
            "mean_field_gamma_phi": gamma_phi,
            "linear_paper_gamma": gamma_linear,
        },
        "convention_identity": "gamma_linear = gamma_population + gamma_phi",
        "phase_identity": "D^R = i (alpha/input) because the PRL cavity equation uses -input",
        "metrics": metrics,
    }
    (data_dir / "linear_vs_semiclassical_crosscheck.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
