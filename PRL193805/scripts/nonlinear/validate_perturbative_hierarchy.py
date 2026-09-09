#!/usr/bin/env python3
"""Validate the genuine n=0...3 perturbative hierarchy against full dynamics."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.nonlinear.mean_field import TwoLevelMeanFieldParameters, propagate_mean_field
from src.nonlinear.perturbative import propagate_perturbative_hierarchy


def _fit_power(amplitude: np.ndarray, error: np.ndarray) -> float:
    valid = np.isfinite(error) & (error > 0)
    return float(np.polyfit(np.log(amplitude[valid]), np.log(error[valid]), 1)[0])


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
    time = np.linspace(0.0, 10.0, 4_001)
    pulse_center = 2.0
    pulse_width = 0.2
    envelope = lambda t: np.exp(-0.5 * ((t - pulse_center) / pulse_width) ** 2)

    hierarchy = propagate_perturbative_hierarchy(
        time,
        parameters=parameters,
        drive_frequency=0.0,
        envelope=envelope,
        max_order=3,
        rtol=1e-11,
        atol=1e-13,
    )

    eta_values = np.array([0.04, 0.06, 0.09, 0.135, 0.20])
    cavity_error = np.empty_like(eta_values)
    coherence_error = np.empty_like(eta_values)
    population_error = np.empty_like(eta_values)
    total_error = np.empty_like(eta_values)
    full_max_population = np.empty_like(eta_values)
    for index, eta in enumerate(eta_values):
        full = propagate_mean_field(
            time,
            parameters=parameters,
            drive_frequency=0.0,
            drive=lambda t, amplitude=float(eta): amplitude * envelope(t),
            rtol=1e-11,
            atol=1e-13,
        )
        reconstruction = hierarchy.reconstruct(float(eta), through_order=3)
        cavity_error[index] = np.max(
            np.abs(full.cavity_amplitude - reconstruction.cavity_amplitude)
        )
        coherence_error[index] = np.max(
            np.abs(full.molecular_coherence - reconstruction.molecular_coherence)
        )
        population_error[index] = np.max(
            np.abs(full.excited_population - reconstruction.excited_population)
        )
        total_error[index] = max(
            cavity_error[index],
            coherence_error[index],
            population_error[index],
        )
        full_max_population[index] = np.max(full.excited_population)

    cavity_power = _fit_power(eta_values, cavity_error)
    population_power = _fit_power(eta_values, population_error)
    total_power = _fit_power(eta_values, total_error)

    figure, axes = plt.subplots(1, 2, figsize=(9.1, 3.55), constrained_layout=True)
    axes[0].plot(time, np.abs(hierarchy.cavity_orders[1]), label=r"$|\alpha^{(1)}|$")
    axes[0].plot(time, np.abs(hierarchy.cavity_orders[3]), label=r"$|\alpha^{(3)}|$")
    axes[0].plot(time, hierarchy.excited_population_orders[2], label=r"$p_e^{(2)}$")
    axes[0].set(xlabel=r"Time $t\kappa$", ylabel="Perturbative coefficient")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].loglog(eta_values, cavity_error, "o-", label=rf"cavity (slope {cavity_power:.2f})")
    axes[1].loglog(
        eta_values,
        population_error,
        "s-",
        label=rf"population (slope {population_power:.2f})",
    )
    axes[1].loglog(eta_values, total_error, "^-", label=rf"total (slope {total_power:.2f})")
    guide = total_error[0] * (eta_values / eta_values[0]) ** 4
    axes[1].loglog(eta_values, guide, "k--", lw=1.0, label=r"$\eta^4$ guide")
    axes[1].set(xlabel=r"Input amplitude $\eta$", ylabel="Max reconstruction error")
    axes[1].legend(frameon=False, fontsize=7)
    figure.savefig(figure_dir / "perturbative_hierarchy_validation.png", dpi=240)
    figure.savefig(figure_dir / "perturbative_hierarchy_validation.pdf")
    plt.close(figure)

    np.savez_compressed(
        data_dir / "perturbative_hierarchy_validation.npz",
        time=time,
        eta_values=eta_values,
        cavity_orders=hierarchy.cavity_orders,
        coherence_orders=hierarchy.coherence_orders,
        excited_population_orders=hierarchy.excited_population_orders,
        cavity_error=cavity_error,
        coherence_error=coherence_error,
        population_error=population_error,
        total_error=total_error,
        full_max_population=full_max_population,
    )
    metrics = {
        "eta_values": eta_values.tolist(),
        "cavity_max_errors": cavity_error.tolist(),
        "coherence_max_errors": coherence_error.tolist(),
        "population_max_errors": population_error.tolist(),
        "total_max_errors": total_error.tolist(),
        "fitted_cavity_error_power": cavity_power,
        "fitted_population_error_power": population_power,
        "fitted_total_error_power": total_power,
        "max_abs_alpha_order_2": float(np.max(np.abs(hierarchy.cavity_orders[2]))),
        "max_abs_population_order_1": float(
            np.max(np.abs(hierarchy.excited_population_orders[1]))
        ),
        "max_abs_population_order_3": float(
            np.max(np.abs(hierarchy.excited_population_orders[3]))
        ),
    }
    metadata = {
        "model": "2LS rotating-frame semiclassical hierarchy, direct order propagation",
        "max_order": 3,
        "parameters": {
            "kappa": parameters.kappa,
            "sqrtN_g": parameters.collective_coupling,
            "gamma": parameters.gamma,
            "gamma_phi": parameters.gamma_phi,
            "pulse_center": pulse_center,
            "pulse_width_standard_deviation": pulse_width,
        },
        "metrics": metrics,
    }
    (data_dir / "perturbative_hierarchy_validation.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

