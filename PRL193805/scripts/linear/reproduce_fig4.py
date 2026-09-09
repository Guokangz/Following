#!/usr/bin/env python3
"""Reproduce linear-paper Fig. 4: energetic disorder and line narrowing."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, simpson
from scipy.signal import find_peaks

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.linear.input_output import spectra_from_susceptibility
from src.linear.susceptibility import (
    gaussian_density,
    gaussian_disorder_susceptibility,
    inhomogeneous_susceptibility_grid,
    lorentzian_density,
    lorentzian_disorder_susceptibility,
)


def _spectra(omega: np.ndarray, chi: np.ndarray, kappa: float):
    return spectra_from_susceptibility(
        omega,
        chi,
        omega_c=0.0,
        kappa_left=kappa / 2.0,
        kappa_right=kappa / 2.0,
    )


def _absolute_half_max_width(
    omega: np.ndarray,
    signal: np.ndarray,
    peak_index: int,
) -> float:
    target = 0.5 * signal[peak_index]
    left = peak_index
    while left > 0 and signal[left] >= target:
        left -= 1
    right = peak_index
    while right < signal.size - 1 and signal[right] >= target:
        right += 1
    if left == 0 or right == signal.size - 1:
        return float("nan")
    left_crossing = np.interp(target, signal[left : left + 2], omega[left : left + 2])
    right_crossing = np.interp(
        target,
        signal[right - 1 : right + 1][::-1],
        omega[right - 1 : right + 1][::-1],
    )
    return float(right_crossing - left_crossing)


def _two_peak_metrics(omega: np.ndarray, signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    peaks, _ = find_peaks(signal, prominence=0.02 * np.max(signal))
    if peaks.size < 2:
        raise RuntimeError("did not resolve two polariton absorption peaks")
    selected = np.sort(peaks[np.argsort(signal[peaks])[-2:]])
    widths = np.array([_absolute_half_max_width(omega, signal, peak) for peak in selected])
    return omega[selected], widths


def _lorentzian_quad_chi(
    omega: float,
    *,
    sigma: float,
    gamma: float,
    collective_coupling: float,
) -> complex:
    def integrand(excitation: float) -> complex:
        density = float(lorentzian_density(excitation, mean=0.0, sigma=sigma))
        return -(collective_coupling**2) * density / (
            omega - excitation + 0.5j * gamma
        )

    real = quad(lambda x: float(np.real(integrand(x))), -np.inf, np.inf, epsabs=2e-10)[0]
    imag = quad(lambda x: float(np.imag(integrand(x))), -np.inf, np.inf, epsabs=2e-10)[0]
    return real + 1j * imag


def _save_panel(
    output_dir: Path,
    name: str,
    omega: np.ndarray,
    curves: list[tuple[np.ndarray, str, str, str]],
    *,
    ylabel: str,
    ylim: tuple[float, float],
) -> None:
    figure, ax = plt.subplots(figsize=(4.2, 3.35), constrained_layout=True)
    for values, color, linestyle, label in curves:
        ax.plot(omega, values, color=color, ls=linestyle, lw=1.45, label=label)
    ax.set(xlim=(-5, 5), ylim=ylim, xlabel=r"Energy $(\omega-\omega_{ph})$", ylabel=ylabel)
    ax.legend(frameon=True, fontsize=8)
    figure.savefig(output_dir / f"{name}.png", dpi=240)
    figure.savefig(output_dir / f"{name}.pdf")
    plt.close(figure)


def main() -> None:
    output_dir = PROJECT_ROOT / "figures" / "linear"
    data_dir = PROJECT_ROOT / "data" / "linear"
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    omega = np.linspace(-5.0, 5.0, 20_001)
    sigma = 1.0
    kappa = gamma = 0.1
    collective_coupling = 1.5

    gaussian_chi = gaussian_disorder_susceptibility(
        omega,
        mean=0.0,
        sigma=sigma,
        gamma=gamma,
        collective_coupling=collective_coupling,
    )
    lorentzian_chi = lorentzian_disorder_susceptibility(
        omega,
        mean=0.0,
        sigma=sigma,
        gamma=gamma,
        collective_coupling=collective_coupling,
    )
    gaussian_spectra = _spectra(omega, gaussian_chi, kappa)
    lorentzian_spectra = _spectra(omega, lorentzian_chi, kappa)

    # Deterministic-grid convergence against the analytic Voigt expression.
    validation_omega = np.linspace(-5.0, 5.0, 401)
    gaussian_grid_sizes = np.array([2_001, 4_001, 8_001, 16_001])
    gaussian_grid_errors = []
    analytic_validation = gaussian_disorder_susceptibility(
        validation_omega,
        mean=0.0,
        sigma=sigma,
        gamma=gamma,
        collective_coupling=collective_coupling,
    )
    for grid_size in gaussian_grid_sizes:
        excitation = np.linspace(-8.0 * sigma, 8.0 * sigma, int(grid_size))
        density = gaussian_density(excitation, mean=0.0, sigma=sigma)
        numeric = inhomogeneous_susceptibility_grid(
            validation_omega,
            excitation,
            density,
            gamma=gamma,
            collective_coupling=collective_coupling,
        )
        gaussian_grid_errors.append(float(np.max(np.abs(numeric - analytic_validation))))

    gaussian_normalization = float(
        simpson(
            gaussian_density(np.linspace(-8, 8, 16_001), mean=0.0, sigma=sigma),
            x=np.linspace(-8, 8, 16_001),
        )
    )
    lorentzian_normalization = float(
        quad(lambda x: float(lorentzian_density(x, mean=0.0, sigma=sigma)), -np.inf, np.inf)[0]
    )
    lorentzian_check_frequencies = np.array([-3.0, 0.0, 3.0])
    lorentzian_quad = np.array(
        [
            _lorentzian_quad_chi(
                value,
                sigma=sigma,
                gamma=gamma,
                collective_coupling=collective_coupling,
            )
            for value in lorentzian_check_frequencies
        ]
    )
    lorentzian_exact_check = lorentzian_disorder_susceptibility(
        lorentzian_check_frequencies,
        mean=0.0,
        sigma=sigma,
        gamma=gamma,
        collective_coupling=collective_coupling,
    )
    lorentzian_quad_error = float(np.max(np.abs(lorentzian_quad - lorentzian_exact_check)))

    coupling_values = np.array([1.3, 3.0])
    narrowed_absorption = np.empty((coupling_values.size, omega.size))
    peak_positions = np.empty((coupling_values.size, 2))
    peak_fwhm = np.empty((coupling_values.size, 2))
    for index, coupling in enumerate(coupling_values):
        chi = gaussian_disorder_susceptibility(
            omega,
            mean=0.0,
            sigma=sigma,
            gamma=gamma,
            collective_coupling=float(coupling),
        )
        absorption = _spectra(omega, chi, kappa).absorption
        narrowed_absorption[index] = absorption / np.max(absorption)
        peak_positions[index], peak_fwhm[index] = _two_peak_metrics(
            omega,
            narrowed_absorption[index],
        )

    bare_absorption = np.imag(
        gaussian_disorder_susceptibility(
            omega,
            mean=0.0,
            sigma=sigma,
            gamma=gamma,
            collective_coupling=1.0,
        )
    )
    bare_absorption /= np.max(bare_absorption)
    bare_peak = int(np.argmax(bare_absorption))
    bare_fwhm = _absolute_half_max_width(omega, bare_absorption, bare_peak)

    panel_a_curves = [
        (gaussian_spectra.transmission, "red", "-", "T"),
        (gaussian_spectra.reflection, "blue", "-", "R"),
        (gaussian_spectra.absorption, "black", "-", "A"),
    ]
    panel_b_curves = [
        (lorentzian_spectra.transmission, "red", "-", "T"),
        (lorentzian_spectra.reflection, "blue", "-", "R"),
        (lorentzian_spectra.absorption, "black", "-", "A"),
    ]
    panel_c_curves = [
        (narrowed_absorption[0], "blue", "-", r"$g\sqrt{N}=1.3$"),
        (narrowed_absorption[1], "red", "-", r"$g\sqrt{N}=3.0$"),
        (bare_absorption, "black", ":", r"bare $\Im\chi$"),
    ]
    _save_panel(output_dir, "fig4a_reproduction", omega, panel_a_curves, ylabel="Intensity", ylim=(-0.02, 1.05))
    _save_panel(output_dir, "fig4b_reproduction", omega, panel_b_curves, ylabel="Intensity", ylim=(-0.02, 1.05))
    _save_panel(output_dir, "fig4c_reproduction", omega, panel_c_curves, ylabel="Normalized absorption", ylim=(-0.02, 1.52))

    figure, axes = plt.subplots(1, 3, figsize=(10.8, 3.25), constrained_layout=True)
    for panel_label, ax, curves, ylabel, ylim in zip(
        ("a", "b", "c"),
        axes,
        (panel_a_curves, panel_b_curves, panel_c_curves),
        ("Intensity", "Intensity", "Normalized absorption"),
        ((-0.02, 1.05), (-0.02, 1.05), (-0.02, 1.52)),
    ):
        for values, color, linestyle, label in curves:
            ax.plot(omega, values, color=color, ls=linestyle, lw=1.3, label=label)
        ax.set(xlim=(-5, 5), ylim=ylim, xlabel=r"Energy $(\omega-\omega_{ph})$", ylabel=ylabel)
        ax.legend(frameon=True, fontsize=7)
        ax.text(-0.18, 1.02, panel_label, transform=ax.transAxes, fontsize=14, fontweight="bold")
    figure.savefig(output_dir / "fig4_reproduction.png", dpi=240)
    figure.savefig(output_dir / "fig4_reproduction.pdf")
    plt.close(figure)

    np.savez_compressed(
        data_dir / "fig4_reproduction.npz",
        omega=omega,
        gaussian_susceptibility=gaussian_chi,
        gaussian_transmission=gaussian_spectra.transmission,
        gaussian_reflection=gaussian_spectra.reflection,
        gaussian_absorption=gaussian_spectra.absorption,
        lorentzian_susceptibility=lorentzian_chi,
        lorentzian_transmission=lorentzian_spectra.transmission,
        lorentzian_reflection=lorentzian_spectra.reflection,
        lorentzian_absorption=lorentzian_spectra.absorption,
        coupling_values=coupling_values,
        normalized_polariton_absorption=narrowed_absorption,
        bare_normalized_absorption=bare_absorption,
        peak_positions=peak_positions,
        peak_fwhm=peak_fwhm,
        bare_fwhm=bare_fwhm,
        gaussian_grid_sizes=gaussian_grid_sizes,
        gaussian_grid_errors=np.asarray(gaussian_grid_errors),
        lorentzian_check_frequencies=lorentzian_check_frequencies,
        lorentzian_quad_values=lorentzian_quad,
    )

    metrics = {
        "gaussian_distribution_integral_on_plus_minus_8_sigma": gaussian_normalization,
        "lorentzian_distribution_integral": lorentzian_normalization,
        "gaussian_grid_sizes": gaussian_grid_sizes.tolist(),
        "gaussian_grid_max_errors": gaussian_grid_errors,
        "lorentzian_quad_max_error": lorentzian_quad_error,
        "max_flux_error_gaussian": float(np.max(np.abs(gaussian_spectra.flux_sum - 1.0))),
        "max_flux_error_lorentzian": float(np.max(np.abs(lorentzian_spectra.flux_sum - 1.0))),
        "polariton_peak_positions": peak_positions.tolist(),
        "polariton_peak_fwhm": peak_fwhm.tolist(),
        "bare_molecular_fwhm": bare_fwhm,
    }
    metadata = {
        "paper": "Yuen-Zhou and Koner, JCP 160, 154107 (2024)",
        "figure": "4(a-c)",
        "equations": [34, 41, 42, 43],
        "parameters": {
            "sigma": sigma,
            "omega_ph_minus_mean_omega_exc": 0.0,
            "kappa": kappa,
            "gamma": gamma,
            "sqrtN_g_fig4ab": collective_coupling,
            "sqrtN_g_fig4c": coupling_values.tolist(),
            "kappa_left": kappa / 2,
            "kappa_right": kappa / 2,
        },
        "eq43_choice": "normalized Cauchy numerator sigma/2; printed (sigma/2)^2 is non-normalized",
        "metrics": metrics,
    }
    (data_dir / "fig4_reproduction.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

