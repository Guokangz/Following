#!/usr/bin/env python3
"""Draw research-useful vector versions of PRL Figs. 1(b,c), 2, and 4(c,d).

The drawings are reconstructed from the perturbative equations and pathway
logic.  They do not contain raster content copied from the article.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.nonlinear.mean_field import TwoLevelMeanFieldParameters
from src.nonlinear.phase_components import propagate_phase_resolved_hierarchy


def _box(ax, xy, width, height, text, *, edge="0.25", face="white", size=10):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.025,rounding_size=0.04",
        linewidth=1.2,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=size)
    return patch


def _arrow(ax, start, end, *, color="0.15", width=1.4, style="-|>", connection="arc3"):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=11,
        linewidth=width,
        color=color,
        connectionstyle=connection,
    )
    ax.add_patch(patch)
    return patch


def _gaussian(value, center, width):
    return np.exp(-0.5 * ((np.asarray(value) - center) / width) ** 2)


def reproduce_fig1_logic(figure_dir: Path, data_dir: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 3.65), constrained_layout=True)

    ax = axes[0]
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    _box(
        ax,
        (0.03, 0.57),
        0.94,
        0.35,
        "free space\n" r"$\dot{\rho}^{(n)}=-i\mathcal{L}_0\rho^{(n)}$"
        "\n"
        r"$\qquad-i\mathcal{L}_{\mu}[\mathcal{E}\rho^{(n-1)}+\mathcal{E}^{*}\rho^{(n-1)}]$",
        edge="0.35",
        face="#f7f7f7",
        size=10,
    )
    ax.text(
        0.5,
        0.54,
        "external laser is linear in the preceding matter order",
        ha="center",
        va="top",
        fontsize=8,
        color="0.35",
    )
    _box(
        ax,
        (0.03, 0.08),
        0.94,
        0.35,
        "cavity\n" r"$\dot{\rho}^{(n)}=-i\mathcal{L}_0\rho^{(n)}$"
        "\n"
        r"$\qquad-i\sum_{j=0}^{n}(\alpha^{(n-j)}+\alpha^{(n-j)*})\mathcal{L}_{\mu}\rho^{(j)}$",
        edge="#d95f02",
        face="#fff5eb",
        size=10,
    )
    ax.text(
        0.5,
        0.04,
        r"matter polarization $P^{(j)}$ feeds back into every $\alpha^{(j)}$",
        ha="center",
        va="top",
        fontsize=8,
        color="#a33d00",
    )
    ax.text(0.01, 0.97, "a", fontsize=14, fontweight="bold", va="top")
    ax.set_title("Free-space vs cavity perturbative coupling", fontsize=11)

    # Use the actual first-order coupled equations to draw the stored fields.
    ax = axes[1]
    parameters = TwoLevelMeanFieldParameters(
        omega_c=0,
        omega_0=0,
        kappa=1,
        collective_coupling=3,
        gamma=0,
        gamma_phi=0.1,
    )
    time = np.linspace(0, 8, 3201)
    pump_center, probe_center, width = 1.4, 3.0, 0.12
    trajectory = propagate_phase_resolved_hierarchy(
        time,
        parameters=parameters,
        drive_frequency=0,
        pump_envelope=lambda value: _gaussian(value, pump_center, width),
        probe_envelope=lambda value: _gaussian(value, probe_center, width),
        rtol=2e-10,
        atol=2e-12,
    )
    pump_field = np.real(trajectory.alpha_pump)
    probe_field = np.real(trajectory.alpha_probe)
    normalization = max(np.max(np.abs(pump_field)), np.max(np.abs(probe_field)))
    pulse_scale = 0.75
    ax.plot(time, pulse_scale * _gaussian(time, pump_center, width) + 1.8, color="orangered", lw=1.4)
    ax.plot(time, pulse_scale * _gaussian(time, probe_center, width) + 1.8, color="dodgerblue", lw=1.4)
    ax.plot(time, pump_field / normalization, color="orangered", lw=1.25, label=r"$\Re\,\alpha^{(1)(0)}$")
    ax.plot(time, probe_field / normalization, color="dodgerblue", lw=1.25, label=r"$\Re\,\alpha^{(0)(1)}$")
    ax.axhline(0, color="0.5", lw=0.7)
    ax.annotate("pump", (pump_center, 2.58), ha="center", color="orangered", fontsize=9)
    ax.annotate("probe", (probe_center, 2.58), ha="center", color="dodgerblue", fontsize=9)
    ax.annotate(
        "stored intracavity tails",
        xy=(4.0, 0.52),
        xytext=(5.15, 1.08),
        arrowprops=dict(arrowstyle="->", color="0.25"),
        ha="center",
        fontsize=8,
    )
    ax.annotate(
        "",
        xy=(probe_center, 2.36),
        xytext=(pump_center, 2.36),
        arrowprops=dict(arrowstyle="<->", color="0.15"),
    )
    ax.text((pump_center + probe_center) / 2, 2.42, r"$\tau_\Delta$", ha="center", fontsize=9)
    ax.set(
        xlim=(0, 8),
        ylim=(-1.15, 2.75),
        xlabel=r"time $t\kappa$",
        yticks=[],
        title="Pulse storage creates reordered interactions",
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.text(0.01, 0.97, "b", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    for suffix in ("png", "pdf"):
        figure.savefig(figure_dir / f"fig1_perturbative_logic.{suffix}", dpi=240 if suffix == "png" else None)
    plt.close(figure)
    np.savez_compressed(
        data_dir / "fig1_storage_fields.npz",
        time=time,
        pump_input=_gaussian(time, pump_center, width),
        probe_input=_gaussian(time, probe_center, width),
        alpha_pump=trajectory.alpha_pump,
        alpha_probe=trajectory.alpha_probe,
    )


def reproduce_fig2_tree(figure_dir: Path) -> None:
    figure, ax = plt.subplots(figsize=(8.8, 4.6), constrained_layout=True)
    ax.set(xlim=(0, 10), ylim=(0, 7.2))
    ax.axis("off")

    nodes = {
        "rn": (5.0, 6.55, r"$\rho^{(n)}$"),
        "r0": (1.0, 5.0, r"$\rho^{(0)}$"),
        "r1": (3.0, 5.0, r"$\rho^{(1)}$"),
        "r2": (5.0, 5.0, r"$\rho^{(2)}$"),
        "rnm1": (8.3, 5.0, r"$\rho^{(n-1)}$"),
        "r10": (2.45, 3.35, r"$\rho^{(0)}$"),
        "r20": (4.35, 3.35, r"$\rho^{(0)}$"),
        "r21": (5.65, 3.35, r"$\rho^{(1)}$"),
        "rnm10": (7.4, 3.35, r"$\rho^{(0)}$"),
        "rnm11": (8.3, 3.35, r"$\rho^{(1)}$"),
        "rnm12": (9.25, 3.35, r"$\rho^{(n-2)}$"),
        "r210": (5.65, 1.7, r"$\rho^{(0)}$"),
        "rnm110": (8.3, 1.7, r"$\rho^{(0)}$"),
        "rnm120": (9.25, 1.7, r"$\rho^{(n-3)}$"),
    }
    conventional_edges = {("rn", "rnm1"), ("rnm1", "rnm12"), ("rnm12", "rnm120")}
    edges = [
        ("rn", "r0", r"$\alpha^{(n)}$"),
        ("rn", "r1", r"$\alpha^{(n-1)}$"),
        ("rn", "r2", r"$\alpha^{(n-2)}$"),
        ("rn", "rnm1", r"$\alpha^{(1)}$"),
        ("r1", "r10", r"$\alpha^{(1)}$"),
        ("r2", "r20", r"$\alpha^{(2)}$"),
        ("r2", "r21", r"$\alpha^{(1)}$"),
        ("rnm1", "rnm10", r"$\alpha^{(n-1)}$"),
        ("rnm1", "rnm11", r"$\alpha^{(n-2)}$"),
        ("rnm1", "rnm12", r"$\alpha^{(1)}$"),
        ("r21", "r210", r"$\alpha^{(1)}$"),
        ("rnm11", "rnm110", r"$\alpha^{(1)}$"),
        ("rnm12", "rnm120", r"$\alpha^{(1)}$"),
    ]
    for parent, child, label in edges:
        x0, y0, _ = nodes[parent]
        x1, y1, _ = nodes[child]
        color = "#e66101" if (parent, child) in conventional_edges else "0.12"
        ax.plot([x0, x1], [y0 - 0.18, y1 + 0.18], color=color, lw=2.0)
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.12, label, color=color, fontsize=9, ha="center")
    for x, y, label in nodes.values():
        ax.text(x, y, label, ha="center", va="center", fontsize=12)
    ax.text(6.67, 5.18, r"$\cdots$", fontsize=18, ha="center")
    ax.text(8.78, 2.55, r"$\cdots$", fontsize=16, ha="center")
    ax.text(5.0, 0.65, "cavity feedback generates all lower-order partitions", ha="center", fontsize=10)
    ax.text(
        5.0,
        0.18,
        r"one orange free-space chain; $2^{n-1}-1$ additional cavity pathways",
        ha="center",
        fontsize=10,
        color="#a33d00",
    )
    ax.set_title("Liouville-space excitation-pathway tree from Eq. (5b)", fontsize=12)
    for suffix in ("png", "pdf"):
        figure.savefig(figure_dir / f"fig2_pathway_tree.{suffix}", dpi=240 if suffix == "png" else None)
    plt.close(figure)


def _draw_liouville_path(ax, x, states, phases, *, title, color, crossed=False):
    bottom, spacing = 0.65, 0.92
    left, right = x - 0.30, x + 0.30
    top = bottom + spacing * (len(states) - 1)
    ax.plot([left, left], [bottom - 0.25, top + 0.35], color="0.15", lw=1.5)
    ax.plot([right, right], [bottom - 0.25, top + 0.35], color="0.15", lw=1.5)
    for index, state in enumerate(states):
        y = bottom + spacing * index
        ax.text(x, y, state, ha="center", va="center", fontsize=9, bbox=dict(facecolor="white", edgecolor="none", pad=0.3))
    sides = ["left", "right", "left"]
    for index, (phase, side) in enumerate(zip(phases, sides)):
        y0 = bottom + spacing * index + 0.18
        y1 = bottom + spacing * (index + 1) - 0.18
        if side == "left":
            start, end = (left - 0.55, y0), (left - 0.02, y1)
            tx = left - 0.62
            ha = "right"
        else:
            start, end = (right + 0.55, y0), (right + 0.02, y1)
            tx = right + 0.62
            ha = "left"
        _arrow(ax, start, end, color=color, width=1.3)
        ax.text(tx, (y0 + y1) / 2, phase, color=color, ha=ha, va="center", fontsize=9)
    ax.annotate(
        "emission",
        xy=(left, top + 0.32),
        xytext=(left - 0.56, top + 0.67),
        arrowprops=dict(arrowstyle="->", linestyle=":", color="0.25"),
        fontsize=7,
        ha="right",
    )
    ax.text(x, 4.30, title, ha="center", fontsize=9, color="0.15")
    if crossed:
        ax.plot([x - 0.58, x + 0.58], [1.2, 3.7], color="crimson", lw=2.0, alpha=0.8)
        ax.plot([x - 0.58, x + 0.58], [3.7, 1.2], color="crimson", lw=2.0, alpha=0.8)


def reproduce_fig4_pathways(figure_dir: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.6), constrained_layout=True)
    for ax in axes:
        ax.set(xlim=(0, 5), ylim=(0, 4.75))
        ax.axis("off")

    ax = axes[0]
    _draw_liouville_path(
        ax,
        1.35,
        [r"$|g\rangle\langle g|$", r"$|e\rangle\langle g|$", r"$|g\rangle\langle g|$", r"$|e\rangle\langle g|$"],
        [r"$\Phi_p$", r"$-\Phi_p$", r"$\Phi_{p'}$"],
        title="ground-state population",
        color="#1b9e77",
    )
    _draw_liouville_path(
        ax,
        3.65,
        [r"$|g\rangle\langle g|$", r"$|g\rangle\langle e|$", r"$|e\rangle\langle e|$", r"$|e\rangle\langle g|$"],
        [r"$-\Phi_p$", r"$\Phi_p$", r"$\Phi_{p'}$"],
        title="excited-state population",
        color="#1b9e77",
    )
    ax.set_title(r"(0,1): $\pm\Phi_p\mp\Phi_p+\Phi_{p'}$  — bright + dark", fontsize=11)
    ax.text(0.02, 0.98, "a", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    ax = axes[1]
    _draw_liouville_path(
        ax,
        1.35,
        [r"$|g\rangle\langle g|$", r"$|e\rangle\langle g|$", r"$|g\rangle\langle g|$", r"$|e\rangle\langle g|$"],
        [r"$\Phi_p$", r"$-\Phi_{p'}$", r"$\Phi_p$"],
        title="stored probe acts on pump coherence",
        color="#7570b3",
    )
    _draw_liouville_path(
        ax,
        3.65,
        [r"$|g\rangle\langle g|$", r"$|g\rangle\langle e|$", r"$|e\rangle\langle e|$", r"$|e\rangle\langle g|$"],
        [r"$-\Phi_{p'}$", r"$\Phi_p$", r"$\Phi_p$"],
        title="probe first: absent from probe DT",
        color="#7570b3",
        crossed=True,
    )
    ax.set_title(r"(2,-1): $\Phi_p+\Phi_p-\Phi_{p'}$  — bright only", fontsize=11)
    ax.text(0.02, 0.98, "b", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    for suffix in ("png", "pdf"):
        figure.savefig(figure_dir / f"fig4_feynman_pathways.{suffix}", dpi=240 if suffix == "png" else None)
    plt.close(figure)


def main() -> None:
    figure_dir = PROJECT_ROOT / "figures" / "nonlinear"
    data_dir = PROJECT_ROOT / "data" / "nonlinear"
    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    reproduce_fig1_logic(figure_dir, data_dir)
    reproduce_fig2_tree(figure_dir)
    reproduce_fig4_pathways(figure_dir)
    metadata = {
        "paper": "Reitz, Koner, and Yuen-Zhou, PRL 134, 193803 (2025)",
        "figures": ["1(b,c)", "2", "4(c,d)"],
        "source_equations": ["PRL 4-6", "PRL 9-10", "SM S.29-S.35", "SM S.39-S.43"],
        "construction": "programmatic matplotlib vector drawings; no paper raster content",
        "fig1_storage_trace": "computed from the first-order coupled cavity-molecule equations with G=3 kappa",
    }
    (data_dir / "conceptual_figures_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
