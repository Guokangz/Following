# Molecular polariton response reproduction

Physics-first numerical reproduction of:

1. Joel Yuen-Zhou and Arghadip Koner, "Linear response of molecular polaritons," J. Chem. Phys. 160, 154107 (2024).
2. Michael Reitz, Arghadip Koner, and Joel Yuen-Zhou, "Nonlinear Semiclassical Spectroscopy of Ultrafast Molecular Polariton Dynamics," Phys. Rev. Lett. 134, 193803 (2025).

The project computes spectra from the published Hamiltonians, response functions, mean-field equations, perturbative hierarchy, and input-output relations. It does not digitize or fit curves from paper images.

## Current status

| Paper / figure | Status |
|---|---|
| Linear Fig. 3 | REPRODUCED |
| Linear Fig. 4 | REPRODUCED |
| Linear Fig. 5 | REPRODUCED |
| Linear Fig. 6 | REPRODUCED with documented inferred Fig. 6 scale/dipoles |
| Nonlinear N0 mean field | IMPLEMENTED and cross-checked against linear response |
| Nonlinear N1 perturbative hierarchy | IMPLEMENTED through third order and scaling-validated |
| Nonlinear N2 pump-probe hierarchy | IMPLEMENTED through `(2)(1)` with analytic FFT benchmark |
| Nonlinear Fig. 3 | REPRODUCED |
| Nonlinear N4 phase hierarchy | IMPLEMENTED and independently phase-cycling validated |
| Nonlinear Fig. 4(a,b) | REPRODUCED; decay-rate convention documented |
| Nonlinear Fig. 1(b,c) | REPRODUCED as a physics-based vector schematic |
| Nonlinear Fig. 2 | REPRODUCED as a pathway-tree vector schematic |
| Nonlinear Fig. 4(c,d) | REPRODUCED as phase-labelled vector diagrams |

See:

- `notes/paper_map.md` for the equation/parameter/figure audit.
- `notes/linear_validation.md` for numerical convergence and quantitative results.
- `notes/linear_nonlinear_connection.md` for the mandatory weak-field interface proof.
- `notes/perturbative_validation.md` for order-by-order/full-solution error scaling.
- `notes/pump_probe_validation.md` for the two-pulse hierarchy and Fourier convention.
- `notes/nonlinear_fig3_validation.md` for the delay-frequency heatmap and pump dynamics.
- `notes/phase_component_validation.md` for explicit hierarchy/phase-cycling agreement.
- `notes/nonlinear_fig4_validation.md` for phase heatmaps and the DQC decay fit.
- `notes/conceptual_figures.md` for the reconstructed pathway schematics.
- `notes/final_reproduction_report.md` for figure-by-figure status and discrepancies.

## Run the linear reproduction

Requirements: Python 3, NumPy, SciPy, and Matplotlib.

```bash
python3 -m unittest discover -s tests -v
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/linear/reproduce_all_linear.py
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/crosscheck_linear_vs_semiclassical.py
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/nonlinear/reproduce_fig3.py
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/nonlinear/reproduce_fig4.py
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/nonlinear/reproduce_all_nonlinear.py
```

Individual workflows:

```bash
python3 scripts/linear/reproduce_fig3.py
python3 scripts/linear/reproduce_fig4.py
python3 scripts/linear/reproduce_fig5.py
python3 scripts/linear/reproduce_fig6.py
```

Each workflow writes PNG/PDF figures to `figures/linear/` and compressed raw arrays plus JSON metadata to `data/linear/`.

## Project layout

```text
ref/                 papers and supplemental material
notes/               paper audit and validation reports
src/linear/          susceptibilities and input-output physics
src/nonlinear/       mean-field, perturbative, pump-probe, and phase hierarchies
scripts/linear/      reproducible figure workflows
scripts/nonlinear/   nonlinear validation and figure workflows
figures/linear/      rendered PNG/PDF outputs
figures/nonlinear/   rendered PNG/PDF outputs
data/linear/         raw arrays and metadata
data/nonlinear/      raw arrays and metadata
tests/               numerical and limiting-case checks
tmp/pdfs/            temporary text extraction/rendering for paper audit
```

The nonlinear implementation was started only after its first-order weak-field
response had been cross-checked against the frozen linear baseline.
