# Final reproduction report

Last updated: 2026-08-23

## Scope and scientific baseline

All spectra in this project are computed from the published susceptibilities,
input-output equations, mean-field equations, and perturbative hierarchies.
No plotted data were digitized from article images.  The local nonlinear
Supplemental Material was used for the Liouville vectorization, two-pulse
expansion, Fourier convention, differential-transmission derivation, and
phase-resolved equations.

The mandatory interface between the two papers is passed: the semiclassical
first-order response gives the same complex photon Green function and
transmission as the linear-response implementation.  The gridwise difference
is zero at floating-point precision, and the full weak-drive nonlinear ODE
differs from the first-order steady solution by at most `7.51e-9`.

## Linear-response paper

### Fig. 3 — REPRODUCED

- Source: Eqs. (34), (37), (38).
- Parameters: resonant photon/2LS ensemble, `kappa=0.1`, `gamma=0.3`,
  `G=sqrt(N)g=2`; balanced mirrors inferred.  Temperature ratios are exactly
  `DeltaE/(kBT)=infinity, 1.1, 0.1`.
- Method: analytic 2LS susceptibility passed to the common photon Green
  function and input-output implementation.
- Script: `scripts/linear/reproduce_fig3.py`.
- Outputs: `figures/linear/fig3_reproduction.{png,pdf}` and
  `data/linear/fig3_reproduction.{npz,json}`.
- Agreement: absorption peaks `+/-1.997`, split `3.994`; pole split `3.99875`;
  exact sampled symmetry; maximum `T+R+A-1` error `2.22e-16`.  The thermal
  contraction follows `G_eff^2/G^2=tanh[DeltaE/(2kBT)]`.
- Remaining difference: only the inferred balanced-mirror convention.

### Fig. 4 — REPRODUCED

- Source: Eqs. (34), (41)-(43).
- Parameters: Gaussian/Lorentzian width `sigma=1`, `kappa=gamma=0.1`,
  `G=1.5`; panel (c) uses `G=1.3,3.0`.
- Method: normalized deterministic distributions, Simpson integration,
  analytic Voigt/Cauchy convolution benchmarks, and numerical FWHM.
- Script: `scripts/linear/reproduce_fig4.py`.
- Outputs: `figures/linear/fig4_reproduction.{png,pdf}` and associated raw
  arrays/metadata.
- Agreement: Gaussian 2001-point integration error `5.59e-9` and converged
  error below `5e-15`; Lorentzian error `1.78e-15`; flux error `2.22e-16`.
  Polariton FWHM contracts from `1.2183` at `G=1.3` to `0.1557` at `G=3`,
  versus bare molecular FWHM `2.4086`.
- Remaining difference: accepted-manuscript Eq. (43) has a non-normalized
  printed numerator.  The code uses the normalized Cauchy density and records
  this apparent typo.

### Fig. 5 — REPRODUCED

- Source: Eq. (44) and Poisson Franck-Condon factors.
- Parameters: resonance, `kappa=gamma=0.1`, `G=1`, `omega_v=0.3`, `S=3`.
- Method: explicit vibronic pole sum with convergence in `m_max`.
- Script: `scripts/linear/reproduce_fig5.py`.
- Outputs: `figures/linear/fig5_reproduction.{png,pdf}`, FC-weight sanity
  plot, and raw arrays/metadata.
- Agreement: selected `m_max=22`, omitted FC weight `2.07e-13`, flux error
  `2.22e-16`; the lower polariton and vibronic comb follow the article.
- Remaining difference: none beyond plotting layout.

### Fig. 6 — REPRODUCED

- Source: Eqs. (45), (46).
- Parameters: the three published populations, `omega_ph=omega_12`,
  `omega_23=2omega_12`, `kappa=0.1`, `gamma=0.3`, `G=1`.
- Method: direct sum of population-difference susceptibilities for all active
  upward transitions.
- Script: `scripts/linear/reproduce_fig6.py`.
- Outputs: `figures/linear/fig6_reproduction.{png,pdf}`, direct `Im chi`
  plot, and raw arrays/metadata.
- Agreement: four, three, and zero molecular-polariton peaks in cases (a),
  (b), and (c); equal populations give `chi=0` and exact empty-cavity recovery.
- Remaining difference: `omega_12=1` and equal relative transition dipoles
  are inferred because the article gives only ratios/a single coupling.

## Nonlinear semiclassical paper

### Fig. 1 — PARTIALLY REPRODUCED

- Source: Eqs. (4)-(6) and computed first-order cavity dynamics.
- Reproduced content: panels (b,c), the free-space/cavity coupling structure
  and pulse-storage mechanism.  The original decorative cavity rendering in
  panel (a) is intentionally not replicated.
- Script: `scripts/nonlinear/reproduce_conceptual_figures.py`.
- Output: `figures/nonlinear/fig1_perturbative_logic.{png,pdf}`.
- Quantitative basis: storage traces are propagated at `G=3 kappa`, not hand
  drawn.
- Remaining difference: panel (a) artwork is not required for the physical
  reproduction and remains omitted.

### Fig. 2 — REPRODUCED

- Source: Eq. (5b).
- Method: programmatic pathway tree for all lower-order light/matter
  partitions; orange branch identifies the free-space chain.
- Script/output: `scripts/nonlinear/reproduce_conceptual_figures.py` and
  `figures/nonlinear/fig2_pathway_tree.{png,pdf}`.
- Agreement: represents the paper's `2^(n-1)-1` additional cavity pathways.
- Remaining difference: vector layout is redesigned rather than pixel-matched.

### Fig. 3 — REPRODUCED

- Source: PRL Eqs. (2), (3), (6), (8); SM Eqs. (S.7)-(S.20).
- Parameters: `G=3 kappa`, `tau_w=0.05 kappa^-1`, resonance,
  `tau_p=2 kappa^-1`, `gamma=0`, and `gamma_phi=0,kappa`.
- Method: direct `(n)(m)` hierarchy through `alpha^(2)(1)`, positive-sign FFT,
  and Eq. (8); pump population/polarization propagated independently.
- Script: `scripts/nonlinear/reproduce_fig3.py`.
- Outputs: `figures/nonlinear/fig3_reproduction.{png,pdf}`, delay-4 cross
  section, and raw arrays/metadata.
- Agreement: raw ranges `[-0.01797,0.03377]` and
  `[-0.01793,0.02871]` naturally match the article scale.  Without dephasing,
  UP/LP coherence oscillates; with `gamma_phi=kappa`, late population is
  `7.48e-3` while `|P|^2=4.28e-6`, producing the contracted long-delay split.
  Timestep-halving errors are below `6.6e-7`.
- Remaining difference: exact delay sampling/window endpoints are absent from
  the paper and were fixed by convergence.

### Fig. 4 — REPRODUCED

- Source: PRL Eqs. (9), (10); SM Eqs. (S.29)-(S.35), (S.39)-(S.43).
- Parameters: `G=3 kappa`, `gamma_phi=0.1 kappa`, `gamma=0`, resonance,
  `tau_w=0.1 kappa^-1`.
- Method: explicit `(0,1)` and `(2,-1)` phase-index hierarchy, independently
  benchmarked by a `6x6` nonlinear phase cycle; Eq. (10) for the two heatmaps;
  programmatic double-sided diagrams.
- Scripts: `scripts/nonlinear/reproduce_fig4.py` and
  `scripts/nonlinear/reproduce_conceptual_figures.py`.
- Outputs: phase heatmaps, DQC decay fit, Feynman pathways, and raw arrays in
  `figures/nonlinear/` and `data/nonlinear/`.
- Agreement: explicit/phase-cycle field differences below `9.6e-11`.
  `(0,1)` persists through delay 8, whereas `(2,-1)` RMS falls to `2.79e-4`.
  Its local-maximum envelope gives `tau=1.812 kappa^-1`, `R^2=0.964`, agreeing
  within `0.32%` with the `2/(kappa+gamma_phi)=1.818` envelope implied by the
  SM amplitude equations.
- Remaining difference: the PRL text quotes the order-of-magnitude scale
  `~1/(kappa+gamma_phi)`.  The factor of two from its explicit amplitude-rate
  convention is documented rather than absorbed by changing linewidths.  The
  paper also omits the absolute perturbative pump amplitude in Eq. (10), so
  Fig. 4 absolute colorbar normalization is not fitted.

## Validation inventory

- 20 unit tests pass.
- Linear spectra independently enforce `T+R+A=1` and pass empty-cavity and
  zero-absorption limits.
- The perturbative reconstruction error scales as `eta^5.002` for the cavity
  field and `eta^4.001` for population/total state, as required after retaining
  third order.
- The FFT matches the analytic SM Gaussian transform to `4.94e-16`.
- The two-pulse hierarchy matches a full nonlinear pump-on minus pump-off
  trajectory to `3.88e-9` in the weak-field validation.

## Reproduction commands

```bash
python3 -m unittest discover -s tests -v
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/linear/reproduce_all_linear.py
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/crosscheck_linear_vs_semiclassical.py
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/nonlinear/reproduce_all_nonlinear.py
```

