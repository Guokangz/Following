# Nonlinear PRL Fig. 3 validation

Last updated: 2026-08-23

## Status

PRL Fig. 3(a-d) is reproduced from the direct `(n)(m)` hierarchy through `alpha^(2)(1)` and the third-order differential-transmission expression Eq. (8).

Run:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/nonlinear/reproduce_fig3.py
```

Outputs:

- `figures/nonlinear/fig3_reproduction.{png,pdf}`
- `figures/nonlinear/fig3b_delay4_cross_section.{png,pdf}`
- `data/nonlinear/fig3_reproduction.npz`
- `data/nonlinear/fig3_reproduction.json`

## Model and parameters

Equations: PRL (2), (3), (6), (8); SM (S.7)-(S.20).

The calculation uses the paper parameters in `kappa=1` units:

- `G=g sqrt(N)=3`;
- `omega_0=omega_c=omega_p=omega_probe` (zero rotating-frame detuning);
- Gaussian standard-deviation width `tau_w=0.05 kappa^-1`, following the transform-defining SM Eq. (S.20);
- pump arrival `tau_p=2 kappa^-1`;
- `gamma=0`;
- `gamma_phi=0` and `gamma_phi=kappa`;
- delays `tau_Delta kappa` from 0 to 6.

The overall `eta_p^2` factor in Eq. (8) is set to one, so the saved heatmaps are the raw third-order coefficient. No visual rescaling is applied to the arrays.

## Numerical settings and convergence

- time step: `0.005 kappa^-1`;
- window: `0 <= t kappa < 20`;
- the pump is 40 pulse widths after the left boundary;
- even at the largest delay, the probe is followed by 12 cavity-time units;
- positive-sign FFT convention from SM Eq. (S.20);
- zero-padding factor: 4;
- plotted frequency grid: `-5 <= omega/kappa <= 5`.

At representative delays `0,2,4,6`, the maximum absolute differences between `dt=0.005` and `dt=0.0025` are:

- `gamma_phi=0`: `6.59e-7, 3.60e-7, 1.68e-7, 8.40e-8`;
- `gamma_phi=kappa`: `5.03e-8, 6.08e-8, 1.57e-7, 3.77e-7`.

The raw signal ranges are:

- `gamma_phi=0`: `[-0.01797, 0.03377]`;
- `gamma_phi=kappa`: `[-0.01793, 0.02871]`.

These values naturally lie on the approximately `+/-0.03` scale shown in the paper; the match does not come from a forced color limit.

## Physical validation

For `gamma_phi=0`, the heatmap shows alternating-sign upper/lower polariton features near `omega/kappa=+/-3` as the delay changes. The pump-only second-order population and squared first-order polarization oscillate together, as expected for coherent bright polariton dynamics without pure dephasing.

For `gamma_phi=kappa`, the polarization decays rapidly while the molecular population remains. At `t kappa=10`:

- pump population coefficient: `7.48e-3`;
- squared polarization coefficient: `4.28e-6`.

The persistent population produces the long-delay contracted-splitting differential spectrum. The separately saved `tau_Delta=4 kappa^-1` cross-section has positive features inside the bare `+/-3 kappa` polariton guides and negative outer lobes, matching the paper's saturation/contraction interpretation.

The late-delay polariton-window RMS at `tau_Delta kappa=6` is `7.13e-4` without dephasing and `6.05e-3` with strong dephasing; the latter is the population-induced stationary response, not surviving optical coherence.

## Remaining presentation differences

- The paper's exact delay sampling and time-window endpoints are not stated; this reproduction selects them by convergence.
- The pump/probe indicators in panels (c,d) are redrawn from the analytic envelopes rather than copied.
- The color scale is chosen symmetrically from the maximum raw magnitude across both computed heatmaps; raw arrays are preserved.

These are plotting/numerical-window conventions, not unresolved physical discrepancies.

