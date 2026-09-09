# Pump-probe hierarchy and Fourier validation

Last updated: 2026-08-23

## Direct `(n)(m)` hierarchy

`src/nonlinear/pump_probe.py` implements the two-pulse expansion of PRL Eq. (6) and SM Eqs. (S.7)-(S.10) directly in reduced 2LS Bloch variables. The propagated coefficient arrays include all orders through pump order 2 and probe order 1, in particular:

- the linear probe field `alpha^(0)(1)`;
- pump population `p_e^(2)(0)`;
- mixed population `p_e^(1)(1)`;
- the required nonlinear field `alpha^(2)(1)`.

No multi-amplitude polynomial fit is used to obtain these coefficients.

## Fourier convention

`src/nonlinear/fourier.py` implements

`f(omega) = (1/sqrt(2pi)) integral dt exp(+i omega t) f(t)`,

which is the convention in SM Eq. (S.20). For

`f(t)=exp[-(t-tau)^2/(2 tau_w^2)]`,

the analytic transform is

`f(omega)=tau_w exp[-omega^2 tau_w^2/2] exp[i omega tau]`.

The zero-padded FFT agrees with this analytic result to maximum absolute error `4.94e-16` for `|omega|<15` in the validation workflow. A separate unit test compares selected FFT bins with direct positive-sign quadrature, locking both the sign and frequency direction.

This also resolves the SM prose ambiguity in favor of Eq. (S.20): `tau_w` is used as the Gaussian standard deviation, i.e. the time-domain exponent contains `1/(2 tau_w^2)`.

## Independent full-solution benchmark

Run:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/nonlinear/validate_pump_probe_hierarchy.py
```

For `kappa=1`, `G=3`, `gamma=0.2`, `gamma_phi=0.1`, pump/probe centers 2 and 3, and standard-deviation width `0.2`, the hierarchy gives:

- `max |alpha^(0)(1)| = 0.30149`;
- `max |alpha^(2)(1)| = 0.06169`;
- forbidden even-total-order field `max |alpha^(1)(1)| = 0`.

As an independent check, full nonlinear dynamics with pump amplitude `0.05` and probe amplitude `0.001` is compared after subtracting the pump-only full trajectory. The directly reconstructed probe-linear hierarchy

`eta_probe [alpha^(0)(1) + eta_pump^2 alpha^(2)(1)]`

agrees with the full probe-induced cavity field to maximum absolute error `3.88e-9`.

Outputs:

- `figures/nonlinear/pump_probe_hierarchy_validation.{png,pdf}`
- `data/nonlinear/pump_probe_hierarchy_validation.npz`
- `data/nonlinear/pump_probe_hierarchy_validation.json`

This validates N2's coefficient bookkeeping and FFT convention before constructing the nonlinear Fig. 3 delay-frequency heatmap.

