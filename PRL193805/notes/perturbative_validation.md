# Perturbative hierarchy validation

Last updated: 2026-08-23

## Implementation

`src/nonlinear/perturbative.py` directly propagates the coefficients in

`alpha(t) = sum_n eta^n alpha^(n)(t)`,

`rho_eg(t) = sum_n eta^n rho_eg^(n)(t)`,

`p_e(t) = sum_n eta^n p_e^(n)(t)`

from order `n=0` through `n=3`. The coefficients are obtained by inserting the expansion into the nonlinear mean-field ODE and collecting powers of `eta`. They are not obtained from polynomial fits to full nonlinear simulations.

For the reduced 2LS variables, the hierarchy is

`dot alpha_n = -(kappa/2+i delta_c) alpha_n - iG s_n - delta_n1 f(t)`,

`dot s_n = -(Gamma_2+i delta_m) s_n - iG alpha_n + 2iG sum_(j=1)^(n-1) alpha_(n-j) p_j`,

`dot p_n = -gamma p_n - 2G Im sum_(j=1)^(n-1) alpha_j^* s_(n-j)`.

The ground-population correction at every positive order is `-p_n`, so the trace of `rho^(0)` is one while the trace of every positive perturbative order is zero. No positivity constraint is imposed on an individual perturbative order.

## Validation workflow

Run:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/nonlinear/validate_perturbative_hierarchy.py
```

Parameters are dimensionless with `kappa=1`, `G=3`, `gamma=0.2`, `gamma_phi=0.1`, resonant drive, and a Gaussian standard-deviation pulse width `0.2 kappa^-1`. The purpose is hierarchy validation, not reproduction of a paper panel.

The full nonlinear mean-field equations are propagated independently for `eta = 0.04, 0.06, 0.09, 0.135, 0.20`, and compared with the reconstruction through third order.

## Results

| `eta` | cavity max error | coherence max error | population max error | total max error |
|---:|---:|---:|---:|---:|
| 0.040 | `5.06e-10` | `4.52e-10` | `3.06e-8` | `3.06e-8` |
| 0.060 | `3.85e-9` | `3.44e-9` | `1.55e-7` | `1.55e-7` |
| 0.090 | `2.92e-8` | `2.61e-8` | `7.84e-7` | `7.84e-7` |
| 0.135 | `2.22e-7` | `1.98e-7` | `3.97e-6` | `3.97e-6` |
| 0.200 | `1.59e-6` | `1.42e-6` | `1.92e-5` | `1.92e-5` |

Log-log fitted error powers:

- cavity field: `eta^5.002`;
- molecular coherence: also fifth order to numerical precision;
- population: `eta^4.001`;
- total state error: `eta^4.001`.

These powers are the correct 2LS selection-rule result. The field/coherence contain odd orders, so after retaining orders 1 and 3 the leading omitted field term is fifth order. The population contains even orders, so after retaining order 2 the leading omitted population term is fourth order.

Directly observed forbidden-order maxima are exactly zero in the implemented hierarchy:

- `max |alpha^(2)| = 0`;
- `max |p_e^(1)| = 0`;
- `max |p_e^(3)| = 0`.

Outputs:

- `figures/nonlinear/perturbative_hierarchy_validation.{png,pdf}`
- `data/nonlinear/perturbative_hierarchy_validation.npz`
- `data/nonlinear/perturbative_hierarchy_validation.json`

This validates N1 and provides the coefficient engine to generalize next to the two-pulse `(n)(m)` hierarchy.

