# Phase-resolved perturbation validation

## Implemented hierarchy

`src/nonlinear/phase_components.py` specializes Supplemental Material Eqs.
(S.39)-(S.43) to a two-level molecule (`g23=0`).  It propagates the first-order
pump and probe fields, the second-order `(0,0)` and `(1,-1)` populations, and
the third-order cavity fields with phase vectors `(0,1)` and `(2,-1)`.

The convention used throughout is

```text
input(Phi_p, Phi_p') = input * exp[-i(Phi_p or Phi_p')]
component(v) = average_Phi alpha(Phi) * exp[+i v.Phi].
```

This is the sign convention in PRL Eq. (9).  Complex-conjugate fields carry
the opposite phase index automatically.

## Independent checks

The validation script is `scripts/nonlinear/validate_phase_components.py`.
It performs two checks which do not share the explicit phase-index
bookkeeping:

1. At zero pulse phases, the sum of the explicit `(0,1)` and `(2,-1)` fields
   must equal the unresolved pump-probe coefficient `alpha^(2)(1)`.
2. Full nonlinear mean-field trajectories are evaluated on a `6 x 6` grid of
   pump/probe phases.  A discrete two-dimensional Fourier transform extracts
   the same two third-order phase components.

The numerical results are:

| check | maximum complex-field error |
|---|---:|
| explicit component sum vs unresolved hierarchy | `1.39e-11` |
| explicit `(0,1)` vs numerical phase cycling | `9.59e-11` |
| explicit `(2,-1)` vs numerical phase cycling | `6.45e-11` |

These errors are at the ODE-tolerance/numerical extraction level.  Thus the
phase signs and the conjugation of the mixed second-order population are not
being inferred from the appearance of Fig. 4.

Raw arrays and metadata are in
`data/nonlinear/phase_component_validation.{npz,json}`.

