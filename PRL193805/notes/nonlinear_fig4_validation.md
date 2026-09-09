# Nonlinear Fig. 4 validation

## Model and parameters

PRL Fig. 4(a,b) is computed by
`scripts/nonlinear/reproduce_fig4.py` from the explicit phase hierarchy in
`src/nonlinear/phase_components.py`, not by phase filtering an image or by
fitting a target heatmap.

The parameters are those in the caption:

```text
kappa = 1
sqrt(N) g = 3
gamma_phi = 0.1
gamma = 0
omega_c = omega_0 = omega_p = omega_p' (rotating-frame detunings are zero)
tau_w = 0.1 kappa^-1
```

The Gaussian convention follows SM Eq. (S.20),
`exp[-(t-tau)^2/(2 tau_w^2)]`.  The pump center is placed at
`2 kappa^-1`; only the relative pump-probe delay is physical.

The phase-resolved differential transmission is evaluated from PRL Eq. (10),
using the positive-sign Fourier transform independently validated in
`notes/pump_probe_validation.md`.

## Qualitative reproduction

- `(0,1)` contains a prompt oscillatory bright response and a long-delay
  contribution at the two polariton frequencies.  This is the expected
  bright-plus-dark/population pathway.
- `(2,-1)` has alternating upper/lower-polariton lobes but becomes negligible
  at long delay.  It therefore isolates the finite-lived bright coherence.
- The peaks track the linear-response guides at approximately
  `omega/kappa = +/-3`.
- At delay `8 kappa^-1`, the polariton-window RMS is `1.77e-2` for `(0,1)` but
  only `2.79e-4` for `(2,-1)`.

The raw ranges, before any plotting color limit, are `[-0.0695, 0.0885]` for
`(0,1)` and `[-0.0622, 0.1153]` for `(2,-1)`.  The absolute normalization is
larger than the published colorbar by an order-one pulse-amplitude factor;
the paper does not state the perturbative amplitude `eta_p` entering Eq. (10).
No amplitude was fitted to the paper.

## Bright-state decay time and the factor-of-two convention

The pointwise polariton-window signal undergoes upper/lower-polariton beating,
so a direct log fit has deep periodic minima and is not a valid exponential
estimator.  The script instead fits the local maxima of the RMS envelope from
delay `0.3` to `6.0 kappa^-1`:

```text
tau_fit = 1.8123 kappa^-1
log-space R^2 = 0.9644
```

The PRL describes the scale as
`~(kappa + gamma_phi)^-1 = 0.9091 kappa^-1`.  The explicit SM equations use
cavity- and molecular-amplitude rates `kappa/2` and `gamma_phi/2`.  On
resonance each polariton amplitude therefore decays at
`(kappa+gamma_phi)/4`; the 2LS DQC source contains two residual first-order
amplitudes and its fitted signal envelope decays at

```text
tau_DQC,SM = 2/(kappa + gamma_phi) = 1.8182 kappa^-1.
```

The numerical fit differs from this equation-level prediction by `-0.32%`.
Thus the paper's statement is correctly reproduced as an order-of-magnitude
bright-state scale, while the factor of two follows from its own explicit
amplitude-decay conventions.  This distinction is retained rather than hidden
by redefining `kappa` or `gamma_phi`.

## Numerical convergence

The production grid uses `dt=0.005 kappa^-1`, end time `24 kappa^-1`, and a
fourfold zero-padded FFT.  Repeating selected delays with `dt=0.0025` gives
maximum spectral differences:

| delay | `(0,1)` | `(2,-1)` |
|---:|---:|---:|
| 0 | `4.18e-7` | `2.09e-7` |
| 1 | `3.98e-7` | `1.93e-7` |
| 2 | `4.75e-7` | `1.61e-7` |
| 4 | `1.13e-6` | `9.31e-8` |
| 8 | `4.42e-6` | `2.33e-8` |

The larger relative late-delay `(0,1)` interpolation error is still small on
its `~1e-2` signal scale.

Outputs:

- `figures/nonlinear/fig4_phase_components_reproduction.{png,pdf}`
- `figures/nonlinear/fig4_dqc_decay_fit.{png,pdf}`
- `data/nonlinear/fig4_phase_components_reproduction.{npz,json}`

