# Linear-response validation

Last updated: 2026-08-23

## Scope and status

The physics-first linear baseline is complete for the main numerical figures of Yuen-Zhou and Koner, J. Chem. Phys. 160, 154107 (2024).

| Figure | Status | Main script | Main outputs |
|---|---|---|---|
| Fig. 3(a,b) | REPRODUCED | `scripts/linear/reproduce_fig3.py` | `figures/linear/fig3_reproduction.{png,pdf}`, separate panel files, thermal-splitting sanity plot |
| Fig. 4(a-c) | REPRODUCED | `scripts/linear/reproduce_fig4.py` | `figures/linear/fig4_reproduction.{png,pdf}`, separate panel files |
| Fig. 5 | REPRODUCED | `scripts/linear/reproduce_fig5.py` | `figures/linear/fig5_reproduction.{png,pdf}`, FC-weight sanity plot |
| Fig. 6(a-c) | REPRODUCED with documented inferred scale/dipoles | `scripts/linear/reproduce_fig6.py` | `figures/linear/fig6_reproduction.{png,pdf}`, susceptibility sanity plot |

All four workflows are run by:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=tmp/matplotlib python3 scripts/linear/reproduce_all_linear.py
```

The raw arrays and JSON metadata are under `data/linear/`.

## Core implementation

The common interface is the molecular susceptibility `chi(omega)`. The retarded photon Green function is

`D^R(omega) = 1 / [omega - omega_ph + i kappa/2 + chi(omega)]`.

The code evaluates all observables independently through Eq. (27):

- `T = kappa_L kappa_R |D^R|^2`
- `R = 1 + 2 kappa_L Im D^R + kappa_L^2 |D^R|^2`
- `A = -kappa_L [kappa |D^R|^2 + 2 Im D^R]`

This makes `T+R+A=1` an independent numerical check. The implementation is in:

- `src/linear/susceptibility.py`
- `src/linear/input_output.py`

## Automated checks

`python3 -m unittest discover -s tests -v` currently runs 10 tests covering:

- balanced empty-cavity Lorentzian and unit resonant transmission;
- no-molecule and no-molecular-damping limits;
- Eq. (34) versus the explicit two-level Eq. (38);
- independent flux conservation and passivity;
- thermal population factor `tanh[DeltaE/(2kBT)]`;
- Gaussian and Lorentzian probability normalization;
- deterministic Gaussian integration versus the analytic Voigt expression;
- Lorentzian disorder linewidth addition;
- Franck-Condon normalization and mean;
- fully saturated 3LS susceptibility `chi=0`.

## Figure-by-figure quantitative validation

### Fig. 3: identical two-level ensemble and thermal saturation

Equations: (34), (37), (38).

Parameters: `omega_ph=omega_exc`, `kappa=0.1`, `gamma=0.3`, `G=sqrt(N)g=2`, with `kappa_L=kappa_R=kappa/2` inferred from the balanced-cavity figure convention.

Results:

- absorption peaks: `-1.997`, `+1.997`;
- absorption peak splitting: `3.994`;
- complex poles: `-1.9993749-0.1i`, `+1.9993749-0.1i`;
- real-part pole splitting: `3.9987498`;
- left/right absorption peak symmetry error: `0` on the chosen grid;
- maximum `|T+R+A-1|`: `2.22e-16`.

Panel (b) uses the actual legend ratios in the paper, `DeltaE/(kBT)=infinity, 1.1, 0.1`. The code scales `G^2` by `tanh[DeltaE/(2kBT)]`. The auxiliary pole-splitting plot verifies the contraction mechanism continuously rather than only at three temperatures.

Agreement: peak locations, symmetry, linewidth scale, peak heights, and high-temperature contraction agree with the paper figure within plotting resolution.

### Fig. 4: disorder and linewidth narrowing

Equations: (34), (41)-(43).

Parameters: `sigma=1`, zero cavity-ensemble detuning, `kappa=gamma=0.1`, `G=1.5` for panels (a,b), and `G=1.3,3.0` for panel (c).

Integration checks:

- Gaussian integral over `[-8 sigma,8 sigma]`: `0.9999999999999987`;
- normalized Lorentzian integral over the real line: `0.9999999999999978`;
- Gaussian deterministic-grid maximum error versus analytic Voigt at `2001,4001,8001,16001` points: `5.59e-9, 4.48e-15, 3.83e-15, 4.83e-15`;
- Lorentzian infinite-range quadrature versus analytic convolution maximum error: `1.78e-15`;
- maximum flux error in both disorder models: `2.22e-16`.

Line narrowing:

| `G` | polariton peak positions | individual absolute half-maximum FWHM |
|---:|---|---|
| 1.3 | `+/-1.4855` | `1.2183` |
| 3.0 | `+/-3.2150` | `0.1557` |

The bare Voigt-broadened Gaussian molecular band has FWHM `2.4086`. Thus the strong-coupling narrowing is quantitative, not inferred from appearance.

Important convention: accepted-manuscript Eq. (43) has `(sigma/2)^2` in its numerator and is not normalized. The reproduction uses the normalized Cauchy numerator `sigma/2`; the test suite explicitly confirms that the printed form integrates to `sigma/2`.

### Fig. 5: vibronic coupling

Equations: (34), (44).

Parameters: `omega_ph=omega_exc`, `kappa=gamma=0.1`, `G=1`, `omega_v=0.3`, `S=3`.

The truncation is chosen algorithmically from an omitted Franck-Condon weight tolerance of `1e-12`:

- selected `m_max=22`;
- retained FC sum: `0.9999999999997933`;
- omitted weight: `2.07e-13`;
- weighted mean vibronic quantum number: `2.999999999995`, as required for a Poisson distribution with mean `S=3`;
- maximum flux error: `2.22e-16`.

Maximum absorption-spectrum errors relative to `m_max=22` are:

| `m_max` | max error |
|---:|---:|
| 6 | `2.64e-1` |
| 9 | `8.85e-4` |
| 12 | `8.61e-6` |
| 15 | `5.47e-8` |
| 18 | `2.10e-10` |
| 22 | `0` |

The reproduced main lower-polariton feature and vibronic comb follow directly from the Eq. (44) transition frequencies and Poisson weights.

### Fig. 6: stationary three-level populations

Equations: (34), (45), (46), evaluated for the three positive-frequency RWA transitions `1-2`, `1-3`, and `2-3`.

Parameters: `omega_ph=omega_12`, `omega_23=2omega_12`, `kappa=0.1`, `gamma=0.3`, `G=1`. Because the paper only states ratios, `omega_12=1` is used as the arbitrary frequency unit. Equal transition dipoles are inferred from the single quoted coupling.

Results:

- case (a), `p=(0.7,0.2,0.1)`: four resolved absorption peaks;
- case (b), `p=(0.48,0.48,0.04)`: three resolved peaks and the `1-2` susceptibility term is exactly zero;
- case (c), equal populations: `max|chi|=0` and the spectrum agrees with the independently evaluated empty cavity with maximum error `0`;
- maximum flux errors across the three cases: `2.22e-16`, `2.22e-16`, `4.44e-16`.

The auxiliary `Im chi` figure directly shows three active population-difference terms, then two, then none.

## Remaining differences and classification

| Difference | Classification | Consequence |
|---|---|---|
| Linear figures do not state `kappa_L/kappa_R` | parameter ambiguity | Balanced mirrors are inferred; changing asymmetry alters absolute T/R but not molecular poles |
| Eq. (43) is non-normalized as printed | unresolved publication typo, resolved by probability normalization | Normalized Cauchy form is used and explicitly documented |
| Fig. 6 omits absolute `omega_12` and individual dipoles | parameter ambiguity | `omega_12=1` and equal relative dipoles are inferred; peak count and saturation physics are robust |
| Exact fonts/colors/panel aspect differ slightly | plotting convention | No physical effect; style was optimized only after validation |
| Bare Fig. 4(c) curve uses the same homogeneous `gamma=0.1`, producing a Voigt rather than a mathematically pure Gaussian | plotting/definition ambiguity | Difference is small; raw Gaussian density is available from the same module if a pure-distribution overlay is desired |

No unresolved numerical approximation remains at a level visible in the reported spectra. The inferred physical conventions above must be carried into the linear-versus-semiclassical cross-check rather than silently retuned.

