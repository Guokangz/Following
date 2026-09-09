# Paper map and reproduction audit

Last updated: 2026-08-23

## 1. Source audit

Local sources inspected in full:

- `ref/Yuen-Zhou和Koner - 2024 - Linear response of molecular polaritons.pdf` (17 pages; accepted manuscript of J. Chem. Phys. 160, 154107, DOI 10.1063/5.0183683).
- `ref/Reitz 等 - 2025 - Nonlinear Semiclassical Spectroscopy of Ultrafast Molecular Polariton Dynamics.pdf` (7 pages; Phys. Rev. Lett. 134, 193803, DOI 10.1103/PhysRevLett.134.193803).
- `ref/supplement_nonlinear_NonMarked.pdf` (11 pages; dated 2025-04-17; sections S1-S12 and Figs. S1-S4).

The local supplemental material contains the Liouville-space mapping, input-output normalization, pump-probe hierarchy, differential-transmission derivation, analytic 2LS calculation, formal integral/pathway expansion, disorder extension, bright/dark population definitions, phase-resolved hierarchy, explicit 3LS equations, and multimode extension. It therefore supplies the details that the PRL main text intentionally omits.

## 2. Logical connection between the papers

### 2.1 Linear-response paper

The linear paper follows this chain:

1. Start from a photon mode coupled to an ensemble of generally anharmonic molecules, Eqs. (1)-(5).
2. Treat the photon as a quantum impurity and the large molecular ensemble as a bath. In the thermodynamic collective limit, higher bath cumulants vanish and an effective harmonic bath with the same two-point correlation function suffices, Eqs. (7)-(20).
3. Identify the effective bath spectral density with the absorptive molecular susceptibility, `J_eff(omega) = hbar Im chi(omega)` for positive frequency, Eqs. (20)-(22).
4. Express transmission, reflection, and absorption through the retarded photon Green function, Eqs. (25)-(27).
5. Sum the molecular self-energy in the effective harmonic problem: `Sigma_M = -chi`, giving

   `D^R(omega) = 1 / [omega - omega_ph + i kappa/2 + chi(omega)]`, Eqs. (32)-(33).

6. Insert this into input-output theory to obtain Eq. (34):

   `T = kappa_L kappa_R / |denominator|^2`,

   `A = 2 kappa_L Im chi / |denominator|^2`,

   `R = 1 - T - A`.

7. Reuse the same susceptibility interface for temperature/saturation (Eqs. 37-38), energetic disorder (Eqs. 41-43), vibronic structure (Eq. 44), and arbitrary stationary three-level populations (Eq. 46).

The scope is broader than thermal equilibrium: the initial molecular state may be any stationary state decoupled from the photon. Population inversion is excluded by the positive-temperature harmonic-bath mapping noted in footnote [40], although the algebraic susceptibility itself can still be evaluated outside that domain with care.

### 2.2 Nonlinear semiclassical paper

The nonlinear paper follows this chain:

1. Start from the Dicke-type Hamiltonian, PRL Eq. (1), and take the collective `N -> infinity` mean-field limit.
2. Propagate a classical intracavity amplitude `alpha(t)` and a quantum single-molecule density matrix `rho(t)` self-consistently, Eqs. (2)-(3):

   `H_MF = H0 + E0 [alpha + alpha*] mu`,

   `d rho/dt = -i[H_MF,rho] + D[rho]`,

   `d alpha/dt = -(kappa/2 + i omega_c) alpha - i N E0 P - input`,

   with `P = Tr(mu rho)`.

3. Expand both light and matter in powers of the incident field, PRL Eq. (4), producing the coupled order-by-order hierarchy Eq. (5). Molecular feedback makes a given density-matrix order depend on all lower cavity-field orders, unlike ordinary free-space spectroscopy.
4. Generalize to pump and probe orders `(n)(m)`, PRL Eq. (6) and SM Eqs. (S.7)-(S.10).
5. Use balanced-cavity input-output relations, SM Eqs. (S.3)-(S.6), and isolate the lowest third-order correction to differential transmission, PRL Eq. (8) / SM Eq. (S.13), through interference of `alpha^(0)(1)` with `alpha^(2)(1)`.
6. Expand in pump/probe phases, PRL Eq. (9), SM Eqs. (S.29)-(S.35), to separate the `(0,1)` pump-probe and `(2,-1)` double-quantum-coherence contributions.
7. Interpret the components as Liouville-space pathways. The `(0,1)` component can retain dark-state population, while `(2,-1)` for a 2LS requires cavity pulse storage and decays with the bright coherence on a scale of order `(kappa + gamma_phi)^-1`.

### 2.3 The required linear/nonlinear interface

The PRL explicitly states after Eq. (5) that first order recovers the standard linear polariton response and cites the 2024 linear-response paper as Ref. [61]. For a resonant 2LS initialized in its ground state, the positive-frequency/RWA linearized mean-field equations are

`d alpha/dt = -(kappa/2 + i omega_c) alpha - i N g rho_eg - input`,

`d rho_eg/dt = -[(gamma + gamma_phi)/2 + i omega_0] rho_eg - i g alpha`.

Fourier transformation gives the same denominator as linear Eqs. (34)/(38), provided the molecular coherence linewidth is identified as

`gamma_linear/2 = (gamma + gamma_phi)/2`.

Thus `gamma_linear = gamma + gamma_phi` under the SM convention. This is a convention identity, not an adjustable fit. The explicit cross-check must also use the same Fourier sign, balanced mirror losses, and collective coupling `G = g sqrt(N)`.

## 3. Figure reproduction audit

### 3.1 Linear-response paper

| Figure | Panel | Observable / purpose | Governing equations | Parameters | Numerical method | Expected feature | Difficulty | SM needed? |
|---|---:|---|---|---|---|---|---|---|
| 1 | - | Conceptual impurity mapping: anharmonic molecular ensemble to an effective harmonic bath | (1)-(20) | `N >> 1`, individual coupling scales as `N^-1/2` | Conceptual redraw only | Same photon reduced dynamics when two-point bath correlations match | Low, after numerics | No |
| 2 | - | Input/output channels and definitions of T/R/A | (25)-(27) | `kappa=kappa_L+kappa_R` | Conceptual redraw plus Eq. (27) checks | Flux partition into transmission, reflection, molecular absorption | Low | No |
| 3 | a | T, R, A for identical 2LS ensemble | (34), (37), (38) | `omega_ph=omega_exc`, `kappa=0.1`, `gamma=0.3`, `G=sqrt(N)g=2`; balanced cavity inferred | Direct complex susceptibility on dense frequency grid | Symmetric UP/LP features near detunings `+/-2`; `T+R+A=1` | Low | No |
| 3 | b | Absorption versus temperature | (37)-(38) | Same as 3a; legend gives `DeltaE/(kBT)=infinity, 1.1, 0.1` | Replace `G^2` by `G_0^2 tanh[DeltaE/(2kBT)]` | Rabi splitting contracts continuously and tends to empty-cavity response at high T | Low | No |
| 4 | a | T, R, A with Gaussian disorder | (34), (41), (42) | Gaussian standard deviation `sigma=1`, zero detuning, `kappa=gamma=0.1`, `G=1.5` | Deterministic quadrature/grid with convergence and distribution-normalization checks | Broad molecular band with two polaritonic features | Medium | No |
| 4 | b | T, R, A with Lorentzian disorder | (34), (41), (43) | Lorentzian width parameter `sigma=1`; otherwise as 4a | Prefer analytic Cauchy convolution, benchmarked by quadrature | Lorentzian tails and broader central response than Gaussian | Medium due Eq. (43) issue | No |
| 4 | c | Normalized cavity absorption and bare Gaussian absorption | (34), (41), (42) | `sigma=1`, `G=1.3` and `3.0`, `kappa=gamma=0.1` | Same Gaussian integral; normalize each displayed curve as caption indicates; compute FWHM separately from raw curves | Strong-coupling polariton linewidth narrowing for `G >> sigma` | Medium | No |
| 5 | - | T, R, A with vibronic progression | (34), (44) | `omega_ph=omega_exc`, `kappa=gamma=0.1`, `G=1`, `omega_v=0.3`, `S=3` | Sum Poisson Franck-Condon weights to converged `m_max` | Main lower polariton plus a vibronic comb; FC sum tends to one | Medium | No |
| 6 | a | T, R, A for populated 3LS | (34), (45), (46) | `p=(0.7,0.2,0.1)`, `omega_ph=omega_12`, `omega_23=2 omega_12`, `kappa=0.1`, `gamma=0.3`, `G=1` | Sum all uphill pair susceptibilities; equal transition dipoles inferred | Three active transitions hybridize with one cavity into four polariton branches | Medium | No |
| 6 | b | One saturated 3LS transition | (34), (46) | `p=(0.48,0.48,0.04)`; other values as 6a | Same | `1<->2` contribution vanishes; three-branch structure remains | Medium | No |
| 6 | c | Fully saturated 3LS | (34), (46) | `p1=p2=p3=1/3` | Same | `chi=0`; exact empty balanced-cavity spectrum | Low | No |

### 3.2 Nonlinear PRL

| Figure | Panel | Observable / purpose | Governing equations | Parameters | Numerical method | Expected feature | Difficulty | SM needed? |
|---|---:|---|---|---|---|---|---|---|
| 1 | a | Cavity/pump/probe schematic | PRL (1)-(3) | Conceptual | SVG/matplotlib redraw | Incident pulses drive cavity, not molecules directly | Low | No |
| 1 | b | Free-space versus cavity perturbative coupling | PRL (4)-(5) | Conceptual | Equation/graph redraw | Cavity feedback couples `rho^(n)` to multiple lower `alpha` orders | Low | Helpful |
| 1 | c | Pulse storage and order reversal | PRL (3), (6) | Storage time `kappa^-1` | Draw from computed first-order fields | Intracavity probe tail can precede/overlap pump interactions | Medium | Yes, S3 |
| 2 | - | Density-matrix pathway tree | PRL (5), formal integrals SM (S.21)-(S.22) | General order `n` | Programmatic graph | `2^(n-1)-1` extra cavity-feedback paths in addition to free-space path | Medium | Yes, S6 |
| 3 | a | Total third-order DT heatmap, no pure dephasing | PRL (2)-(8), SM (S.7)-(S.20) | `G=3 kappa`, `gamma_phi=0`, `tau_w=0.05 kappa^-1`, resonance, `tau_p=2 kappa^-1`; `gamma=0` supported by Fig. S1 | Propagate `(n)(m)` hierarchy through `(2)(1)`, FFT, store raw arrays | Long-lived oscillatory UP/LP coherence; features at rotating-frame `omega/kappa approximately +/-3` | High | Essential |
| 3 | b | Total third-order DT heatmap, strong dephasing | Same | `gamma_phi=kappa`; otherwise 3a | Same | Coherence decays quickly; persistent molecular population contracts splitting; cross-section at `tau_Delta=4 kappa^-1` | High | Essential |
| 3 | c | Pump-only population and `|P(t)|^2`, no dephasing | PRL (2)-(5) | As 3a | Pump-only perturbative/full mean-field propagation | Rabi oscillations decay with cavity leakage | Medium | Helpful |
| 3 | d | Pump-only population and `|P(t)|^2`, dephasing | PRL (2)-(5) | As 3b | Same | Polarization decays, molecular population remains | Medium | Helpful |
| 4 | a | Phase component `(0,1)` of DT | PRL (9)-(10), SM (S.29)-(S.35) | `G=3 kappa`, `gamma_phi=0.1 kappa`, resonance, `tau_w=0.1 kappa^-1` | Explicit phase-index hierarchy, benchmark with phase cycling | Bright early-time and dark long-time contributions | Very high | Essential |
| 4 | b | Phase component `(2,-1)` of DT | Same | Same | Same plus delay-envelope fit | Bright-only signal decays on order `(kappa+gamma_phi)^-1` | Very high | Essential |
| 4 | c | Double-sided diagrams for `(0,1)` | SM phase rules | Conceptual | Programmatic diagram | Pump phases cancel; population/zero-quantum paths | Medium | Essential |
| 4 | d | Double-sided diagrams for `(2,-1)` | SM phase rules | Conceptual | Programmatic diagram | Cavity storage allows probe to act on pump coherence; reversed timing is excluded from probe DT | Medium | Essential |

### 3.3 Supplemental figures

| Figure | Observable | Equations / parameters | Role in this project |
|---|---|---|---|
| S1 | Analytic Fourier-domain versus time-domain DT for 2LS | (S.14)-(S.20), `G=3kappa`, `tau_w=0.05kappa^-1`, resonance, `gamma_phi=gamma=0` | Mandatory benchmark for FFT and hierarchy implementation |
| S2 | Narrowband lower-polariton pump-probe | S7 discussion, `G=3kappa`, `tau_w=kappa^-1` | Optional extension after main Fig. 3; verifies disappearance of UP/LP beating |
| S3 | Disorder dependence of phase components and bright lifetime | two equal subensembles, `delta=0,0.3,0.6 kappa`, `G=3kappa`, `gamma_phi=0.2kappa` | Optional validation of generalized ensemble code |
| S4 | 3LS total and phase-resolved DT | (S.36)-(S.43), `G_12=3kappa`, `g_23=2g_12`, `omega_3-omega_2=1.8omega_2`, `gamma_phi=0.1kappa` | Later extension; demonstrates that 3LS DQC can persist beyond the 2LS bright lifetime |

## 4. Unified parameter and convention table

| Quantity | Linear paper | Nonlinear paper / SM | Status and implementation note |
|---|---|---|---|
| `omega_ph`, `omega_c` | Cavity frequency | Cavity frequency | Same physical quantity; spectra may use lab frame or detuning/rotating frame |
| `omega_exc`, `omega_0` | 2LS transition frequency | 2LS transition frequency | Same physical quantity |
| `kappa` | Total photon escape `kappa_L+kappa_R`; denominator `+i kappa/2` | Field equation decays as `kappa/2`; lifetime discussion uses `kappa^-1` | Compatible |
| `kappa_L`, `kappa_R` | General in Eqs. (27),(34) | SM assumes balanced cavity for displayed T: `kappa_L=kappa_R=kappa/2` | Linear figures also appear balanced; mark as inferred there |
| `gamma` | Molecular full coherence linewidth; denominator `+i gamma/2` | Spontaneous emission parameter; coherence decay `(gamma+gamma_phi)/2` in SM (S.14) | Cross-check uses `gamma_linear=gamma+gamma_phi` |
| `gamma_phi` | Not separated | Pure dephasing parameter | Fig. 3: 0 and `kappa`; Fig. 4: `0.1kappa`; spontaneous `gamma=0` unless stated |
| `g sqrt(N)` | Collective coupling `G`; Fig. 3 `2`, Fig. 4 `1.5` / `1.3` / `3`, Fig. 5/6 `1` | `3kappa` in main nonlinear figures | Polariton splitting is approximately `2G` at resonance for weak damping |
| `tau_w` | N/A | `0.05kappa^-1` Fig. 3; `0.1kappa^-1` Fig. 4; `kappa^-1` Fig. S2 | The SM has a notational inconsistency: prose below (S.7) prints `exp[-(t-tau)^2/tau_w^2]`, while Fourier Eq. (S.20) uses `exp[-(t-tau)^2/(2tau_w^2)]`. Use Eq. (S.20) as the transform-defining convention and document it in nonlinear validation. |
| Pump/probe carrier | N/A | `omega_p=omega_p'=omega_0=omega_c` in main figures | Results plotted in frame rotating at carrier |
| Pulse timing | N/A | Fig. 3 pump `tau_p=2kappa^-1`, probe delay `tau_Delta=tau_p'-tau_p>=0` | Absolute start/end windows missing; infer only after convergence studies |
| Disorder `sigma` | Fig. 4 `sigma=1` | SM S3 uses a two-subensemble shift `delta` | Gaussian `sigma` is standard deviation |
| Huang-Rhys `S` | Fig. 5 `S=3` | N/A | FC weights are Poisson with mean 3 |
| Vibrational frequency | Fig. 5 `omega_v=0.3` | N/A | Vibronic transition `omega_exc-S omega_v+m omega_v` |
| 3LS populations | Fig. 6: `(0.7,0.2,0.1)`, `(0.48,0.48,0.04)`, `(1/3,1/3,1/3)` | Ground-state 2LS in main nonlinear figures | Stationary linear susceptibility weights are population differences |
| Frequency units | Arbitrary throughout figures | Frequencies scaled by `kappa`; times by `kappa^-1` | Do not assign SI units |
| Fourier convention | Linear Eq. (23): `f(omega)=-i integral exp(i omega t) f(t) dt` for Green functions | SM Eq. (S.20): symmetric `1/sqrt(2pi) integral exp(i omega t) f(t) dt` for pulse transform | The extra Green-function `-i` is not a generic FFT prefactor; cross-check signs at the response level |
| T normalization | Eq. (34), incident from left | SM Eq. (S.6), balanced cavity | Identical after drive-amplitude normalization |

## 5. Explicit ambiguities and suspected typographical issues

1. **Linear Eq. (43) is not normalized as printed.** It shows

   `p(x) = (1/pi) (sigma/2)^2 / [(x-xbar)^2 + (sigma/2)^2]`,

   whose integral is `sigma/2`, not one. Because the surrounding text calls `p` a probability distribution and the numerical quality requirements demand normalization, the implementation should use the normalized Cauchy density with numerator `sigma/2`. Both the printed and normalized forms should be testable so the choice is explicit.

2. **Balanced mirrors in linear figures are not stated in their captions.** Their unit-height empty-cavity transmission behavior and the plotting conventions strongly indicate `kappa_L=kappa_R=kappa/2`. This is an inference and will be recorded in every data metadata file.

3. **Absolute transition scale in linear Fig. 6 is not stated.** Only the ratios `omega_ph=omega_12` and `omega_23=2omega_12` are given. Setting `omega_12=1` is the natural arbitrary-unit choice but must be marked inferred.

4. **Transition dipoles in Fig. 6 are not separately specified.** Use equal `|lambda mu_yz|=g` for all three transitions, consistent with the single collective coupling quoted in the caption; mark inferred.

5. **Nonlinear pulse-envelope convention differs between SM prose and Eq. (S.20).** The latter is analytically self-consistent with its stated transform and is the preferred convention.

6. **Nonlinear absolute pulse amplitude and time window are absent from captions.** Perturbative spectra have an overall amplitude factor; raw arrays should preserve the chosen perturbative normalization. Time window, timestep, padding, and pulse truncation must be selected by convergence, not visual matching.

7. **SM Eq. (S.16) appears to label the Fourier-transformed quantity with a time derivative on the left in extracted text.** The intended algebraic result is clear from the coupled Fourier-domain equations; it will be rederived before implementation rather than copied literally.

## 6. Reproduction dependency roadmap

1. Implement susceptibilities and Eq. (34), with independent Eq. (27) flux checks.
2. Reproduce linear Figs. 3-6 and quantify poles, peaks, FWHM, integral/truncation convergence.
3. Freeze a linear validation baseline and one-command reproduction script.
4. Implement 2LS mean-field dynamics in a rotating frame, then prove its weak-field denominator matches Eq. (34)/(38).
5. Implement genuine order hierarchy through third order and validate perturbative scaling against the full nonlinear solution.
6. Implement pump-probe hierarchy and validate time-domain FFT against SM Fig. S1 / Eqs. (S.14)-(S.20).
7. Reproduce nonlinear Fig. 3, then add explicit phase indices and phase-cycling benchmark for Fig. 4.
8. Draw conceptual/pathway figures only after the numerical pathway structure has been validated.

