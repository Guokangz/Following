# Nonlinear conceptual figures

The script `scripts/nonlinear/reproduce_conceptual_figures.py` reconstructs
the logical content of PRL Figs. 1(b,c), 2, and 4(c,d) as Matplotlib vector
graphics.  No article image, pixel data, or traced geometry is used.

## Fig. 1(b,c)

The first panel contrasts the perturbative equations:

- free space couples `rho^(n)` only to the preceding matter order through the
  external laser field;
- the cavity hierarchy contains the convolution over all lower light and
  matter orders because molecular polarization feeds back into the field.

The pulse-storage panel does not use an illustrative damped sinusoid.  Its
red and blue intracavity traces are the real parts of the first-order pump and
probe fields computed from the resonant 2LS cavity equations with
`G=3 kappa`, `gamma_phi=0.1 kappa`.  They explicitly demonstrate that a short
external pulse produces a long, oscillatory intracavity tail.

## Fig. 2

The pathway tree is generated from the partitions in PRL Eq. (5b).  The
orange branch is the conventional chain of first-order laser interactions.
The black branches show the additional lower-order partitions introduced by
cavity feedback.  The stated `2^(n-1)-1` extra paths follows the paper's
counting.

## Fig. 4(c,d)

The double-sided diagrams label each intermediate 2LS density-matrix state
and each pulse phase.  They show:

- `(0,1) = +/-Phi_p -/+Phi_p + Phi_p'`: ground- and excited-population
  pathways, hence bright and dark contributions;
- `(2,-1) = Phi_p + Phi_p - Phi_p'`: the stored probe acts on pump-created
  coherence, hence a bright-only contribution;
- the reversed probe-first ordering is crossed out because it is absent from
  the differential probe transmission.

Outputs:

- `figures/nonlinear/fig1_perturbative_logic.{png,pdf}`
- `figures/nonlinear/fig2_pathway_tree.{png,pdf}`
- `figures/nonlinear/fig4_feynman_pathways.{png,pdf}`
- `data/nonlinear/fig1_storage_fields.npz`
- `data/nonlinear/conceptual_figures_metadata.json`

