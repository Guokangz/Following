(* Core analytic checks for the molecular-polariton group-talk manuscript. *)

ClearAll["Global`*"];
$Assumptions = Element[{omega0, omegac, G, gamma, kappa, field,
    populationDifference, input}, Reals] && G > 0 && gamma >= 0 &&
    kappa >= 0;

Print["CHECK 1: coupled-mode eigenvalues"];
hamiltonian = {{omegac, G}, {G, omega0}};
eigenvalues = FullSimplify[Eigenvalues[hamiltonian]];
resonantEigenvalues = FullSimplify[
  eigenvalues /. omegac -> omega0,
  Assumptions -> G > 0
];
Print["general = ", eigenvalues];
Print["resonance = ", resonantEigenvalues];
Print["resonance check = ",
  FullSimplify[Sort[resonantEigenvalues] == Sort[{omega0 - G, omega0 + G}],
    Assumptions -> G > 0]
];

Print["CHECK 2: weak-drive 2LS optical coherence"];
coherenceEquation =
  -I omega coherence == -(gamma/2 + I omega0) coherence +
    I dipole field populationDifference;
coherenceSolution = FullSimplify[
  coherence /. First[Solve[coherenceEquation, coherence]]
];
coherenceTarget =
  -dipole field populationDifference/(omega - omega0 + I gamma/2);
Print["coherence = ", coherenceSolution];
Print["denominator/sign check = ",
  FullSimplify[coherenceSolution == coherenceTarget]
];

Print["CHECK 3: insert the 2LS susceptibility into the photon Green function"];
chi2LS = -G^2/(omega - omega0 + I gamma/2);
green = 1/(omega - omegac + I kappa/2 + chi2LS);
greenTarget = (omega - omega0 + I gamma/2)/(
   (omega - omegac + I kappa/2) (omega - omega0 + I gamma/2) - G^2
);
Print["D^R = ", Factor[Together[green]]];
Print["rational-form check = ", FullSimplify[green == greenTarget]];

Print["CHECK 4: polariton poles and the resonant lossless limit"];
polePolynomial =
  (omega - omegac + I kappa/2) (omega - omega0 + I gamma/2) - G^2;
generalPoles = FullSimplify[omega /. Solve[polePolynomial == 0, omega]];
losslessResonantPoles = FullSimplify[
  generalPoles /. {omegac -> omega0, kappa -> 0, gamma -> 0},
  Assumptions -> G > 0
];
Print["general poles = ", generalPoles];
Print["lossless resonance = ", losslessResonantPoles];
Print["pole check = ",
  FullSimplify[
    Sort[losslessResonantPoles] == Sort[{omega0 - G, omega0 + G}],
    Assumptions -> G > 0
  ]
];

Print["CHECK 5: eliminate the weak-field molecular coherence from mean field"];
meanFieldEquations = {
  -I omega alpha == -(kappa/2 + I omegac) alpha - I G sigma - input,
  -I omega sigma == -(gamma/2 + I omega0) sigma - I G alpha
};
alphaSolution = FullSimplify[
  alpha /. First[Solve[meanFieldEquations, {alpha, sigma}]]
];
greenFromMeanField = FullSimplify[I alphaSolution/input];
Print["alpha/input = ", Factor[Together[alphaSolution/input]]];
Print["i alpha/input = ", Factor[Together[greenFromMeanField]]];
Print["linear-paper denominator check = ",
  FullSimplify[greenFromMeanField == green]
];

allChecks = {
  FullSimplify[Sort[resonantEigenvalues] == Sort[{omega0 - G, omega0 + G}],
    Assumptions -> G > 0],
  FullSimplify[coherenceSolution == coherenceTarget],
  FullSimplify[green == greenTarget],
  FullSimplify[
    Sort[losslessResonantPoles] == Sort[{omega0 - G, omega0 + G}],
    Assumptions -> G > 0],
  FullSimplify[greenFromMeanField == green]
};
Print["ALL CHECKS PASSED = ", And @@ allChecks];
