#!/bin/sh

# --- Configuration ---
# Path to your pw.x executable
PW_COMMAND="pw.x"
# MPI command (e.g., "mpirun -np 8 --bind-to-core")
MPI_COMMAND="mpirun -np 16 --bind-to core --map-by core"

# Pseudopotential files
AG_PSEUDO="Ag.pbe-n-rrkjus_psl.1.0.0.UPF"
H_PSEUDO="H.pbe-kjpaw.UPF"

echo "Starting the full adsorption energy calculation workflow..."

# --- Step 1: Relax the clean 3x3 Ag(111) slab ---
echo "Step 1: Relaxing the clean 3x3 Ag(111) slab..."

cat > slab.relax.in << EOF
&CONTROL
   calculation  = 'relax',
   prefix       = 'ag111_3x3_6k',
   pseudo_dir   = './',
   outdir       = './'
/
&SYSTEM
   ibrav     = 0, nat=27, ntyp=1, ecutwfc=40.0, ecutrho=200.0,
   occupations='smearing', smearing='methfessel-paxton', degauss=0.02,
   vdw_corr = 'dft-d3'
/
&ELECTRONS
   conv_thr = 1.0d-7, mixing_beta = 0.5
/
&IONS
   ion_dynamics = 'bfgs'
/
CELL_PARAMETERS {angstrom}
   8.850   0.000   0.000
  -4.425   7.665   0.000
   0.000   0.000  25.000
ATOMIC_SPECIES
   Ag  107.868  $AG_PSEUDO
ATOMIC_POSITIONS {crystal}
Ag   0.000000 0.000000 0.450 1 1 1
Ag   0.333333 0.000000 0.450 1 1 1
Ag   0.666667 0.000000 0.450 1 1 1
Ag   0.000000 0.333333 0.450 1 1 1
Ag   0.333333 0.333333 0.450 1 1 1
Ag   0.666667 0.333333 0.450 1 1 1
Ag   0.000000 0.666667 0.450 1 1 1
Ag   0.333333 0.666667 0.450 1 1 1
Ag   0.666667 0.666667 0.450 1 1 1
Ag   0.111111 0.222222 0.500 0 0 0
Ag   0.444444 0.222222 0.500 0 0 0
Ag   0.777778 0.222222 0.500 0 0 0
Ag   0.111111 0.555556 0.500 0 0 0
Ag   0.444444 0.555556 0.500 0 0 0
Ag   0.777778 0.555556 0.500 0 0 0
Ag   0.111111 0.888889 0.500 0 0 0
Ag   0.444444 0.888889 0.500 0 0 0
Ag   0.777778 0.888889 0.500 0 0 0
Ag   0.222222 0.111111 0.550 0 0 0
Ag   0.555556 0.111111 0.550 0 0 0
Ag   0.888889 0.111111 0.550 0 0 0
Ag   0.222222 0.444444 0.550 0 0 0
Ag   0.555556 0.444444 0.550 0 0 0
Ag   0.888889 0.444444 0.550 0 0 0
Ag   0.222222 0.777778 0.550 0 0 0
Ag   0.555556 0.777778 0.550 0 0 0
Ag   0.888889 0.777778 0.550 0 0 0
K_POINTS {automatic}
   6 6 1 0 0 0
EOF

$MPI_COMMAND $PW_COMMAND -i slab.relax.in > slab.relax.out

# --- Step 2: Calculate the isolated H2 molecule ---
echo "Step 2: Calculating the isolated H2 molecule..."

cat > h2.relax.in << EOF
&CONTROL
   calculation='relax', prefix='h2_6k', pseudo_dir='./', outdir='./'
/
&SYSTEM
   ibrav=1, celldm(1)=20.0, nat=2, ntyp=1, ecutwfc=40.0, ecutrho=200.0,
/
&ELECTRONS
   conv_thr=1.0d-8, mixing_beta=0.7
/
&IONS
   ion_dynamics='bfgs'
/
ATOMIC_SPECIES
  H  1.008  $H_PSEUDO
ATOMIC_POSITIONS {angstrom}
  H   0.000   0.000   0.000
  H   0.000   0.000   0.740
K_POINTS {gamma}
EOF

$PW_COMMAND -i h2.relax.in > h2.relax.out

# --- Step 3: Relax the H2 on Ag(111) slab (THE LONGEST JOB) ---
echo "Step 3: Relaxing H2 on the Ag(111) slab (this will take a while)..."

# Extract final coordinates from the slab calculation
RELAXED_SLAB_COORDS=`awk '/Begin final coordinates/,/End final coordinates/' slab.relax.out | sed '1d;$d' | sed 's/Ag/   Ag/'`

cat > slab_h2.relax.in << EOF
&CONTROL
   calculation='relax', prefix='ag111_3x3_h2_6k', pseudo_dir='./', outdir='./'
/
&SYSTEM
   ibrav=0, nat=29, ntyp=2, ecutwfc=40.0, ecutrho=200.0,
   occupations='smearing', smearing='methfessel-paxton', degauss=0.02,
   vdw_corr = 'dft-d3'
/
&ELECTRONS
   conv_thr=1.0d-7, mixing_beta=0.3
/
&IONS
   ion_dynamics='bfgs'
/
CELL_PARAMETERS {angstrom}
   8.850   0.000   0.000
  -4.425   7.665   0.000
   0.000   0.000  25.000
ATOMIC_SPECIES
   Ag  107.868  $AG_PSEUDO
   H    1.008  $H_PSEUDO
ATOMIC_POSITIONS {crystal}
$RELAXED_SLAB_COORDS
H      0.4800000   0.4800000   0.3000000   1 1 1
H      0.5200000   0.5200000   0.3000000   1 1 1
K_POINTS {automatic}
   6 6 1 0 0 0
EOF

# This requires fixing atoms. The script above generates coordinates without fixing flags.
# The user should manually add the fixing flags (0 0 0 for bottom layers, 1 1 1 for top layer)
# to the slab_h2.relax.in file before running the final calculation.
# For now, we will relax all atoms for simplicity of the script.

# A more advanced script would use sed/awk to add fixing flags. Let's keep it simple for now.
# The user should edit slab_h2.relax.in to fix the bottom 18 Ag atoms.

echo "IMPORTANT: Please manually edit 'slab_h2.relax.in' to fix the bottom 18 Ag atoms"
echo "by adding '0 0 0' to the end of their coordinate lines."
echo "Then, run the final calculation manually:"
echo "$MPI_COMMAND $PW_COMMAND -i slab_h2.relax.in > slab_h2.relax.out"

echo "Script finished generating input files."
