#!/bin/bash

# --- User Configuration ---
# This section contains parameters you should check and modify.

# Path to your pw.x executable.
PW_COMMAND="pw.x"

# MPI command for parallel execution.
MPI_COMMAND="mpirun -np 16 --bind-to core --map-by core"

# Pseudopotential filenames (must be present in this directory).
AG_PSEUDO="Ag.pbe-n-rrkjus_psl.1.0.0.UPF"
H_PSEUDO="H.pbe-kjpaw.UPF"

# --- Scan Configuration ---
# This section defines the scan range.

# Total number of points to calculate on the curve.
N_POINTS=10
# Starting distance in Angstrom.
START_DIST=2.0
# Ending distance in Angstrom.
END_DIST=5.0

# --- Physical Constants ---
# These values are derived from your PREVIOUSLY relaxed clean slab calculation.

# The length of the c-vector of your cell in Angstrom.
CELL_Z=25.0

# IMPORTANT: The average fractional z-coordinate of your TOPMOST Ag layer.
# PLEASE UPDATE THIS VALUE from your 6k relaxed slab .out file!
Z_AG_TOP=0.344028

# --- Main Script Logic (No need to change anything below) ---

echo "Starting Potential Energy Surface (PES) scan..."
echo "Scanning from $START_DIST A to $END_DIST A in $N_POINTS steps."

# Loop from 0 to N_POINTS-1
for i in $(seq 0 $((N_POINTS-1)))
do
    # Calculate the current H2 distance from the surface in Angstrom.
    CURRENT_DIST=$(echo "scale=10; $START_DIST + ($END_DIST - $START_DIST) * $i / ($N_POINTS - 1)" | bc -l)
    
    # Calculate the corresponding fractional Z coordinate for the H2 molecule's center.
    Z_H2_CENTER=$(echo "scale=10; $Z_AG_TOP + $CURRENT_DIST / $CELL_Z" | bc -l)

    # Calculate Z for each H atom assuming a 0.75 Angstrom bond length perpendicular to the surface.
    H1_Z=$(echo "scale=10; $Z_H2_CENTER - (0.75 / 2 / $CELL_Z)" | bc -l)
    H2_Z=$(echo "scale=10; $Z_H2_CENTER + (0.75 / 2 / $CELL_Z)" | bc -l)
    
    # Create a clean filename label from the distance (e.g., 2.00 -> 2p00).
    DIST_LABEL=$(printf "%.2f" $CURRENT_DIST | tr '.' 'p')
    
    PREFIX="ag_h2_${DIST_LABEL}A_6k"
    INPUT_FILE="${PREFIX}.relax.in"
    OUTPUT_FILE="${PREFIX}.relax.out"

    echo "-----------------------------------------------------"
    echo "Step $((i+1))/$N_POINTS: Calculating for distance = $CURRENT_DIST A"
    echo "-----------------------------------------------------"

    # Create the input file for the current distance point using a "here document".
    cat > $INPUT_FILE << EOF
&CONTROL
   calculation  = 'relax',
   prefix       = '$PREFIX',
   pseudo_dir   = './',
   outdir       = './'
/
&SYSTEM
   ibrav=0, nat=29, ntyp=2, ecutwfc=40.0, ecutrho=200.0,
   occupations='smearing', smearing='methfessel-paxton', degauss=0.02,
   vdw_corr = 'dft-d3'
/
&ELECTRONS
   conv_thr    = 1.0d-8,
   mixing_beta = 0.1,
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
   H    1.008  $H_PSEUDO
ATOMIC_POSITIONS {crystal}
Ag   0.2222220000   0.1111110000   0.2000000000   0 0 0
Ag   0.5555560000   0.1111110000   0.2000000000   0 0 0
Ag   0.8888890000   0.1111110000   0.2000000000   0 0 0
Ag   0.2222220000   0.4444440000   0.2000000000   0 0 0
Ag   0.5555560000   0.4444440000   0.2000000000   0 0 0
Ag   0.8888890000   0.4444440000   0.2000000000   0 0 0
Ag   0.2222220000   0.7777780000   0.2000000000   0 0 0
Ag   0.5555560000   0.7777780000   0.2000000000   0 0 0
Ag   0.8888890000   0.7777780000   0.2000000000   0 0 0
Ag   0.1111110000   0.2222220000   0.2500000000   0 0 0
Ag   0.4444440000   0.2222220000   0.2500000000   0 0 0
Ag   0.7777780000   0.2222220000   0.2500000000   0 0 0
Ag   0.1111110000   0.5555560000   0.2500000000   0 0 0
Ag   0.4444440000   0.5555560000   0.2500000000   0 0 0
Ag   0.7777780000   0.5555560000   0.2500000000   0 0 0
Ag   0.1111110000   0.8888890000   0.2500000000   0 0 0
Ag   0.4444440000   0.8888890000   0.2500000000   0 0 0
Ag   0.7777780000   0.8888890000   0.2500000000   0 0 0
Ag   0.0000029568   0.0000059136   0.3440280909   0 0 0
Ag   0.3333370134   0.0000060547   0.3440276273   0 0 0
Ag   0.6666690413   0.0000060547   0.3440276273   0 0 0
Ag   0.0000030712   0.3333395735   0.3440283118   0 0 0
Ag   0.3333365023   0.3333395735   0.3440283118   0 0 0
Ag   0.6666702255   0.3333394510   0.3440279663   0 0 0
Ag   0.0000043613   0.6666747760   0.3440281131   0 0 0
Ag   0.3333371438   0.6666752877   0.3440280446   0 0 0
Ag   0.6666704147   0.6666747760   0.3440281131   0 0 0
H      0.500000   0.500000   $H1_Z   1 1 0
H      0.500000   0.500000   $H2_Z   1 1 0
K_POINTS {automatic}
   2 2 1 0 0 0
EOF

    # Run the calculation for the current point
    $MPI_COMMAND $PW_COMMAND -i $INPUT_FILE > $OUTPUT_FILE

done

echo "-----------------------------------------------------"
echo "All calculations are finished."
echo "You can now extract the total energy from each .out file."
echo "Example: grep '!' *.out"
echo "-----------------------------------------------------"