#!/bin/bash

# --- User Configuration ---
# Please double-check these paths and filenames before running.

PW_COMMAND="pw.x"
MPI_COMMAND="mpirun -np 16 --bind-to core --map-by core"
AG_PSEUDO="Ag.pbe-n-rrkjus_psl.1.0.0.UPF"
H_PSEUDO="H.pbe-kjpaw.UPF"

# --- Scan Configuration ---
N_POINTS=20
START_DIST=2.0
END_DIST=5.0

# --- Physical Constants (Based on your final, successful slab relaxation) ---
CELL_Z=25.0
# This Z_AG_TOP value is calculated from your provided final slab coordinates.
Z_AG_TOP=0.344028

# --- Main Script Logic ---
echo "======================================================"
echo "Starting a fresh Potential Energy Surface (PES) scan."
echo "Total points to calculate: $N_POINTS"
echo "K-point mesh: 2x2x1"
echo "======================================================"

# Loop from 0 to N_POINTS-1
for i in $(seq 0 $((N_POINTS-1)))
do
    # Calculate the current H2 distance from the surface in Angstrom.
    CURRENT_DIST=$(echo "scale=10; $START_DIST + ($END_DIST - $START_DIST) * $i / ($N_POINTS - 1)" | bc -l)
    
    # Calculate the corresponding fractional Z coordinate for the H2 molecule's plane.
    Z_H2=$(echo "scale=10; $Z_AG_TOP + $CURRENT_DIST / $CELL_Z" | bc -l)

    # Create a clean filename label from the distance (e.g., 2.15 -> 2p15).
    DIST_LABEL=$(printf "%.2f" $CURRENT_DIST | tr '.' 'p')
    
    PREFIX="ag_h2_${DIST_LABEL}A_2k"
    INPUT_FILE="${PREFIX}.relax.in"
    OUTPUT_FILE="${PREFIX}.relax.out"

    echo ""
    echo "--> Step $((i+1))/$N_POINTS: Calculating for distance = $CURRENT_DIST A"

    # Create the input file for the current distance point.
    cat > $INPUT_FILE << EOF
&CONTROL
   calculation  = 'relax',
   prefix       = '$PREFIX',
   pseudo_dir   = './',
   outdir       = './',
   verbosity    = 'low'
/
&SYSTEM
   ibrav=0, nat=29, ntyp=2, ecutwfc=40.0, ecutrho=200.0,
   occupations='smearing', smearing='methfessel-paxton', degauss=0.02,
   vdw_corr = 'dft-d3'
/
&ELECTRONS
   conv_thr    = 1.0d-8,
   mixing_beta = 0.1
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
H      0.480000   0.500000   $Z_H2   1 1 0
H      0.520000   0.500000   $Z_H2   1 1 0
K_POINTS {automatic}
   2 2 1 0 0 0
EOF

    # Run the calculation and redirect output.
    echo "Running calculation... Output will be in $OUTPUT_FILE"
    $MPI_COMMAND $PW_COMMAND -i $INPUT_FILE > $OUTPUT_FILE

    # Check if the job finished successfully
    if grep -q "JOB DONE." "$OUTPUT_FILE"; then
        echo "Calculation for $CURRENT_DIST A finished successfully."
    else
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "ERROR: Calculation for $CURRENT_DIST A FAILED."
        echo "Please check the output file: $OUTPUT_FILE"
        echo "Aborting script."
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        exit 1 # Stop the script if any calculation fails
    fi
done

echo ""
echo "======================================================"
echo "All 20 calculations are finished successfully!"
echo "You can now extract the total energy from all .out files"
echo "and plot your high-resolution potential energy curve."
echo "======================================================"