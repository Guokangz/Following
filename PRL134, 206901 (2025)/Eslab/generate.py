import numpy as np

# ================= INPUT DATA =================
# 1. Final relaxed CELL_PARAMETERS (in Angstrom) from your .out file
#    This defines the a, b, c vectors of our crystal box.
cell_vectors = np.array([
    [8.850, 0.000, 0.000],
    [-4.425, 7.665, 0.000],
    [0.000, 0.000, 25.000]
])

# 2. Final relaxed ATOMIC_POSITIONS (in crystal/fractional coordinates)
#    This is the block of 27 atoms you provided.
fractional_coords = [
    # Layer 3 (Bottom, FIXED)
    ['Ag', 0.2222220000, 0.1111110000, 0.2000000000],
    ['Ag', 0.5555560000, 0.1111110000, 0.2000000000],
    ['Ag', 0.8888890000, 0.1111110000, 0.2000000000],
    ['Ag', 0.2222220000, 0.4444440000, 0.2000000000],
    ['Ag', 0.5555560000, 0.4444440000, 0.2000000000],
    ['Ag', 0.8888890000, 0.4444440000, 0.2000000000],
    ['Ag', 0.2222220000, 0.7777780000, 0.2000000000],
    ['Ag', 0.5555560000, 0.7777780000, 0.2000000000],
    ['Ag', 0.8888890000, 0.7777780000, 0.2000000000],
    # Layer 2 (Middle, FIXED)
    ['Ag', 0.1111110000, 0.2222220000, 0.2500000000],
    ['Ag', 0.4444440000, 0.2222220000, 0.2500000000],
    ['Ag', 0.7777780000, 0.2222220000, 0.2500000000],
    ['Ag', 0.1111110000, 0.5555560000, 0.2500000000],
    ['Ag', 0.4444440000, 0.5555560000, 0.2500000000],
    ['Ag', 0.7777780000, 0.5555560000, 0.2500000000],
    ['Ag', 0.1111110000, 0.8888890000, 0.2500000000],
    ['Ag', 0.4444440000, 0.8888890000, 0.2500000000],
    ['Ag', 0.7777780000, 0.8888890000, 0.2500000000],
    # Layer 1 (Top, RELAXED)
    ['Ag', 0.0000029568, 0.0000059136, 0.3440280909],
    ['Ag', 0.3333370134, 0.0000060547, 0.3440276273],
    ['Ag', 0.6666690413, 0.0000060547, 0.3440276273],
    ['Ag', 0.0000030712, 0.3333395735, 0.3440283118],
    ['Ag', 0.3333365023, 0.3333395735, 0.3440283118],
    ['Ag', 0.6666702255, 0.3333394510, 0.3440279663],
    ['Ag', 0.0000043613, 0.6666747760, 0.3440281131],
    ['Ag', 0.3333371438, 0.6666752877, 0.3440280446],
    ['Ag', 0.6666704147, 0.6666747760, 0.3440281131]
]
# ===============================================

# --- SCRIPT LOGIC (No need to change anything below) ---

output_filename = "ag_slab_final_clean.xyz"
num_atoms = len(fractional_coords)

with open(output_filename, 'w') as f:
    # Write the header of the XYZ file
    f.write(f"{num_atoms}\n")
    f.write("Final relaxed structure from QE, converted to Cartesian coordinates\n")

    # Loop through each atom, convert coordinates, and write to file
    for atom in fractional_coords:
        symbol = atom[0]
        frac_vec = np.array(atom[1:4])
        
        # The core conversion: Position = x_f * a + y_f * b + z_f * c
        cart_vec = np.dot(frac_vec, cell_vectors)
        
        x, y, z = cart_vec[0], cart_vec[1], cart_vec[2]
        
        # Write the line in XYZ format
        f.write(f"  {symbol:<4s} {x:15.8f} {y:15.8f} {z:15.8f}\n")

print(f"SUCCESS: A clean XYZ file named '{output_filename}' has been created with exactly {num_atoms} atoms.")