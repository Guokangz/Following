import matplotlib.pyplot as plt
import numpy as np

# --- Input Your Final Data Here ---
# This is the data from your successful 2k scan.

# Distances in Angstrom
distances = np.array([
    2.00, 2.33, 2.67, 3.00, 3.33, 3.67, 4.00, 4.33, 4.67, 5.00
])

# Corresponding adsorption energies in eV 
# (I have converted your meV data to eV)
energies_eV = np.array([
    0.1919, 0.0500, -0.0261, -0.0538, -0.0571, -0.0510, -0.0412, -0.0314, -0.0233, -0.0172
])

# --- Plotting Logic (No need to change anything below) ---

# Create a new figure and axis for the plot
fig, ax = plt.subplots(figsize=(10, 7))

# Plot the calculated data points with a connecting line
# We'll use an orange color to distinguish from the paper's blue line
ax.plot(distances, energies_eV, marker='o', linestyle='-', color='orangered', label='Your PBE+D3 (2k) Calculation')

# Find and highlight the minimum energy point
min_energy_index = np.argmin(energies_eV)
min_dist = distances[min_energy_index]
min_energy = energies_eV[min_energy_index]

# Add annotation for the minimum point
ax.annotate(f'Equilibrium\nDist: {min_dist:.2f} Å\nE_ads: {min_energy*1000:.1f} meV',
            xy=(min_dist, min_energy),
            xytext=(min_dist + 0.4, min_energy + 0.03),
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
            fontsize=12, ha='center')

# Add a horizontal line at y=0
ax.axhline(y=0, color='gray', linestyle='-', linewidth=1)

# --- Formatting and Labels ---
ax.set_xlabel('Surface Distance (Å)', fontsize=14)
ax.set_ylabel('Adsorption Energy (eV)', fontsize=14)
ax.set_title('Potential Energy Curve for H₂ on Ag(111)', fontsize=16)

# Improve tick label size
ax.tick_params(axis='both', which='major', labelsize=12)

# Add grid lines
ax.grid(True)

# Add a legend
ax.legend(fontsize=12)

# --- Set Axis Limits (as per your request) ---
ax.set_xlim(2.0, 5.0)
ax.set_ylim(-0.1, 0.04) # Set Y-axis from -0.1 eV to +0.04 eV

# Invert the y-axis if you want deeper energies to be lower (optional)
# ax.invert_yaxis() # Uncomment this line if you want to flip the y-axis

# --- Display and Save the Plot ---
output_filename = "H2_on_Ag111_PES_Final.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Plot saved successfully as '{output_filename}'")

plt.show()