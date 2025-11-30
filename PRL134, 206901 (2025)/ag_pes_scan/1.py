
import matplotlib.pyplot as plt
import numpy as np
# If you want a smoother curve, you can use scipy for interpolation
# from scipy.interpolate import make_interp_spline

# --- Input Your Final 20-Point Data Here ---
distances = np.array([
    2.00, 2.16, 2.32, 2.47, 2.63, 2.79, 2.95, 3.11, 3.26, 3.42, 
    3.58, 3.74, 3.89, 4.05, 4.21, 4.37, 4.53, 4.68, 4.84, 5.00
])

energies_meV = np.array([
    191.9, 117.2, 55.8, 10.6, -20.7, -40.6, -51.6, -56.4, -57.4, -56.1, 
    -53.1, -49.1, -44.4, -39.6, -34.8, -30.4, -26.4, -23.0, -19.9, -17.2
])

# Convert meV to eV for plotting
energies_eV = energies_meV / 1000.0

# --- Optional: Create a smooth curve for plotting ---
# x_smooth = np.linspace(distances.min(), distances.max(), 300)
# spl = make_interp_spline(distances, energies_eV, k=3)  # k=3 for cubic spline
# y_smooth = spl(x_smooth)

# --- Plotting Logic ---
fig, ax = plt.subplots(figsize=(10, 7))

# Plot the smooth curve
# ax.plot(x_smooth, y_smooth, linestyle='-', color='orangered', label='Your PBE+D3 (2k, 20 pts)')
# Plot the original data points
ax.plot(distances, energies_eV, marker='o', linestyle='-', color='orangered', label='Your PBE+D3 (2k, 20 pts)')


# Find and highlight the minimum energy point
min_energy_index = np.argmin(energies_eV)
min_dist = distances[min_energy_index]
min_energy = energies_eV[min_energy_index]

ax.annotate(f'Equilibrium\nDist: {min_dist:.2f} Å\nE_ads: {min_energy*1000:.1f} meV',
            xy=(min_dist, min_energy),
            xytext=(min_dist + 0.5, min_energy + 0.04),
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
            fontsize=12, ha='center')

# --- Formatting and Labels ---
ax.axhline(y=0, color='gray', linestyle='-', linewidth=1)
ax.set_xlabel('Surface Distance (Å)', fontsize=14)
ax.set_ylabel('Adsorption Energy (eV)', fontsize=14)
ax.set_title('High-Resolution Potential Energy Curve for H₂ on Ag(111)', fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.grid(True)
ax.legend(fontsize=12)

# Set Axis Limits to match the paper's plot
ax.set_xlim(2.0, 5.0)
ax.set_ylim(-0.1, 0.04)

# --- Display and Save the Plot ---
output_filename = "H2_on_Ag111_PES_20pts.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Plot saved successfully as '{output_filename}'")

plt.show()