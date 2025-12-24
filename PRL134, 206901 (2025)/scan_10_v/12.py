import matplotlib.pyplot as plt
import numpy as np

# --- Input Your Data Here ---
# This is the data you extracted from your calculations.

# Distances in Angstrom
distances = np.array([
    2.00, 2.33, 2.67, 3.00, 3.33, 3.67, 4.00, 4.33, 4.67, 5.00
])

# Corresponding adsorption energies in meV
# (I have calculated these for you based on your provided total energies)
energies_meV = np.array([
    202.3, 33.0, -39.7, -66.3, -68.1, -58.3, -46.7, -35.7, -26.6, -19.8
])

# --- Plotting Logic (No need to change anything below) ---

# Create a new figure and axis for the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the calculated data points as blue circles with a connecting line
ax.plot(distances, energies_meV, marker='o', linestyle='-', color='b', label='Your PBE+D3 Calculation')

# Find and highlight the minimum energy point (the most stable adsorption site)
min_energy_index = np.argmin(energies_meV)
min_dist = distances[min_energy_index]
min_energy = energies_meV[min_energy_index]

# Add a vertical line and annotation for the minimum point
ax.axvline(x=min_dist, color='r', linestyle='--', linewidth=1, label=f'Equilibrium: {min_dist:.2f} Å')
ax.axhline(y=min_energy, color='r', linestyle='--', linewidth=1)
ax.annotate(f'E_ads = {min_energy:.1f} meV',
            xy=(min_dist, min_energy),
            xytext=(min_dist + 0.3, min_energy + 15),
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
            fontsize=12)

# Add a horizontal line at y=0 to represent infinite separation
ax.axhline(y=0, color='gray', linestyle='-', linewidth=1)

# --- Formatting and Labels ---
ax.set_xlabel('Surface Distance (Å)', fontsize=14)
ax.set_ylabel('Adsorption Energy (meV)', fontsize=14)
ax.set_title('Potential Energy Curve for H₂ on Ag(111)', fontsize=16)

# Add grid lines for better readability
ax.grid(True)

# Add a legend to identify the data series
ax.legend(fontsize=12)

# Set the limits of the x and y axes to match the article's plot for easy comparison
ax.set_xlim(2.0, 5.0)
# Adjust ylim to nicely fit the data
ax.set_ylim(np.min(energies_meV) - 20, np.max(energies_meV) + 20)


# --- Display and Save the Plot ---

# Save the figure to a file
output_filename = "H2_on_Ag111_PES.png"
plt.savefig(output_filename, dpi=300)
print(f"Plot saved successfully as '{output_filename}'")

# Show the plot in a new window
plt.show()