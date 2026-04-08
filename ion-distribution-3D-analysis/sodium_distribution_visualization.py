import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------
# 1. Load trajectory
# ----------------------
u = mda.Universe("step5_input.gro", "step7_filtered_skip50_Set_0_0_1.xtc")

# Select sodium atoms and D72 sidechain oxygens
sodium = u.select_atoms("resname SOD")  # or 'name NA'
d72_atoms = u.select_atoms("resid 72 and (name OD1 or name OD2)")

# Store sodium and D72 coordinates
all_na_coords = []
all_d72_coords = []

for ts in u.trajectory:
    all_na_coords.append(sodium.positions.copy())
    
    # Compute the average position of OD1 and OD2 to represent D72
    d72_center = d72_atoms.positions.mean(axis=0)
    all_d72_coords.append(d72_center)
    
    print(f"Frame {ts.frame}: {len(sodium)} Na⁺ atoms")

# Combine all coordinates across frames
all_na_coords = np.vstack(all_na_coords)
all_d72_coords = np.vstack(all_d72_coords)
print(f"\nTotal sodium coordinates accumulated: {all_na_coords.shape[0]}")
print(f"Total D72 center positions accumulated: {all_d72_coords.shape[0]}")

# ----------------------
# 2. Center everything around sodium cloud mean
# ----------------------
na_mean = all_na_coords.mean(axis=0)
na_centered = all_na_coords - na_mean
d72_centered = all_d72_coords - na_mean

# ----------------------
# 3. Plotting
# ----------------------
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
#ax.set_title("Na⁺ Distribution and D72 (Centered)", fontsize=14)

# Plot Na⁺ positions
ax.scatter(na_centered[:, 0], na_centered[:, 1], na_centered[:, 2],
           c='dodgerblue', s=10, alpha=0.4, label='Na⁺')

# Plot D72 as a single red dot per frame
ax.scatter(d72_centered[:, 0], d72_centered[:, 1], d72_centered[:, 2],
           c='red', s=30, alpha=0.8, label=r'D72$^{2.50}$')

# Axis labels
ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")
ax.legend()
ax.grid(True)

# Ensure equal aspect ratio
def set_axes_equal(ax):
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d()
    ])
    spans = limits[:, 1] - limits[:, 0]
    centers = np.mean(limits, axis=1)
    radius = 0.5 * max(spans)
    for center, axis in zip(centers, [ax.set_xlim3d, ax.set_ylim3d, ax.set_zlim3d]):
        axis(center - radius, center + radius)

set_axes_equal(ax)

# Set diagonal view between x and y axes
ax.view_init(elev=0.5, azim=45)

# Save the figure
plt.tight_layout()
output_filename = "sodium_distribution_with_D72.png"
plt.savefig(output_filename, dpi=300)
print(f"\n3D plot saved as: {output_filename}")

plt.show()
