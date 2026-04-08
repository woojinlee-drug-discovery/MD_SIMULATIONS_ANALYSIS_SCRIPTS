import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
import os

# Define file names and parameters
trajectory_files = [f'step7_filtered_skip50_Set_0_0_{i}.xtc' for i in range(9)]
structure_file = 'step5_input.gro'

# Prepare to store combined distances and time
combined_distances_ca72_ca76 = []
combined_time_ns = []

# Loop over all trajectory files
for i, traj_file in enumerate(trajectory_files):
    if not os.path.exists(traj_file):
        print(f"Trajectory file {traj_file} not found. Skipping.")
        continue

    print(f"Processing {traj_file}...")

    # Load the universe for each trajectory
    u = mda.Universe(structure_file, traj_file)

    # Select C-alpha atoms of residues 72 (D72) and 76 (K76)
    ca_d72 = u.select_atoms('resid 72 and name CA')
    ca_k76 = u.select_atoms('resid 76 and name CA')

    # Prepare to store distances for the current trajectory
    distances_ca72_ca76 = []
    time_ns = []

    # Calculate distances between selected atoms over time
    for ts in u.trajectory:
        dist_matrix = distance_array(ca_d72.positions, ca_k76.positions, box=ts.dimensions)

        # Since these are single atoms, get the distance directly
        distance = dist_matrix[0, 0]

        distances_ca72_ca76.append(distance)
        time_ns.append(u.trajectory.time / 1000)  # Convert time from ps to ns

    # Append the current trajectory's data to the combined lists
    combined_distances_ca72_ca76.extend(distances_ca72_ca76)
    combined_time_ns.extend(time_ns)

    # Save the individual results for this trajectory
    np.savetxt(
        f"CA72_CA76_distances_Set_{i}.txt",
        distances_ca72_ca76,
        fmt="%0.4f",
        header=f"Distances between C-alpha atoms of residue 72 (D72) and 76 (K76) for trajectory Set {i}"
    )

# Convert combined lists to numpy arrays
combined_distances_ca72_ca76 = np.array(combined_distances_ca72_ca76)
combined_time_ns = np.array(combined_time_ns)

# Save the combined results to a single text file
np.savetxt(
    "CA72_CA76_distances_combined.txt",
    combined_distances_ca72_ca76,
    fmt="%0.4f",
    header="Distances between C-alpha atoms of residue 72 (D72) and residue 76 (K76) over all frames and trajectories"
)

# Plot the C-alpha distance over time for the combined data
plt.figure(figsize=(12, 8))
plt.plot(
    combined_time_ns,
    combined_distances_ca72_ca76,
    label='C-alpha Distance (Residue 72 - Residue 76)',
    color='green'
)

plt.xlabel('Time (ns)')
plt.ylabel('Distance (Å)')
plt.title('C-alpha Distance between Residue 72 (D72) and Residue 76 (K76) over Time (All Trajectories)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('ca72_ca76_distance_combined.png')
plt.show()

print("Distance analysis complete. "
      "Combined results saved to CA72_CA76_distances_combined.txt and ca72_ca76_distance_combined.png.")
