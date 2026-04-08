import os
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# List to hold all loaded trajectories, RMSD data, and time data
all_trajs = []
all_rmsd = []
all_time_ns = []

# Colors for each set
colors = plt.cm.viridis(np.linspace(0, 1, 9))  # Using viridis colormap for distinct colors

# Create the plot figure
plt.figure(figsize=(10, 6))

# Loop over the sets (from Set_0_0_0 to Set_0_0_8)
for i in range(9):
    # Generate the directory and file names
    traj_file = f'./Set_0_0_{i}/step7_filtered_skip50_Set_0_0_{i}.xtc'
    top_file = f'./Set_0_0_{i}/step5_input.gro'

    # Load the trajectory and topology files
    print(f"Loading trajectory and topology files for Set_0_0_{i}...")
    traj = md.load(traj_file, top=top_file)

    # Apply PBC corrections by wrapping the molecules inside the unit cell
    print("Applying periodic boundary conditions (PBC) corrections...")
    traj = traj.image_molecules()

    # Select the C-alpha atoms of the protein for residues 1 to 342
    c_alpha_atoms = traj.topology.select('name CA and resid 0 to 341')

    # Align the trajectory to the first frame using the C-alpha atoms of residues 1 to 330
    print("Aligning the trajectory to the first frame using C-alpha atoms of residues 1 to 330...")
    traj.superpose(traj, 0, atom_indices=c_alpha_atoms)

    # Calculate RMSD relative to the first frame using the selected C-alpha atoms
    print("Calculating RMSD relative to the first frame...")
    rmsd = md.rmsd(traj, traj, frame=0, atom_indices=c_alpha_atoms)

    # Convert RMSD from nm to Å
    rmsd_angstrom = rmsd * 10

    # Convert time from ps to ns
    time_ns = traj.time / 1000

    # Store the trajectory, RMSD, and time values
    all_trajs.append(traj)
    all_rmsd.append(rmsd_angstrom)
    all_time_ns.append(time_ns)

    # Plot RMSD for the current trajectory with a distinct color
    plt.plot(time_ns, rmsd_angstrom, color=colors[i], label=f'Set_0_0_{i}')

# Save all RMSD values (in Å) to a text file
print("Saving all RMSD values (in Å) to rmsd_values_ca_aligned_res1_330_all_trajectories_angstrom.txt...")
all_rmsd_concat = np.concatenate(all_rmsd)
np.savetxt('rmsd_values_ca_aligned_res1_330_all_trajectories_angstrom.txt', all_rmsd_concat)

# Finalize the plot
print("Plotting RMSD over time (in ns) for all trajectories...")
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (Å)')

plt.ylim(0, 10)  # Set y-axis limit to 10 Å
plt.grid(False)  # Remove grid from the plot
plt.legend()  # This now correctly captures the labels from the loop
plt.tight_layout()
plt.savefig('rmsd_plot_ca_aligned_res1_330_all_trajectories_angstrom.png')
plt.show()

print("RMSD analysis complete. Results saved to rmsd_values_ca_aligned_res1_330_all_trajectories_angstrom_new.txt and rmsd_plot_ca_aligned_res1_330_all_trajectories_angstrom_new.png.")

