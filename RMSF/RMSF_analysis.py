import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Load the input structure and trajectory
traj = md.load('step7_filtered_skip50_Set_0_0_0.xtc', top='step5_input.gro')

# Select the alpha carbons in Chain A, residues 1-342
atom_indices = traj.topology.select('(chainid == 0) and (resid >= 0 and resid < 342) and (name CA)')

# Align the trajectory to the first frame using the selected atoms
traj.superpose(traj, 0, atom_indices=atom_indices)

# Calculate the mean structure (average structure)
mean_structure = traj.xyz.mean(axis=0)

# Create a new Trajectory object for the mean structure
mean_traj = md.Trajectory(mean_structure[np.newaxis, :, :], traj.topology)

# Calculate RMSF using the mean structure as the reference
rmsf = md.rmsf(traj, mean_traj, atom_indices=atom_indices)

# Convert RMSF to Ångströms (MDTraj outputs in nanometers by default, multiply by 10)
rmsf_angstroms = rmsf * 10

# Save RMSF values to a text file
residue_numbers = np.arange(1, len(rmsf_angstroms) + 1)
rmsf_data = np.column_stack((residue_numbers, rmsf_angstroms))

np.savetxt('rmsf_per_residue_edit.txt', rmsf_data, fmt='%d %.5f', header='Residue_Number RMSF(Å)')

# Plot the RMSF
plt.figure(figsize=(10, 6))
plt.plot(residue_numbers, rmsf_angstroms)
plt.xlabel('Residue number (Chain A)')
plt.ylabel('RMSF (Å)')
plt.title('RMSF per Residue for Chain A Receptor (1-342)')
plt.show()
