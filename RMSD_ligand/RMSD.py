import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os

# Loop over each Set_0_0_i folder
for i in range(10):
    folder = f'./Set_0_0_{i}'
    topology_file = os.path.join(folder, 'step5_input.gro')
    trajectory_file = os.path.join(folder, f'step7_filtered_skip50_Set_0_0_{i}.xtc')
    output_txt = os.path.join(folder, 'acho_rmsd.txt')
    output_plot = os.path.join(folder, 'acho_rmsd_plot.png')

    # Load trajectory
    traj = md.load(trajectory_file, top=topology_file)

    # Select ACHO ligand atoms
    acho_atoms = traj.topology.select('resname ACHO')
    acho_traj = traj.atom_slice(acho_atoms)

    # Reference = first frame
    reference_frame = acho_traj[0]

    # RMSD calculation
    rmsd = md.rmsd(acho_traj, reference_frame)

    # Save RMSD to text file
    np.savetxt(output_txt, rmsd, header='RMSD (nm)')

    # Plot
    time = traj.time / 1000  # ns
    plt.figure()
    plt.plot(time, rmsd)
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (nm)')
    plt.title(f'RMSD of ACHO Ligand Over Time (Set_0_0_{i})')
    plt.savefig(output_plot)
    plt.close()

print("RMSD calculation completed for all folders.")
