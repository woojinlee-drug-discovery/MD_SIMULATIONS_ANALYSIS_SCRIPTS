#!/usr/bin/env python3
"""
com_distance_mdtraj.py

Compute the center of mass distance between receptor (resid 1–434)
and beta-arrestin (resid 435–796) using MDTraj, and save the result
in Ångströms to ./com_distance/.
"""

import mdtraj as md
import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

# Load trajectory and topology
gro_file = "step5_input.gro"
xtc_file = "step7_filtered_skip50_Set_0_0_1.xtc"
traj = md.load(xtc_file, top=gro_file)

# Apply PBC correction
traj = traj.image_molecules()

# Select atom indices for receptor and arrestin
receptor_atoms = traj.topology.select("resid 0 to 433 and name CA")
arrestin_atoms = traj.topology.select("resid 434 to 795 and name CA")

# Atom masses from MDTraj topology
masses = np.array([atom.element.mass for atom in traj.topology.atoms])

# Center of mass calculation function
def compute_com(xyz, indices):
    sub_xyz = xyz[indices]
    sub_masses = masses[indices]
    com = np.average(sub_xyz, axis=0, weights=sub_masses)
    return com

# Compute COM distance for each frame
com_distances = []
for i in tqdm(range(traj.n_frames)):
    frame_xyz = traj.xyz[i]
    com_receptor = compute_com(frame_xyz, receptor_atoms)
    com_arrestin = compute_com(frame_xyz, arrestin_atoms)
    dist = np.linalg.norm(com_arrestin - com_receptor)
    com_distances.append(dist)

# Convert distances to Ångströms
com_distances = np.array(com_distances) * 10.0  # nm → Å
time_ns = np.arange(len(com_distances)) * 0.5  # Adjust timestep if needed

# Create output directory
output_dir = "com_distance"
os.makedirs(output_dir, exist_ok=True)

# Save results
output_txt = os.path.join(output_dir, "com_distance_receptor_arrestin.txt")
np.savetxt(output_txt,
           np.column_stack((time_ns, com_distances)),
           header="Time (ns)\tCOM Distance (Å)")

# Plot and save the image
plt.figure(figsize=(8, 5))
plt.plot(time_ns, com_distances, lw=2)
plt.xlabel("Time (ns)")
plt.ylabel("COM Distance (Å)")
plt.title("Center of Mass Distance: Receptor vs Arrestin")
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "com_distance_mdtraj_plot.png"), dpi=300)
plt.show()
