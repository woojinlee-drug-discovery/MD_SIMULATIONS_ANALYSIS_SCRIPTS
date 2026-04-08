## rmsd_all_trajectories_analysis.py

This script calculates C-alpha RMSD across nine molecular dynamics trajectories using MDTraj and generates a combined RMSD plot.

### Description

The script:

- loads trajectory and topology files from `Set_0_0_0` to `Set_0_0_8`
- applies periodic boundary condition (PBC) corrections with `image_molecules()`
- selects C-alpha atoms for residues 1–342
- aligns each trajectory to its first frame
- calculates RMSD relative to the first frame
- converts RMSD values from nanometers to Ångströms
- plots all trajectories on a single RMSD vs time graph
- saves concatenated RMSD values to a text file
- saves the combined RMSD plot as a PNG image

### Input Files

For each set:
- `./Set_0_0_i/step7_filtered_skip50_Set_0_0_i.xtc`
- `./Set_0_0_i/step5_input.gro`

where `i = 0, 1, 2, ..., 8`

### Output Files

- `rmsd_values_ca_aligned_res1_330_all_trajectories_angstrom.txt`
- `rmsd_plot_ca_aligned_res1_330_all_trajectories_angstrom.png`

### Requirements

- Python
- MDTraj
- NumPy
- Matplotlib

### Usage

```bash
python rmsd_all_trajectories_analysis.py
