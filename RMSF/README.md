# RMSF Analysis per Residue

This script calculates the root mean square fluctuation (RMSF) for C-alpha atoms in Chain A of a receptor trajectory using MDTraj.

## Description

The script:

- loads the structure file `step5_input.gro` and trajectory file `step7_filtered_skip50_Set_0_0_0.xtc`
- selects C-alpha atoms in Chain A for residues 1–342
- aligns the trajectory to the first frame using the selected C-alpha atoms
- computes the mean structure from the aligned trajectory
- calculates RMSF using the mean structure as the reference
- converts RMSF values from nanometers to Ångströms
- saves the per-residue RMSF values to a text file
- generates an RMSF plot

## Requirements
- Python
- MDTraj
- NumPy
- Matplotlib

## Input Files

- `step5_input.gro`
- `step7_filtered_skip50_Set_0_0_0.xtc`

## Output Files

- `rmsf_per_residue_edit.txt`

## Output Format

The output text file contains two columns:

1. Residue number
2. RMSF value in Å

Example:

```text
Residue_Number RMSF(Å)
1 0.84215
2 0.79103
3 0.76544
