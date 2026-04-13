# ACHO Ligand RMSD Analysis

This script calculates the RMSD of the ACHO ligand over time for multiple MD trajectory folders using MDTraj.

## Description

The script:

- loops over trajectory folders from `Set_0_0_0` to `Set_0_0_9`
- loads the topology and trajectory files from each folder
- selects atoms belonging to the ligand with residue name `ACHO`
- extracts the ligand-only trajectory
- uses the first frame as the reference structure
- calculates the ligand RMSD over time
- saves RMSD values to a text file in each folder
- generates and saves an RMSD plot for each folder

## Input Files

For each folder `Set_0_0_i`, the script expects:

- `step5_input.gro`
- `step7_filtered_skip50_Set_0_0_i.xtc`

where `i = 0, 1, 2, ..., 9`

## Output Files

For each folder `Set_0_0_i`, the script generates:

- `acho_rmsd.txt`
- `acho_rmsd_plot.png`

## Output Format

### `acho_rmsd.txt`
The output text file contains RMSD values in nanometers.

Example:

```text
# RMSD (nm)
0.000000000000000000e+00
1.234567890000000012e-01
1.456789010000000034e-01
