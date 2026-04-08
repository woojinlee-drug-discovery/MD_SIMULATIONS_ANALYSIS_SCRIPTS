# MD Simulation Scripts

Molecular dynamics analysis scripts used for GPCR simulations, including C-alpha distance analysis with MDAnalysis.

## Script
- `calpha_distance_analysis.py`

## Description
This script calculates the C-alpha–C-alpha distance between residue 72 (D72) and residue 76 (K76) across multiple molecular dynamics trajectories using MDAnalysis.

## Input Files
- `step5_input.gro`
- `step7_filtered_skip50_Set_0_0_0.xtc`
- `step7_filtered_skip50_Set_0_0_1.xtc`
- `step7_filtered_skip50_Set_0_0_2.xtc`
- `step7_filtered_skip50_Set_0_0_3.xtc`
- `step7_filtered_skip50_Set_0_0_4.xtc`
- `step7_filtered_skip50_Set_0_0_5.xtc`
- `step7_filtered_skip50_Set_0_0_6.xtc`
- `step7_filtered_skip50_Set_0_0_7.xtc`
- `step7_filtered_skip50_Set_0_0_8.xtc`

## Output Files
- `CA72_CA76_distances_Set_0.txt` to `CA72_CA76_distances_Set_8.txt`
- `CA72_CA76_distances_combined.txt`
- `ca72_ca76_distance_combined.png`

## Requirements
- Python
- MDAnalysis
- NumPy
- Matplotlib

## Usage
```bash
python calpha_distance_analysis.py
