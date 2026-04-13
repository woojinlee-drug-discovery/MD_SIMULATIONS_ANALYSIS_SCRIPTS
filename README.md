# MD Simulation Scripts

This repository contains molecular dynamics analysis scripts used for GPCR simulation workflows.

## Included analysis categories

- C-alpha distance analysis
- RMSD analysis
- Ligand RMSD analysis
- Center-of-mass (COM) distance analysis
- RMSF analysis
- Ion interaction analysis
- 3D ion distribution visualization
- Crosslinking pair distance analysis

## Repository structure

Scripts are organized by analysis type.  
See the README file inside each folder for script-specific details, inputs, outputs, and usage.

## Common requirements

Depending on the script, the following Python packages may be required:

- MDAnalysis
- MDTraj
- NumPy
- Matplotlib

## Notes

- These scripts were developed for practical post-processing of MD simulations.
- Input filenames and directory layouts may reflect the original simulation environment.
- Trajectory and topology files are not included in this repository.
