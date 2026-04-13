# Inter-Chain Residue Distance Analysis from PDB Structures

This script measures selected inter-chain residue distances from a set of PDB structures using MDAnalysis.

## Description

The script:

- loads multiple PDB files in a specified order
- processes each structure individually
- measures distances between predefined residue pairs across **Chain A** and **Chain B**
- selects either `CA` or `CB` atoms depending on the residue
- converts residue names from three-letter codes to one-letter codes for reporting
- writes all measured distances to a separate text file for each input PDB

## Input Files

The script expects the following PDB files:

- `m134_final.pdb`
- `m162_final.pdb`
- `m164_final.pdb`
- `m181_final.pdb`

## Output Files

The script generates one output text file per PDB:

- `distance_values_1.txt`
- `distance_values_2.txt`
- `distance_values_3.txt`
- `distance_values_4.txt`

Each output file contains the distance measurements for the corresponding input PDB structure.

## What the Script Measures

The script calculates distances between selected residue pairs defined in `atom_pairs`.

- Residues are measured across:
  - `segid A`
  - `segid B`

- Atom selection rule:
  - use `CA` for residues `5`, `286`, and `313`
  - use `CB` for all other residues

## Output Format

Each output file includes:

- the name of the processed PDB file
- frame number
- residue identifiers
- one-letter amino acid codes
- atom names
- chain information
- measured distance in Å

Example output:

```text
Processing m134_final.pdb

===== Frame 0 =====
Resid 281 (X-CB, Chain A) -> Resid 152 (Y-CB, Chain B) : 8.532 Å
Resid 282 (X-CB, Chain A) -> Resid 152 (Y-CB, Chain B) : 7.941 Å
