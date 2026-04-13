import sys
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# Mapping three-letter to one-letter amino acid codes
residue_map = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "HIE": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

# List all your PDB files in the order you want to process them
pdb_files = ["m134_final.pdb", "m162_final.pdb", "m164_final.pdb", "m181_final.pdb"

]

# Define atom pairs to measure distances (Chain A to Chain B)
atom_pairs = [
    (281, 152), (282, 152), (283, 152), (284, 152), (285, 152), (287, 152), (289, 152), (290, 152),
    (290, 153),
    (294, 157),
    (293, 159),
    (296, 160), (297, 160), (298, 160), (299, 160), (300, 160), (301, 160), (302, 160), (303, 160), 
    (304, 160), (305, 160),
    (309, 5), (311, 7),
    (298, 14), (299, 14), (300, 14), (301, 14), (302, 14), (303, 14), (304, 14),
    (313, 99), (314, 99), (315, 99), (316, 99), (317, 99), (319, 99), (320, 99), (321, 99), (322, 99),
    (313, 103), (314, 103), (315, 103), (316, 103), (317, 103), (319, 103), (320, 103), (321, 103), (322, 103),
    (314, 106), (315, 106), (316, 106), (317, 106), (319, 106), (320, 106), (321, 106), (322, 106), (323, 106),
    (319, 110), (320, 110), (321, 110), (322, 110), (323, 110), (325, 110),
    # Added glycine 286 pairs:
    (127, 286),
    (132, 286)
]

for i, pdb_file in enumerate(pdb_files, start=1):
    # Output file for each PDB
    output_file = f"distance_values_{i}.txt"
    orig_stdout = sys.stdout  # keep track of the original stdout

    with open(output_file, "w") as f:
        sys.stdout = f  # now any 'print' writes to the file

        # Load the PDB structure
        print(f"Processing {pdb_file}")
        u = mda.Universe(pdb_file)

        # Loop through trajectory frames
        for ts in u.trajectory:
            print(f"\n===== Frame {ts.frame} =====")

            for resA, resB in atom_pairs:
                # Use CA if residue in [5, 286, 313], else CB
                atomA = u.select_atoms(f"segid A and resid {resA} and name {'CA' if resA in [5, 286, 313] else 'CB'}")
                atomB = u.select_atoms(f"segid B and resid {resB} and name {'CA' if resB in [5, 286, 313] else 'CB'}")

                if len(atomA) == 0 or len(atomB) == 0:
                    print(f"WARNING: Missing atoms for resid {resA} (Chain A) or {resB} (Chain B) (frame {ts.frame}).")
                    continue

                # Calculate distance
                dist_matrix = distances.distance_array(atomA.positions, atomB.positions, box=u.dimensions)

                for ia, atomAi in enumerate(atomA):
                    for ib, atomBj in enumerate(atomB):
                        dist_val = dist_matrix[ia, ib]

                        # Convert three-letter to one-letter code
                        resA_one_letter = residue_map.get(atomAi.resname, "?")
                        resB_one_letter = residue_map.get(atomBj.resname, "?")

                        print(
                            f"Resid {atomAi.resid} ({resA_one_letter}-{atomAi.name}, Chain A) -> "
                            f"Resid {atomBj.resid} ({resB_one_letter}-{atomBj.name}, Chain B) : {dist_val:.3f} Å"
                        )

    # Restore the original stdout
    sys.stdout = orig_stdout
    print(f"Saved distances from {pdb_file} to {output_file}")
