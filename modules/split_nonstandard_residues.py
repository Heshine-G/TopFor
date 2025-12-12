# === FILE: split_nonstandard_residues.py ===
import os
import sys
from Bio.PDB import PDBParser, Select, PDBIO

# List of 20 standard amino acids
standard_residues = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

class NonStandardSelect(Select):
    def __init__(self, target_residue):
        self.target_residue = target_residue

    def accept_residue(self, residue):
        return residue == self.target_residue

def extract_nonstandard_residues(input_pdb, output_dir):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("peptide", input_pdb)

    os.makedirs(output_dir, exist_ok=True)
    io = PDBIO()

    extracted = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname not in standard_residues:
                    res_id = f"{resname}_{chain.id}_{residue.id[1]}"
                    output_path = os.path.join(output_dir, f"{res_id}.pdb")
                    io.set_structure(structure)
                    io.save(output_path, select=NonStandardSelect(residue))
                    extracted.append(output_path)

    return extracted

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python split_nonstandard_residues.py <input.pdb> <output_dir>")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_dir = sys.argv[2]

    paths = extract_nonstandard_residues(input_pdb, output_dir)
    if paths:
        print("Extracted non-standard residues:")
        for p in paths:
            print(f" - {p}")
    else:
        print("No non-standard residues found.")
