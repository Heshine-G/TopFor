# modules/split_nonstandard_residues.py
from __future__ import annotations

import os
import sys
from pathlib import Path
from Bio.PDB import PDBParser, Select, PDBIO


STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

# Optionally treat common histidine variants as standard if present in PDBs
HIS_VARIANTS = {"HID", "HIE", "HIP"}


class NonStandardSelect(Select):
    def __init__(self, target_residue):
        self.target_residue = target_residue

    def accept_residue(self, residue):
        return residue == self.target_residue


def extract_nonstandard_residues(input_pdb: str, output_dir: str, *, include_his_variants: bool = True) -> list[str]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("peptide", input_pdb)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    extracted: list[str] = []

    std = set(STANDARD_RESIDUES)
    if include_his_variants:
        std |= HIS_VARIANTS

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname and resname not in std:
                    res_id = f"{resname}_{chain.id}_{residue.id[1]}"
                    output_path = out / f"{res_id}.pdb"
                    io.set_structure(structure)
                    io.save(str(output_path), select=NonStandardSelect(residue))
                    extracted.append(str(output_path))

    return extracted


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: python split_nonstandard_residues.py <input.pdb> <output_dir>")
        return 1

    input_pdb = sys.argv[1]
    output_dir = sys.argv[2]

    paths = extract_nonstandard_residues(input_pdb, output_dir)
    if paths:
        print("Extracted non-standard residues:")
        for p in paths:
            print(f" - {p}")
    else:
        print("No non-standard residues found.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())