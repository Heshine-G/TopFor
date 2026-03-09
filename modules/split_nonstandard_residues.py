from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List

STANDARD_RESIDUES = {
    "ALA","ARG","ASN","ASP","CYS",
    "GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL",
}

HIS_VARIANTS = {"HID","HIE","HIP"}


def _clean_resname(name: str) -> str:
    """
    Convert OGU1 → OGU
    Convert SER7 → SER
    """
    name = name.split(".")[0]
    name = re.sub(r"\d+$", "", name)
    name = name.split("_")[0]
    return name.upper()


def _parse_sections(mol2_path: str):

    sections = {}
    current = None

    with open(mol2_path, "r", encoding="utf-8", errors="replace") as f:

        for line in f:

            if line.startswith("@<TRIPOS>"):
                current = line.strip()
                sections[current] = []
                continue

            if current:
                sections[current].append(line.rstrip())

    return sections


def _parse_atoms(lines: List[str]):

    atoms = []

    for line in lines:

        if not line.strip():
            continue

        p = line.split()

        if len(p) < 8:
            continue

        atoms.append(
            {
                "id": int(p[0]),
                "name": p[1],
                "x": p[2],
                "y": p[3],
                "z": p[4],
                "type": p[5],
                "subst_id": int(p[6]),
                "subst_name": p[7],
                "charge": p[8] if len(p) > 8 else "0.0",
            }
        )

    return atoms


def _parse_bonds(lines: List[str]):

    bonds = []

    for line in lines:

        if not line.strip():
            continue

        p = line.split()

        if len(p) < 4:
            continue

        bonds.append(
            {
                "id": int(p[0]),
                "a1": int(p[1]),
                "a2": int(p[2]),
                "type": p[3],
            }
        )

    return bonds


def extract_nonstandard_residues_from_mol2(input_mol2: str, output_dir: str):

    sections = _parse_sections(input_mol2)

    atoms = _parse_atoms(sections.get("@<TRIPOS>ATOM", []))
    bonds = _parse_bonds(sections.get("@<TRIPOS>BOND", []))

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    std = STANDARD_RESIDUES | HIS_VARIANTS

    residues: Dict[int, List[dict]] = {}

    for atom in atoms:
        residues.setdefault(atom["subst_id"], []).append(atom)

    extracted = []
    seen_residue_types = set()

    for subst_id, residue_atoms in residues.items():

        raw_name = residue_atoms[0]["subst_name"]

        resname = _clean_resname(raw_name)

        # skip standard residues
        if resname in std:
            continue

        # skip duplicates
        if resname in seen_residue_types:
            continue

        seen_residue_types.add(resname)

        output_file = out / f"{resname}.mol2"

        atom_ids = {a["id"] for a in residue_atoms}

        residue_bonds = [
            b for b in bonds
            if b["a1"] in atom_ids and b["a2"] in atom_ids
        ]

        id_map = {old: new for new, old in enumerate(atom_ids, start=1)}

        with open(output_file, "w") as f:

            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{resname}\n")
            f.write(f"{len(atom_ids)} {len(residue_bonds)} 1\n")
            f.write("SMALL\n")
            f.write("USER_CHARGES\n\n")

            f.write("@<TRIPOS>ATOM\n")

            for atom in residue_atoms:

                new_id = id_map[atom["id"]]

                f.write(
                    f"{new_id:>6} {atom['name']:<6} "
                    f"{float(atom['x']):>10.4f} {float(atom['y']):>10.4f} {float(atom['z']):>10.4f} "
                    f"{atom['type']:<6} 1 {resname:<6} {float(atom['charge']):>10.6f}\n"
                )

            f.write("@<TRIPOS>BOND\n")

            for i, bond in enumerate(residue_bonds, 1):

                f.write(
                    f"{i:>6} {id_map[bond['a1']]:>4} {id_map[bond['a2']]:>4} {bond['type']}\n"
                )

            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write(f"1 {resname} 1\n")

        extracted.append(str(output_file))

    return extracted


def extract_nonstandard_residues(input_path: str, output_dir: str):

    suffix = Path(input_path).suffix.lower()

    if suffix == ".mol2":
        return extract_nonstandard_residues_from_mol2(input_path, output_dir)

    raise ValueError("Peptide splitting currently supports only MOL2 input.")


def main():

    if len(sys.argv) != 3:
        print("Usage: split_nonstandard_residues.py peptide.mol2 output_dir")
        sys.exit(1)

    files = extract_nonstandard_residues(sys.argv[1], sys.argv[2])

    print(f"Generated {len(files)} unique NSAA residues")


if __name__ == "__main__":
    main()