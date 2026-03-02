# modules/amber_to_gromacs.py
from __future__ import annotations

import os
from pathlib import Path
from typing import List, Dict, Tuple


def _parse_mol2_atoms_coords(mol2_file: str) -> List[Dict]:
    atoms: List[Dict] = []
    in_atoms = False
    with open(mol2_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if s.startswith("@<TRIPOS>") and in_atoms:
                break
            if not in_atoms or not s:
                continue

            parts = line.split()
            # id name x y z type subst_id subst_name charge
            if len(parts) < 9:
                continue
            atom_id = int(parts[0])
            name = parts[1]
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            atype = parts[5]
            charge = float(parts[8])
            atoms.append({"id": atom_id, "name": name, "type": atype, "x": x, "y": y, "z": z, "charge": charge})
    if not atoms:
        raise ValueError("No atoms parsed from MOL2 (missing @<TRIPOS>ATOM?)")
    return atoms


def _parse_mol2_bonds(mol2_file: str) -> List[Tuple[int, int]]:
    bonds: List[Tuple[int, int]] = []
    in_bonds = False
    with open(mol2_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>BOND"):
                in_bonds = True
                continue
            if s.startswith("@<TRIPOS>") and in_bonds:
                break
            if not in_bonds or not s:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            a1 = int(parts[1])
            a2 = int(parts[2])
            bonds.append((a1, a2))
    return bonds


def convert_to_gromacs_best_effort(mol2_file: str, output_dir: str, residue_name: str) -> None:
    """
    Best-effort structure export:
    - Writes .gro coordinates from MOL2 coordinates (no RDKit re-embedding!)
    - Writes a minimal .itp with atom list + bonds ONLY.
      (This is structure-level, not a full parameterized topology.)

    Thesis-friendly honesty:
    - Full parameter export should be done via prmtop->gromacs (e.g., ParmEd),
      which you can add as an optional dependency.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    atoms = _parse_mol2_atoms_coords(mol2_file)
    bonds = _parse_mol2_bonds(mol2_file)

    gro_path = out / f"{residue_name}.gro"
    itp_path = out / f"{residue_name}.itp"

    # MOL2 coords are in Angstrom; GROMACS uses nm
    coords_nm = [(a["x"] / 10.0, a["y"] / 10.0, a["z"] / 10.0) for a in atoms]
    xs, ys, zs = zip(*coords_nm)
    padding = 1.2  # nm

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    box_x = (max_x - min_x) + 2 * padding
    box_y = (max_y - min_y) + 2 * padding
    box_z = (max_z - min_z) + 2 * padding

    shift_x = -min_x + padding
    shift_y = -min_y + padding
    shift_z = -min_z + padding

    # Write .gro
    with gro_path.open("w", encoding="utf-8") as gro:
        gro.write(f"{residue_name}\n")
        gro.write(f"{len(atoms):5d}\n")
        for i, (a, (x, y, z)) in enumerate(zip(atoms, coords_nm), start=1):
            atom_name = a["name"][:5]
            gro.write(
                f"{1:5d}{residue_name:<5}{atom_name:>5}{i:5d}"
                f"{x + shift_x:8.3f}{y + shift_y:8.3f}{z + shift_z:8.3f}\n"
            )
        gro.write(f"{box_x:10.5f}{box_y:10.5f}{box_z:10.5f}\n")

    # Write minimal .itp (structure only)
    with itp_path.open("w", encoding="utf-8") as itp:
        itp.write("; nsaa-paramgen: structure-only ITP (no parameters)\n")
        itp.write("; For full parameters: export from AMBER prmtop using ParmEd (recommended).\n\n")

        itp.write("[ moleculetype ]\n")
        itp.write(f"{residue_name} 3\n\n")

        itp.write("[ atoms ]\n")
        for i, a in enumerate(atoms, start=1):
            itp.write(
                f"{i:<5} {a['type']:<8} 1 {residue_name:<6} {a['name']:<6} {i:<5} {a['charge']:10.6f} 0.0000\n"
            )

        itp.write("\n[ bonds ]\n")
        for a1, a2 in bonds:
            itp.write(f"{a1:<5} {a2:<5} 1\n")

    print(f"GROMACS files written: {gro_path}, {itp_path}")