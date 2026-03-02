# modules/itp_generator.py
from __future__ import annotations

from .utils import (
    get_atomtypes,
    get_mass,
    extract_charges,
    get_bonds_angles_dihedrals,
    get_impropers,
    get_pairs,
    fix_duplicates,
)


def generate_gromacs_itp(mol2_path: str, frcmod_paths, output_path: str, resname: str) -> None:
    """
    Generates a *structure-level* ITP (atoms/bonds/angles/dihedrals/pairs).
    It does not invent LJ/dihedral parameters; those should come from a GROMACS forcefield
    or from an Amber-export route (ParmEd recommended).
    """
    atoms, atomtypes = get_atomtypes(mol2_path)
    charges = extract_charges(mol2_path)
    masses = [get_mass(a["type"]) for a in atoms]

    bonds, angles, dihedrals = get_bonds_angles_dihedrals(mol2_path)
    impropers = get_impropers(mol2_path)
    pairs = get_pairs(dihedrals)

    # Fix duplicate issues
    angles = fix_duplicates(angles)
    dihedrals = fix_duplicates(dihedrals)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("; nsaa-paramgen structure-level itp\n")
        f.write("; atomtypes parameters must be provided by the selected forcefield\n\n")

        f.write("[ atomtypes ]\n")
        for at in atomtypes:
            f.write(f"{at}\n")

        f.write("\n[ moleculetype ]\n")
        f.write(f"{resname:<15} 3\n")

        f.write("\n[ atoms ]\n")
        for i, atom in enumerate(atoms, 1):
            f.write(
                f"{i:>5} {atom['type']:<8} 1 {resname:<6} {atom['name']:<6} {i:>5} "
                f"{charges[i-1]:>10.6f} {masses[i-1]:>10.5f}\n"
            )

        f.write("\n[ bonds ]\n")
        for b in bonds:
            f.write(f"{b}\n")

        f.write("\n[ pairs ]\n")
        for p in pairs:
            f.write(f"{p}\n")

        f.write("\n[ angles ]\n")
        for a in angles:
            f.write(f"{a}\n")

        f.write("\n[ dihedrals ]\n")
        for d in dihedrals:
            f.write(f"{d}\n")

        f.write("\n[ dihedrals ] ; impropers\n")
        for imp in impropers:
            f.write(f"{imp}\n")