# modules/capping.py
from __future__ import annotations

import os
import sys
from pathlib import Path

# PyMOL is optional but required for this step
try:
    import pymol  # type: ignore
except Exception:
    pymol = None


def _safe_unique_rename(selection: str, prefix: str) -> None:
    """
    Rename atoms in 'selection' to a stable prefix-based scheme to reduce collisions.
    We keep names <= 4 chars where possible (GROMACS-friendly too).
    """
    # Collect atom names in selection
    names = pymol.cmd.get_model(selection).atom  # type: ignore[attr-defined]
    used = set()

    # First pass: rename by element-like first char, plus index
    for i, atom in enumerate(names, start=1):
        newname = f"{prefix}{i}"
        # keep <=4 chars if possible (e.g. AC1, AC2, NM1 ...)
        if len(newname) > 4:
            newname = newname[:4]
        # ensure uniqueness (very defensive)
        while newname in used:
            newname = (newname[:3] + str((i % 10))).ljust(4)[:4]
        used.add(newname)
        pymol.cmd.alter(f"{selection} and id {atom.id}", f"name='{newname}'")  # type: ignore[attr-defined]


def main() -> int:
    if pymol is None:
        print("ERROR: PyMOL is not available in this Python environment. Install/enable PyMOL to use capping.")
        return 2

    pymol.finish_launching(["pymol", "-cq"])

    if len(sys.argv) < 2:
        print("Usage: python capping.py <residue_folder>")
        return 1

    residue_folder = Path(sys.argv[1]).resolve()
    input_file = residue_folder / "residue.mol2"
    output_file = residue_folder / "residue_capped.mol2"

    if not input_file.exists():
        print(f"ERROR: {input_file} not found!")
        return 1

    pymol.cmd.reinitialize()  # type: ignore[attr-defined]
    pymol.cmd.load(str(input_file), "prot")  # type: ignore[attr-defined]

    # remove existing terminal markers/hydrogens then re-add later
    pymol.cmd.remove("hydro")  # type: ignore[attr-defined]
    pymol.cmd.remove("name OXT")  # type: ignore[attr-defined]

    # Attach ACE to N terminus if N exists
    pymol.cmd.select("pk1", "name N")  # type: ignore[attr-defined]
    if pymol.cmd.count_atoms("pk1") > 0:  # type: ignore[attr-defined]
        pymol.cmd.editor.attach_amino_acid("pk1", "ace")  # type: ignore[attr-defined]
        # Rename ACE atoms to AC*
        _safe_unique_rename("resn ACE", "AC")

    # Attach NME to C terminus if carbonyl C exists (not the ACE carbon)
    pymol.cmd.select("pk1", "name C and not resn ACE")  # type: ignore[attr-defined]
    if pymol.cmd.count_atoms("pk1") > 0:  # type: ignore[attr-defined]
        pymol.cmd.editor.attach_amino_acid("pk1", "nme")  # type: ignore[attr-defined]
        _safe_unique_rename("resn NME", "NM")

    pymol.cmd.h_add("prot")  # type: ignore[attr-defined]
    pymol.cmd.save(str(output_file), "prot")  # type: ignore[attr-defined]

    if output_file.exists() and output_file.stat().st_size > 0:
        print(f"Capping successful: {output_file}")
        return 0

    print(f"ERROR: PyMOL failed to create {output_file}!")
    return 3


if __name__ == "__main__":
    raise SystemExit(main())