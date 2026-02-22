# modules/remove.py
from __future__ import annotations

from pathlib import Path
from collections import Counter
from typing import List, Tuple, Optional


def _parse_mol2_atoms(mol2_path: str) -> List[Tuple[str, str, str]]:
    atoms: List[Tuple[str, str, str]] = []
    in_atoms = False
    with open(mol2_path, "r", encoding="utf-8", errors="replace") as f:
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
            if len(parts) < 8:
                continue
            atom_name = parts[1]
            atom_type = parts[5]
            subst_name = parts[7]
            atoms.append((atom_name, subst_name, atom_type))
    return atoms


def process_mol2_file(
    file_path: str,
    output_file: str,
    *,
    central_resname: Optional[str] = None,
    cap_resnames: Tuple[str, ...] = ("ACE", "NME"),
) -> str:
    atoms = _parse_mol2_atoms(file_path)
    if not atoms:
        raise ValueError(f"No MOL2 ATOM records parsed from: {file_path}")

    cap_set = {c.upper() for c in cap_resnames}

    subst_counts = Counter([subst for _, subst, _ in atoms if subst.upper() not in cap_set])
    if central_resname is None:
        if not subst_counts:
            raise ValueError("Could not infer central residue name (MOL2 seems to contain only caps).")
        central_resname = subst_counts.most_common(1)[0][0]

    omit_atom_names: List[str] = []
    kept = 0
    for atom_name, subst_name, _ in atoms:
        if subst_name.upper() in cap_set:
            omit_atom_names.append(atom_name)
        else:
            kept += 1

    if kept < 6:
        raise ValueError(
            f"After omitting caps {cap_resnames}, too few atoms remain ({kept}). "
            f"Inferred central residue: {central_resname}"
        )

    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", encoding="utf-8") as fh:
        fh.write("HEAD_NAME N\n")
        fh.write("TAIL_NAME C\n")
        fh.write("MAIN_CHAIN CA\n")
        for name in omit_atom_names:
            fh.write(f"OMIT_NAME {name}\n")
        fh.write("PRE_HEAD_TYPE C\n")
        fh.write("POST_TAIL_TYPE N\n")
        fh.write("CHARGE 0.0\n")

    return str(out)
