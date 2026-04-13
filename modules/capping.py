from __future__ import annotations

import json
import sys
from pathlib import Path

try:
    import pymol  # type: ignore
except Exception:
    pymol = None


MAX_ATOM_NAME_LEN = 4


def _all_atom_names(selection: str = "all") -> set[str]:
    atoms = pymol.cmd.get_model(selection).atom  
    return {str(atom.name).strip() for atom in atoms}


def _make_unique_name(base: str, used: set[str]) -> str:
    base = (base or "X")[:MAX_ATOM_NAME_LEN]
    if base not in used:
        return base

    for i in range(1, 10000):
        suffix = str(i)
        prefix = base[: MAX_ATOM_NAME_LEN - len(suffix)]
        candidate = f"{prefix}{suffix}"
        if candidate not in used:
            return candidate

    raise RuntimeError(f"Could not generate a unique atom name from base '{base}'")


def _safe_unique_rename(selection: str, prefix: str) -> None:
    atoms = pymol.cmd.get_model(selection).atom  
    if not atoms:
        return

    used = _all_atom_names("all")

    for i, atom in enumerate(atoms, start=1):
        old_name = str(atom.name).strip()
        used.discard(old_name)
        newname = _make_unique_name(f"{prefix}{i}", used)
        used.add(newname)
        pymol.cmd.alter(f"{selection} and id {atom.id}", f"name='{newname}'")  


def _normalize_name_arg(value: str | None) -> str | None:
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.upper() in {"NONE", "NULL", "0"}:
        return None
    return s


def main() -> int:
    if pymol is None:
        print("ERROR: PyMOL is not available in this Python environment. Install/enable PyMOL to use capping.")
        return 2

    pymol.finish_launching(["pymol", "-cq"])

    if len(sys.argv) < 2:
        print("Usage: python capping.py <residue_folder> [head_name|NONE] [tail_name|NONE]")
        return 1

    residue_folder = Path(sys.argv[1]).resolve()
    head_name = _normalize_name_arg(sys.argv[2] if len(sys.argv) >= 3 else "N")
    tail_name = _normalize_name_arg(sys.argv[3] if len(sys.argv) >= 4 else "C")

    input_file = residue_folder / "residue.mol2"
    output_file = residue_folder / "residue_capped.mol2"
    meta_file = residue_folder / "residue_capping_meta.json"

    if not input_file.exists():
        print(f"ERROR: {input_file} not found!")
        return 1

    pymol.cmd.reinitialize()  
    pymol.cmd.load(str(input_file), "prot")  

    pymol.cmd.remove("hydro")  
    pymol.cmd.remove("name OXT")  

    applied_caps: list[str] = []
    found_head = False
    found_tail = False

    if head_name:
        head_sel = f"name {head_name}"
        pymol.cmd.select("pk1", head_sel)  
        found_head = pymol.cmd.count_atoms("pk1") > 0  
        if found_head:
            pymol.cmd.editor.attach_amino_acid("pk1", "ace")  
            _safe_unique_rename("resn ACE", "AC")
            applied_caps.append("ACE")

    if tail_name:
        tail_sel = f"name {tail_name} and not resn ACE"
        pymol.cmd.select("pk1", tail_sel)  
        found_tail = pymol.cmd.count_atoms("pk1") > 0  
        if found_tail:
            pymol.cmd.editor.attach_amino_acid("pk1", "nme")  
            _safe_unique_rename("resn NME", "NM")
            applied_caps.append("NME")

    pymol.cmd.h_add("prot")  
    pymol.cmd.save(str(output_file), "prot")  

    meta = {
        "requested_head_name": head_name,
        "requested_tail_name": tail_name,
        "has_head": bool(found_head),
        "has_tail": bool(found_tail),
        "applied_caps": applied_caps,
    }
    meta_file.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    if output_file.exists() and output_file.stat().st_size > 0:
        print(f"Capping successful: {output_file}")
        print(json.dumps(meta))
        return 0

    print(f"ERROR: PyMOL failed to create {output_file}!")
    return 3


if __name__ == "__main__":
    raise SystemExit(main())