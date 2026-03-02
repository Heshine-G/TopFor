# modules/utils.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Set, Optional


@dataclass(frozen=True)
class AtomRec:
    atom_id: int
    name: str
    atype: str
    charge: float


def get_mass(atom_type: str) -> float:
    # Lightweight mass map (fallback)
    masses = {
        "H": 1.008, "HO": 1.008, "HC": 1.008, "H1": 1.008, "H2": 1.008, "H3": 1.008, "HN": 1.008,
        "C": 12.011, "CA": 12.011, "CT": 12.011, "CM": 12.011, "C2": 12.011, "C3": 12.011,
        "N": 14.007, "NA": 14.007, "N3": 14.007, "NB": 14.007,
        "O": 15.999, "O2": 15.999, "OH": 15.999, "OS": 15.999,
        "S": 32.06, "SH": 32.06,
    }
    return masses.get(atom_type, 12.011)


def _parse_mol2_atoms(mol2_path: str) -> List[AtomRec]:
    atoms: List[AtomRec] = []
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
            if len(parts) < 9:
                continue
            atom_id = int(parts[0])
            name = parts[1]
            atype = parts[5]
            charge = float(parts[-1])
            atoms.append(AtomRec(atom_id, name, atype, charge))
    if not atoms:
        raise ValueError(f"No atoms parsed from MOL2: {mol2_path}")
    return atoms


def _parse_mol2_bonds(mol2_path: str) -> List[Tuple[int, int]]:
    bonds: List[Tuple[int, int]] = []
    in_bonds = False
    with open(mol2_path, "r", encoding="utf-8", errors="replace") as f:
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


def extract_charges(mol2_path: str) -> List[float]:
    atoms = _parse_mol2_atoms(mol2_path)
    atoms_sorted = sorted(atoms, key=lambda a: a.atom_id)
    return [a.charge for a in atoms_sorted]


def get_atomtypes(mol2_path: str):
    """
    Returns:
      atoms: list of dicts with keys {id, name, type}
      atomtypes: unique list of types (as strings) for an [ atomtypes ] block placeholder
    """
    atoms = _parse_mol2_atoms(mol2_path)
    atoms_sorted = sorted(atoms, key=lambda a: a.atom_id)

    atom_dicts = [{"id": a.atom_id, "name": a.name, "type": a.atype} for a in atoms_sorted]
    unique_types = sorted({a.atype for a in atoms_sorted})
    atomtypes_block = [f"; {t} (parameters from forcefield)" for t in unique_types]
    return atom_dicts, atomtypes_block


def detect_formal_charge_from_mol2(mol2_path: str) -> Optional[int]:
    """
    Detect *formal* net charge from MOL2 using RDKit.
    Returns:
      int  -> detected net formal charge
      None -> RDKit not available OR RDKit failed to parse mol2

    Why Optional?
    - Your pipeline should still run on clusters where RDKit isn't installed.
    - In that case you MUST provide charge via JSON/CLI/default.
    """
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return None

    try:
        mol = Chem.MolFromMol2File(mol2_path, removeHs=False, sanitize=True)
        if mol is None:
            return None
        return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
    except Exception:
        return None


def renormalize_mol2_partial_charges_to_integer(mol2_path: str, target_integer_charge: int) -> bool:
    """
    Ensures the *sum of partial charges* in a MOL2 equals target_integer_charge.
    We distribute the correction uniformly across ALL atoms (simple + stable).
    Returns True if a correction was applied.
    """
    charges = extract_charges(mol2_path)
    current_sum = float(sum(charges))
    diff = float(target_integer_charge) - current_sum

    if abs(diff) <= 1e-6:
        return False

    per_atom = diff / float(len(charges))

    lines = open(mol2_path, "r", encoding="utf-8", errors="replace").read().splitlines()
    new_lines: List[str] = []

    in_atoms = False
    atom_seen = 0

    for line in lines:
        s = line.strip()
        if s.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            new_lines.append(line)
            continue
        if s.startswith("@<TRIPOS>") and in_atoms and not s.startswith("@<TRIPOS>ATOM"):
            in_atoms = False
            new_lines.append(line)
            continue

        if in_atoms and s:
            parts = line.split()
            if len(parts) >= 9:
                atom_seen += 1
                q = float(parts[-1]) + per_atom
                parts[-1] = f"{q:.6f}"
                new_lines.append(" ".join(parts))
                continue

        new_lines.append(line)

    # Defensive: only overwrite if we touched expected number of atoms
    if atom_seen != len(charges):
        raise RuntimeError(
            f"Charge renorm mismatch: saw {atom_seen} atom lines but parsed {len(charges)} charges."
        )

    with open(mol2_path, "w", encoding="utf-8") as f:
        f.write("\n".join(new_lines) + "\n")

    return True


def _build_graph(num_atoms: int, bonds: List[Tuple[int, int]]) -> Dict[int, Set[int]]:
    g: Dict[int, Set[int]] = {i: set() for i in range(1, num_atoms + 1)}
    for a, b in bonds:
        if a in g and b in g:
            g[a].add(b)
            g[b].add(a)
    return g


def get_bonds_angles_dihedrals(mol2_path: str):
    atoms = _parse_mol2_atoms(mol2_path)
    atoms_sorted = sorted(atoms, key=lambda a: a.atom_id)
    num_atoms = len(atoms_sorted)

    bonds = _parse_mol2_bonds(mol2_path)
    g = _build_graph(num_atoms, bonds)

    bond_lines = [f"{a:>5} {b:>5} 1" for a, b in sorted({tuple(sorted(x)) for x in bonds})]

    angles_set: Set[Tuple[int, int, int]] = set()
    for j in range(1, num_atoms + 1):
        neigh = sorted(g[j])
        for idx_i in range(len(neigh)):
            for idx_k in range(idx_i + 1, len(neigh)):
                i = neigh[idx_i]
                k = neigh[idx_k]
                angles_set.add((i, j, k) if i < k else (k, j, i))
    angle_lines = [f"{i:>5} {j:>5} {k:>5} 1" for (i, j, k) in sorted(angles_set)]

    dihed_set: Set[Tuple[int, int, int, int]] = set()
    for j, k in sorted({tuple(sorted(b)) for b in bonds}):
        for i in g[j]:
            if i == k:
                continue
            for l in g[k]:
                if l == j or l == i:
                    continue
                tup = (i, j, k, l)
                rev = (l, k, j, i)
                dihed_set.add(tup if tup < rev else rev)

    dihed_lines = [f"{i:>5} {j:>5} {k:>5} {l:>5} 1" for (i, j, k, l) in sorted(dihed_set)]
    return bond_lines, angle_lines, dihed_lines


def get_impropers(mol2_path: str):
    return []


def get_pairs(dihedrals: List[str]) -> List[str]:
    pairs: Set[Tuple[int, int]] = set()
    for d in dihedrals:
        parts = d.split()
        if len(parts) < 4:
            continue
        i = int(parts[0])
        l = int(parts[3])
        a, b = (i, l) if i < l else (l, i)
        pairs.add((a, b))
    return [f"{a:>5} {b:>5} 1" for (a, b) in sorted(pairs)]


def fix_duplicates(entries: List[str]) -> List[str]:
    seen = set()
    out = []
    for e in entries:
        if e not in seen:
            seen.add(e)
            out.append(e)
    return out