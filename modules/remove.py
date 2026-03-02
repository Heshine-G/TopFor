# modules/remove.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections import Counter, deque
from typing import Dict, List, Tuple, Optional, Set


@dataclass(frozen=True)
class Mol2Atom:
    atom_id: int
    name: str
    atom_type: str
    subst_name: str


def _parse_mol2_atoms(mol2_path: str) -> List[Mol2Atom]:
    atoms: List[Mol2Atom] = []
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
            atom_id = int(parts[0])
            atom_name = parts[1]
            atom_type = parts[5]
            subst_name = parts[7]
            atoms.append(Mol2Atom(atom_id, atom_name, atom_type, subst_name))
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
            # bond_id = int(parts[0])
            a1 = int(parts[1])
            a2 = int(parts[2])
            bonds.append((a1, a2))
    return bonds


def _shortest_path(graph: Dict[int, Set[int]], start: int, goal: int) -> Optional[List[int]]:
    if start == goal:
        return [start]
    q = deque([start])
    prev: Dict[int, int] = {}
    seen = {start}
    while q:
        v = q.popleft()
        for n in graph.get(v, set()):
            if n in seen:
                continue
            seen.add(n)
            prev[n] = v
            if n == goal:
                # reconstruct
                path = [goal]
                cur = goal
                while cur != start:
                    cur = prev[cur]
                    path.append(cur)
                path.reverse()
                return path
            q.append(n)
    return None


def process_mol2_file(
    file_path: str,
    output_file: str,
    *,
    head_name: str = "N",
    tail_name: str = "C",
    main_chain: Optional[List[str]] = None,
    charge: float = 0.0,
    central_resname: Optional[str] = None,
    cap_resnames: Tuple[str, ...] = ("ACE", "NME"),
    pre_head_type: str = "C",
    post_tail_type: str = "N",
    infer_mainchain_from_connectivity: bool = True,
) -> str:
    """
    Write an AmberTools 'mc' file for prepgen, omitting cap atoms and defining head/tail.

    Thesis-grade improvement:
    - If infer_mainchain_from_connectivity=True, MAIN_CHAIN atoms are inferred as the shortest
      bonded path between head and tail within the central residue (caps removed).
    - That automatically supports beta/gamma residues without hardcoding CB/CG/CE.
    """
    atoms = _parse_mol2_atoms(file_path)
    if not atoms:
        raise ValueError(f"No MOL2 ATOM records parsed from: {file_path}")

    bonds = _parse_mol2_bonds(file_path)
    cap_set = {c.upper() for c in cap_resnames}

    subst_counts = Counter([a.subst_name for a in atoms if a.subst_name.upper() not in cap_set])
    if central_resname is None:
        if not subst_counts:
            raise ValueError("Could not infer central residue name (MOL2 seems to contain only caps).")
        central_resname = subst_counts.most_common(1)[0][0]

    central_atoms = [a for a in atoms if a.subst_name.upper() not in cap_set and a.subst_name == central_resname]
    if len(central_atoms) < 4:
        raise ValueError(
            f"After omitting caps {cap_resnames}, too few atoms remain ({len(central_atoms)}). "
            f"Inferred central residue: {central_resname}"
        )

    # Determine omit atoms (all cap atoms)
    omit_atom_names: List[str] = [a.name for a in atoms if a.subst_name.upper() in cap_set]

    # Build id<->atom mapping for central residue
    id_to_atom: Dict[int, Mol2Atom] = {a.atom_id: a for a in central_atoms}
    name_to_ids: Dict[str, List[int]] = {}
    for a in central_atoms:
        name_to_ids.setdefault(a.name, []).append(a.atom_id)

    # Build graph restricted to central residue
    graph: Dict[int, Set[int]] = {aid: set() for aid in id_to_atom.keys()}
    for a1, a2 in bonds:
        if a1 in graph and a2 in graph:
            graph[a1].add(a2)
            graph[a2].add(a1)

    # Infer MAIN_CHAIN list if asked
    mainchain_names: List[str] = []
    if main_chain and isinstance(main_chain, list) and len(main_chain) > 0:
        mainchain_names = list(main_chain)
    elif infer_mainchain_from_connectivity:
        # choose first matching ids for head/tail names
        if head_name not in name_to_ids:
            raise ValueError(f"head_name '{head_name}' not found in central residue atoms.")
        if tail_name not in name_to_ids:
            raise ValueError(f"tail_name '{tail_name}' not found in central residue atoms.")

        head_id = name_to_ids[head_name][0]
        tail_id = name_to_ids[tail_name][0]
        path = _shortest_path(graph, head_id, tail_id)
        if not path or len(path) < 2:
            raise ValueError(f"Could not infer a bonded path between head '{head_name}' and tail '{tail_name}'.")

        # MAIN_CHAIN expects names, keep order from head->tail
        mainchain_names = [id_to_atom[i].name for i in path]

    else:
        # safe fallback
        mainchain_names = ["CA"]

    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", encoding="utf-8") as fh:
        fh.write(f"HEAD_NAME {head_name}\n")
        fh.write(f"TAIL_NAME {tail_name}\n")
        for mc in mainchain_names:
            fh.write(f"MAIN_CHAIN {mc}\n")
        for name in omit_atom_names:
            fh.write(f"OMIT_NAME {name}\n")
        fh.write(f"PRE_HEAD_TYPE {pre_head_type}\n")
        fh.write(f"POST_TAIL_TYPE {post_tail_type}\n")
        fh.write(f"CHARGE {float(charge)}\n")

    return str(out)