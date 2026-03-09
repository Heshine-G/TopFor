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

            atoms.append(
                Mol2Atom(
                    atom_id=int(parts[0]),
                    name=parts[1],
                    atom_type=parts[5],
                    subst_name=parts[7],
                )
            )

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

            bonds.append((int(parts[1]), int(parts[2])))

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
                path = [goal]
                cur = goal

                while cur != start:
                    cur = prev[cur]
                    path.append(cur)

                path.reverse()
                return path

            q.append(n)

    return None


def _looks_like_c_or_n_atom(atom: Mol2Atom) -> bool:
    """
    Use BOTH atom name and atom type, because different MOL2 writers vary.
    Accept common carbon/nitrogen names/types such as:
      C, CA, CB, CG, C1, C2, N, N1, ND, NE, NT, etc.
    """
    name = atom.name.upper()
    atype = atom.atom_type.upper()

    return (
        name.startswith("C")
        or name.startswith("N")
        or atype.startswith("C")
        or atype.startswith("N")
    )


def _dedupe_preserve_order(items: List[str]) -> List[str]:
    out: List[str] = []
    seen: Set[str] = set()

    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)

    return out


def _select_anchor_atom(atoms: List[Mol2Atom]) -> Optional[str]:
    """
    Last-resort fallback if no C/N atoms are found.

    Priority:
    1. atom name starts with C
    2. atom type starts with C
    3. any heavy atom except oxygen
    4. any heavy atom
    """
    for a in atoms:
        if a.name.upper().startswith("C"):
            return a.name

    for a in atoms:
        if a.atom_type.upper().startswith("C"):
            return a.name

    for a in atoms:
        if not a.atom_type.upper().startswith("H") and not a.name.upper().startswith("O"):
            return a.name

    for a in atoms:
        if not a.atom_type.upper().startswith("H"):
            return a.name

    return None


def _collect_cn_mainchain_atoms(
    atoms: List[Mol2Atom],
    *,
    exclude_names: Set[str],
) -> List[str]:
    """
    For residues where head and/or tail are absent:
    MAIN_CHAIN should contain all C and N atoms other than head/tail.
    """
    names: List[str] = []

    for atom in atoms:
        if atom.name in exclude_names:
            continue
        if _looks_like_c_or_n_atom(atom):
            names.append(atom.name)

    return _dedupe_preserve_order(names)


def process_mol2_file(
    file_path: str,
    output_file: str,
    *,
    head_name: str | None = "N",
    tail_name: str | None = "C",
    main_chain: Optional[List[str]] = None,
    charge: float = 0.0,
    central_resname: Optional[str] = None,
    cap_resnames: Tuple[str, ...] = ("ACE", "NME"),
    pre_head_type: str = "C",
    post_tail_type: str = "N",
    infer_mainchain_from_connectivity: bool = True,
) -> str:
    atoms = _parse_mol2_atoms(file_path)
    bonds = _parse_mol2_bonds(file_path)

    cap_set = {c.upper() for c in cap_resnames}

    subst_counts = Counter(a.subst_name for a in atoms if a.subst_name.upper() not in cap_set)
    if central_resname is None:
        if not subst_counts:
            raise ValueError("Could not infer central residue name from MOL2.")
        central_resname = subst_counts.most_common(1)[0][0]

    central_atoms = [
        a for a in atoms
        if a.subst_name.upper() not in cap_set and a.subst_name == central_resname
    ]

    if not central_atoms:
        raise ValueError(f"No central residue atoms found for residue '{central_resname}'.")

    omit_atom_names = [a.name for a in atoms if a.subst_name.upper() in cap_set]

    id_to_atom = {a.atom_id: a for a in central_atoms}

    name_to_ids: Dict[str, List[int]] = {}
    for a in central_atoms:
        name_to_ids.setdefault(a.name, []).append(a.atom_id)

    graph: Dict[int, Set[int]] = {a.atom_id: set() for a in central_atoms}
    for a1, a2 in bonds:
        if a1 in graph and a2 in graph:
            graph[a1].add(a2)
            graph[a2].add(a1)

    has_head = bool(head_name) and head_name in name_to_ids
    has_tail = bool(tail_name) and tail_name in name_to_ids

    mainchain_names: List[str] = []

    # 1) Explicit user-defined mainchain always wins
    if main_chain:
        mainchain_names = [x for x in main_chain if x in name_to_ids]

    # 2) Standard polymeric case: both head and tail exist
    elif infer_mainchain_from_connectivity and has_head and has_tail:
        head_id = name_to_ids[head_name][0]
        tail_id = name_to_ids[tail_name][0]

        path = _shortest_path(graph, head_id, tail_id)

        if path and len(path) >= 3:
            path_names = [id_to_atom[i].name for i in path]

            # HEAD_NAME and TAIL_NAME must NOT appear in MAIN_CHAIN
            mainchain_names = [
                n for n in path_names
                if n != head_name and n != tail_name
            ]

    # 3) If head and/or tail are absent:
    #    MAIN_CHAIN = all C/N atoms except head/tail
    if not mainchain_names and ((not has_head) or (not has_tail)):
        exclude_names: Set[str] = set()
        if has_head and head_name:
            exclude_names.add(head_name)
        if has_tail and tail_name:
            exclude_names.add(tail_name)

        mainchain_names = _collect_cn_mainchain_atoms(
            central_atoms,
            exclude_names=exclude_names,
        )

    # 4) Final fallback
    if not mainchain_names:
        anchor = _select_anchor_atom(central_atoms)
        if anchor:
            mainchain_names = [anchor]

    mainchain_names = _dedupe_preserve_order(mainchain_names)

    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", encoding="utf-8") as fh:
        if has_head:
            fh.write(f"HEAD_NAME {head_name}\n")

        if has_tail:
            fh.write(f"TAIL_NAME {tail_name}\n")

        for mc in mainchain_names:
            fh.write(f"MAIN_CHAIN {mc}\n")

        for name in omit_atom_names:
            fh.write(f"OMIT_NAME {name}\n")

        if has_head:
            fh.write(f"PRE_HEAD_TYPE {pre_head_type}\n")

        if has_tail:
            fh.write(f"POST_TAIL_TYPE {post_tail_type}\n")

        fh.write(f"CHARGE {float(charge)}\n")

    return str(out)