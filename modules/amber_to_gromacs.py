from modules.gromacs_templates import write_top_mdp_files
import os
from rdkit import Chem
from rdkit.Chem import AllChem


def amber_type_to_element(amber_type):
    amber_type = amber_type.upper()

    # Explicit hydrogen types
    if amber_type.startswith('H'):
        return 'H'

    # Explicit carbon types
    if amber_type.startswith('C'):
        return 'C'

    # Explicit nitrogen types
    if amber_type.startswith('N'):
        return 'N'

    # Explicit oxygen types
    if amber_type.startswith('O'):
        return 'O'

    # Sulfur
    if amber_type.startswith('S'):
        return 'S'

    # Safe fallback
    return 'C'



def parse_mol2_atoms_and_write_clean(mol2_file, cleaned_file):
    atoms = []

    with open(mol2_file) as f:
        lines = f.readlines()

    in_atoms = False
    with open(cleaned_file, 'w') as out:
        for line in lines:
            if line.strip().startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                out.write(line)
                continue
            elif line.strip().startswith("@<TRIPOS>BOND"):
                in_atoms = False
                out.write(line)
                continue

            if in_atoms:
                parts = line.split()
                if len(parts) < 9:
                    continue

                atom_id, atom_name = parts[0], parts[1]
                x, y, z = parts[2:5]
                amber_type = parts[5]
                subst_id, subst_name = parts[6], parts[7]
                charge = parts[8]

                element = amber_type_to_element(amber_type)

                out.write(
                    f"{atom_id:>7} {atom_name:<4} {x:>10} {y:>10} {z:>10} "
                    f"{element:<6} {subst_id} {subst_name} {charge}\n"
                )

                atoms.append({
                    "id": int(atom_id),
                    "name": atom_name,
                    "type": amber_type,
                    "charge": float(charge),
                })
            else:
                out.write(line)

    return atoms


def get_pairs(mol):
    pairs = set()
    for atom in mol.GetAtoms():
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        for i in neighbors:
            for j in mol.GetAtomWithIdx(i).GetNeighbors():
                if j.GetIdx() == atom.GetIdx():
                    continue
                for k in j.GetNeighbors():
                    if k.GetIdx() in neighbors or k.GetIdx() == atom.GetIdx():
                        continue
                    pairs.add(tuple(sorted([atom.GetIdx() + 1, k.GetIdx() + 1])))
    return sorted(pairs)


def convert_to_gromacs(mol2_file, frcmod_file, prepin_file, output_dir):
    base = os.path.splitext(os.path.basename(mol2_file))[0]

    cleaned_mol2 = os.path.join(output_dir, f"{base}_cleaned_for_rdkit.mol2")
    atoms_info = parse_mol2_atoms_and_write_clean(mol2_file, cleaned_mol2)

    mol = Chem.MolFromMol2File(cleaned_mol2, removeHs=False)
    if mol is None:
        raise RuntimeError("RDKit failed to load MOL2")

    AllChem.EmbedMolecule(mol)
    conf = mol.GetConformer()

    itp_path = os.path.join(output_dir, f"{base}.itp")
    gro_path = os.path.join(output_dir, f"{base}.gro")

    
    coords = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        coords.append((p.x / 10.0, p.y / 10.0, p.z / 10.0))

    xs, ys, zs = zip(*coords)
    padding = 1.2

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    box_x = (max_x - min_x) + 2 * padding
    box_y = (max_y - min_y) + 2 * padding
    box_z = (max_z - min_z) + 2 * padding

    shift_x = -min_x + padding
    shift_y = -min_y + padding
    shift_z = -min_z + padding

    
    with open(gro_path, "w") as gro:
        gro.write(f"{base}\n")
        gro.write(f"{mol.GetNumAtoms():5d}\n")

        for i, (x, y, z) in enumerate(coords):
            atom_name = atoms_info[i]["name"][:5]

            gro.write(
                f"{1:5d}{base:<5}{atom_name:>5}{i+1:5d}"
                f"{x+shift_x:8.3f}{y+shift_y:8.3f}{z+shift_z:8.3f}\n"
            )

        gro.write(f"{box_x:10.5f}{box_y:10.5f}{box_z:10.5f}\n")

    
    with open(itp_path, "w") as itp:
        itp.write("; Generated from AMBER parameters\n")
        itp.write("; Requires AMBER-compatible force field\n\n")

        itp.write("[ moleculetype ]\n")
        itp.write(f"{base} 3\n\n")

        itp.write("[ atoms ]\n")
        for i, atom in enumerate(mol.GetAtoms()):
            data = atoms_info[i]
            itp.write(
                f"{i+1:<5} {data['type']:<6} 1 {base:<4} "
                f"{data['name']:<4} {i+1:<5} "
                f"{data['charge']:8.4f} {atom.GetMass():.4f}\n"
            )

        itp.write("\n[ bonds ]\n")
        for b in mol.GetBonds():
            itp.write(
                f"{b.GetBeginAtomIdx()+1:<5} {b.GetEndAtomIdx()+1:<5} 1\n"
            )

        

        itp.write("\n[ dihedrals ]\n")
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            for n1 in a1.GetNeighbors():
                if n1.GetIdx() == a2.GetIdx():
                    continue
                for n2 in a2.GetNeighbors():
                    if n2.GetIdx() in (a1.GetIdx(), n1.GetIdx()):
                        continue
                    itp.write(
                        f"{n1.GetIdx()+1:<5} {a1.GetIdx()+1:<5} "
                        f"{a2.GetIdx()+1:<5} {n2.GetIdx()+1:<5} 1\n"
                    )

        itp.write("\n[ pairs ]\n")
        itp.write("; 1-4 interactions\n")
        for i, j in get_pairs(mol):
            itp.write(f"{i:<5} {j:<5}\n")

    write_top_mdp_files(output_dir, base)
    print(f"GROMACS files written: {gro_path}, {itp_path}")
