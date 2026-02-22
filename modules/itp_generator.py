from .utils import get_atomtypes, get_mass, extract_charges, get_bonds_angles_dihedrals, get_impropers, get_pairs, fix_duplicates

def generate_gromacs_itp(mol2_path, frcmod_paths, output_path, resname):
    atoms, atomtypes = get_atomtypes(mol2_path)
    charges = extract_charges(mol2_path)
    masses = [get_mass(a['type']) for a in atoms]
    
    bonds, angles, dihedrals = get_bonds_angles_dihedrals(mol2_path)
    impropers = get_impropers(mol2_path)
    pairs = get_pairs(dihedrals)

    # Fix duplicate issues
    angles = fix_duplicates(angles)
    dihedrals = fix_duplicates(dihedrals)

    # Write ITP file
    with open(output_path, 'w') as f:
        f.write("[ atomtypes ]\n")
        for at in atomtypes:
            f.write(f"{at}\n")

        f.write("\n[ moleculetype ]\n")
        f.write(f"{resname:<15} 3\n")

        f.write("\n[ atoms ]\n")
        for i, atom in enumerate(atoms, 1):
            f.write(f"{i:>5} {atom['type']:<5} 1 {resname:<5} {atom['name']:<5} {i:>5} {charges[i-1]:>10.6f} {masses[i-1]:>10.5f}\n")

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
