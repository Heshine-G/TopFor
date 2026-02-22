def get_mass(atom_type):
    
    masses = {
        # Hydrogens
        "H": 1.008,
        "HO": 1.008,
        "HC": 1.008,
        "H1": 1.008,
        "H2": 1.008,
        "H3": 1.008,
        "HN": 1.008,

        # Carbons
        "C": 12.011,
        "CA": 12.011,
        "CT": 12.011,
        "CM": 12.011,
        "C2": 12.011,
        "C3": 12.011,

        # Nitrogens
        "N": 14.007,
        "NA": 14.007,
        "N3": 14.007,
        "NB": 14.007,

        # Oxygens
        "O": 15.999,
        "O2": 15.999,
        "OH": 15.999,
        "OS": 15.999,

        # Sulfur 
        "S": 32.06,
        "SH": 32.06,
    }

    return masses.get(atom_type, 12.011)  


def extract_charges(mol2_path):
    charges = []
    with open(mol2_path) as f:
        read_atoms = False
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                read_atoms = True
                continue
            elif "@<TRIPOS>" in line and read_atoms:
                break
            if read_atoms:
                parts = line.split()
                charges.append(float(parts[-1]))
    return charges

def get_atomtypes(mol2_path):
    
    return [], []  

def get_bonds_angles_dihedrals(mol2_path):
    
    return [], [], []

def get_impropers(mol2_path):
    
    return []

def get_pairs(dihedrals):
    
    return []

def fix_duplicates(entries):
    
    seen = set()
    result = []
    for e in entries:
        key = tuple(sorted(e[:3])) if len(e) == 3 else tuple(sorted(e[:4]))
        if key not in seen:
            seen.add(key)
            result.append(e)
    return result
