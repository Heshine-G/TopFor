def process_mol2_file(file_path, output_file):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_start_index = None
    bond_index = None
    for i, line in enumerate(lines):
        if "@<TRIPOS>ATOM" in line:
            atom_start_index = i + 1
        elif "@<TRIPOS>BOND" in line:
            bond_index = i
            break

    if atom_start_index is None or bond_index is None:
        raise ValueError("Missing ATOM or BOND section in MOL2 file.")

    atom_lines = lines[atom_start_index:bond_index]
    omit_lines = []
    for atom_line in atom_lines[:6] + atom_lines[-6:]:
        parts = atom_line.split()
        if len(parts) >= 2:
            omit_lines.append(f"OMIT_NAME {parts[1]}")

    with open(output_file, 'w') as out_file:
        out_file.write("HEAD_NAME N\nTAIL_NAME C\nMAIN_CHAIN CA\n")
        for omit in omit_lines:
            out_file.write(omit + '\n')
        out_file.write("PRE_HEAD_TYPE C\nPOST_TAIL_TYPE N\nCHARGE 0.0\n")