# === FILE: pdb_to_mol2.py ===
import pymol
import os
import sys

pymol.finish_launching(['pymol', '-cq'])

if len(sys.argv) < 3:
    print("Usage: python pdb_to_mol2.py <input_pdb_path> <output_mol2_path>")
    sys.exit(1)

input_pdb = sys.argv[1]
output_mol2 = sys.argv[2]

if not os.path.exists(input_pdb):
    print(f"ERROR: {input_pdb} not found!")
    sys.exit(1)

pymol.cmd.load(input_pdb, "prot")
pymol.cmd.h_add("prot")
pymol.cmd.save(output_mol2, "prot")
print(f"Conversion successful: {output_mol2}")
