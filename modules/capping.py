import os
import pymol
import sys

pymol.finish_launching(['pymol', '-cq'])
if len(sys.argv) < 2:
    print("Usage: python capping.py <residue_folder>")
    sys.exit(1)

residue_folder = sys.argv[1]
input_file = os.path.join(residue_folder, "residue.mol2")
output_file = os.path.join(residue_folder, "residue_capped.mol2")

if not os.path.exists(input_file):
    print(f"ERROR: {input_file} not found!")
    sys.exit(1)

pymol.cmd.load(input_file, "prot")
pymol.cmd.remove("hydro")
pymol.cmd.select("pk1", "name N")
if pymol.cmd.count_atoms("pk1") > 0:
    pymol.cmd.editor.attach_amino_acid("pk1", "ace")
    pymol.cmd.alter("resn ACE and name C", "name = 'C1'")
    pymol.cmd.alter("resn ACE and not name C1", "name = '1' + name")

pymol.cmd.select("pk1", "name C and not resn ace")
if pymol.cmd.count_atoms("pk1") > 0:
    pymol.cmd.editor.attach_amino_acid("pk1", "nme")
    pymol.cmd.alter("resn NME and name N", "name = 'N1'")

pymol.cmd.h_add("prot")
pymol.cmd.save(output_file, "prot")

if os.path.exists(output_file):
    print(f"Capping successful: {output_file}")
else:
    print(f"ERROR: PyMOL failed to create {output_file}!")