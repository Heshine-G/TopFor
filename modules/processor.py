import os
import sys
import subprocess

class NonStandardAminoAcidProcessor:
    def __init__(self, input_file):
        self.input_file = input_file
        self.output_dir = "non_standard_residues"

    def create_output_directory(self):
        os.makedirs(self.output_dir, exist_ok=True)

    def extract_residue_from_mol2(self, mol2_path):
        resname = os.path.splitext(os.path.basename(mol2_path))[0]
        residue_dir = os.path.join(self.output_dir, resname)
        os.makedirs(residue_dir, exist_ok=True)

        residue_file = os.path.join(residue_dir, "residue.mol2")
        with open(mol2_path, 'r') as fin, open(residue_file, 'w') as fout:
            for line in fin:
                fout.write(line)

        return residue_dir

    def add_terminal_groups_with_pymol(self, residue_folder):
        try:
            script_path = os.path.join(os.path.dirname(__file__), "capping.py")
            result = subprocess.run([sys.executable, script_path, residue_folder], capture_output=True, text=True)
            if result.returncode != 0:
                print(f"[Capping Error] stderr:\n{result.stderr}")
                return None

            capped_mol2 = os.path.join(residue_folder, "residue_capped.mol2")
            if not os.path.exists(capped_mol2) or os.path.getsize(capped_mol2) == 0:
                print("[Capping Error] Output mol2 file is missing or empty.")
                return None

            return capped_mol2
        except Exception as e:
            print(f"Exception during capping: {e}")
            return None

    def run_antechamber(self, input_mol2, residue_name, output_dir):
        output_mol2 = os.path.join(output_dir, f"{residue_name}_charged.mol2")
        cmd = [
            "antechamber", "-fi", "mol2", "-i", input_mol2, "-bk", residue_name,
            "-fo", "mol2", "-o", output_mol2, "-c", "gas", "-at", "amber"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[Antechamber Error] stderr:\n{result.stderr}")
            return None

        if not os.path.exists(output_mol2) or os.path.getsize(output_mol2) == 0:
            print("[Antechamber Error] Output mol2 file is missing or empty.")
            return None

        print(f"[Antechamber] Successfully generated: {output_mol2}")
        return output_mol2

    def process(self):
        self.create_output_directory()
        input_paths = []

        if self.input_file.endswith(".txt"):
            with open(self.input_file, 'r') as f:
                input_paths = [line.strip() for line in f if line.strip()]
        elif self.input_file.endswith(".mol2"):
            input_paths = [self.input_file]
        else:
            print("ERROR: Input must be a .mol2 file or a .txt file with paths to .mol2 files.")
            return []

        final_mol2_files = []
        for mol2_path in input_paths:
            print(f"\nProcessing {mol2_path}")
            residue_folder = self.extract_residue_from_mol2(mol2_path)
            resname = os.path.splitext(os.path.basename(mol2_path))[0].upper()

            capped_file = self.add_terminal_groups_with_pymol(residue_folder)
            if not capped_file:
                continue

            charged_file = self.run_antechamber(capped_file, resname, residue_folder)
            if charged_file:
                final_mol2_files.append(charged_file)

        return final_mol2_files
