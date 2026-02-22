import os
import sys
import subprocess

class NonStandardAminoAcidProcessor:
    def __init__(self, input_file, charge_model='bcc'):
        self.input_file = input_file
        self.charge_model = charge_model
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

    def convert_pdb_to_mol2(self, pdb_path):
        mol2_path = os.path.splitext(pdb_path)[0] + ".mol2"
        script_path = os.path.join(os.path.dirname(__file__), "pdb_to_mol2.py")
        result = subprocess.run([sys.executable, script_path, pdb_path, mol2_path], capture_output=True, text=True)

        if result.returncode != 0:
            print(f"[Conversion Error] stderr:\n{result.stderr}")
            return None
        if not os.path.exists(mol2_path):
            print("[Conversion Error] Output .mol2 not found after PDB conversion.")
            return None
        return mol2_path

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
        output_mol2 = os.path.join(output_dir, f"{residue_name}.mol2")
        cmd = [
            "antechamber", "-fi", "mol2", "-i", input_mol2, "-bk", residue_name,
            "-fo", "mol2", "-o", output_mol2, "-c", self.charge_model, "-at", "amber"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[Antechamber Error] stderr:\n{result.stderr}")
            return None

        if not os.path.exists(output_mol2) or os.path.getsize(output_mol2) == 0:
            print("[Antechamber Error] Output mol2 file is missing or empty.")
            return None

        return output_mol2

    def _process_mol2_files(self, input_paths):
        final_mol2_files = []
        for mol2_path in input_paths:
            print(f"Processing {mol2_path}")
            residue_folder = self.extract_residue_from_mol2(mol2_path)
            resname = os.path.splitext(os.path.basename(mol2_path))[0].upper()

            capped_file = self.add_terminal_groups_with_pymol(residue_folder)
            if not capped_file:
                continue
            print(f"Capping successful: {capped_file}")

            charged_file = self.run_antechamber(capped_file, resname, residue_folder)
            if charged_file:
                print(f"Charge generation successful: {charged_file}")
                final_mol2_files.append(charged_file)

        return final_mol2_files

    def process_single_residue(self):
        self.create_output_directory()

        if self.input_file.endswith(".mol2"):
            mol2_files = [self.input_file]
        elif self.input_file.endswith(".pdb"):
            mol2_path = self.convert_pdb_to_mol2(self.input_file)
            if mol2_path:
                mol2_files = [mol2_path]
            else:
                return []
        else:
            print("ERROR: Input must be .mol2 or .pdb for -i option.")
            return []

        return self._process_mol2_files(mol2_files)

    def process_peptide(self):
        self.create_output_directory()

        split_script = os.path.join(os.path.dirname(__file__), "split_nonstandard_residues.py")
        input_name = os.path.splitext(os.path.basename(self.input_file))[0]
        residue_dir = os.path.join(self.output_dir, "split_residues", input_name)
        os.makedirs(residue_dir, exist_ok=True)

        result = subprocess.run([sys.executable, split_script, self.input_file, residue_dir], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[Split Error] {result.stderr}")
            return []

        extracted_pdbs = [os.path.join(residue_dir, f) for f in os.listdir(residue_dir) if f.endswith(".pdb")]
        if not extracted_pdbs:
            print("No non-standard residues found in peptide.")
            return []

        mol2_converted = []
        seen_residues = set()

        for pdb_path in extracted_pdbs:
            basename = os.path.basename(pdb_path)
            
            resname = os.path.splitext(basename)[0].split('_')[0].upper()

            if resname in seen_residues:
                print(f"Skipping duplicate residue type: {resname}")
                continue

            seen_residues.add(resname)

            mol2_path = self.convert_pdb_to_mol2(pdb_path)
            if mol2_path:
                mol2_converted.append(mol2_path)

        print(f"\nProcessed unique non-standard residues: {', '.join(sorted(seen_residues))}")


        return self._process_mol2_files(mol2_converted)
