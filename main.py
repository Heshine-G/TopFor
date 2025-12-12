<<<<<<< HEAD
=======
# === FILE: main.py ===
# CHANGES: Added argparse flags for --backbone, --sidechain, --charge and passed them to processor and run_antechamber
>>>>>>> 5e491f921586e2c42c0ba801adfd91fffc7c84ea

import argparse
import os
import glob
from modules.processor import NonStandardAminoAcidProcessor
from modules.run_antechamber import run_antechamber_for_all

def main():
    parser = argparse.ArgumentParser(description="Full pipeline for processing non-standard amino acids")
    group = parser.add_mutually_exclusive_group(required=True)
<<<<<<< HEAD
    group.add_argument('-i', '--input', help="Single non-standard residue (PDB or MOL2)")
    group.add_argument('-p', '--peptide', help="Full peptide PDB file (extracts non-standard residues)")
    group.add_argument('-b', '--batch', help="Batch mode: a .txt file with paths to .mol2 or .pdb files")

    parser.add_argument('--backbone', '-bb', choices=['ff14SB', 'ff19SB', 'ff99SB'], default='ff19SB')
    parser.add_argument('--sidechain', '-sc', choices=['gaff', 'gaff2'], default='gaff2')
    parser.add_argument('--charge', '-c', choices=['gas', 'bcc'], default='bcc')

    args = parser.parse_args()

    charged_files = []

    if args.input:
        processor = NonStandardAminoAcidProcessor(args.input, args.charge)
        charged_files.extend(processor.process_single_residue())

    elif args.peptide:
        processor = NonStandardAminoAcidProcessor(args.peptide, args.charge)
        charged_files.extend(processor.process_peptide())

    elif args.batch:
        with open(args.batch, 'r') as f:
            input_files = [line.strip() for line in f if line.strip()]
        for file_path in input_files:
            processor = NonStandardAminoAcidProcessor(file_path, args.charge)
            charged_files.extend(processor.process_single_residue())
=======
    group.add_argument('-i', '--input', help="Single MOL2 file")
    group.add_argument('-b', '--batch', help="Batch mode: a .txt file listing MOL2 paths or a glob pattern")

    parser.add_argument('--backbone', '-bb', choices=['ff14SB', 'ff19SB', 'ff99SB'], default='ff19SB', help="Backbone forcefield")
    parser.add_argument('--sidechain', '-sc', choices=['gaff', 'gaff2'], default='gaff2', help="Sidechain forcefield")
    parser.add_argument('--charge', '-c', choices=['gas', 'bcc'], default='bcc', help="Charge model")

    args = parser.parse_args()

    if args.input:
        input_files = [args.input]
    else:
        batch_arg = args.batch
        if batch_arg.endswith('.txt'):
            with open(batch_arg, 'r') as f:
                input_files = [line.strip() for line in f if line.strip()]
        else:
            input_files = glob.glob(batch_arg)
        if not input_files:
            print(f"No input files found for batch argument: {batch_arg}")
            return

    charged_files = []
    for file_path in input_files:
        if not os.path.isfile(file_path):
            print(f"ERROR: File not found: {file_path}")
            continue
        processor = NonStandardAminoAcidProcessor(file_path, args.charge)
        charged_files.extend(processor.process())
>>>>>>> 5e491f921586e2c42c0ba801adfd91fffc7c84ea

    if charged_files:
        run_antechamber_for_all(charged_files, args.backbone, args.sidechain, args.charge)
    else:
        print("No charged .mol2 files were generated. Skipping parameter generation.")

if __name__ == "__main__":
    main()
