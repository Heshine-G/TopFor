import argparse
import os
import glob
from modules.processor import NonStandardAminoAcidProcessor
from modules.run_antechamber import run_antechamber_for_all

def main():
    parser = argparse.ArgumentParser(
        description="Full pipeline for processing non-standard amino acids"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-i', '--input',
        help="Single MOL2 file"
    )
    group.add_argument(
        '-b', '--batch',
        help="Batch mode: a .txt file listing MOL2 paths or a glob pattern"
    )
    args = parser.parse_args()

    # Determine input files based on mode
    if args.input:
        input_files = [args.input]
    else:
        # Batch mode
        batch_arg = args.batch
        if batch_arg.endswith('.txt'):
            # Read file paths from text file
            with open(batch_arg, 'r') as f:
                input_files = [line.strip() for line in f if line.strip()]
        else:
            # Treat as glob pattern
            input_files = glob.glob(batch_arg)
        if not input_files:
            print(f"No input files found for batch argument: {batch_arg}")
            return

    # Process each file
    charged_files = []
    for file_path in input_files:
        if not os.path.isfile(file_path):
            print(f"ERROR: File not found: {file_path}")
            continue
        processor = NonStandardAminoAcidProcessor(file_path)
        charged_files.extend(processor.process())

    if charged_files:
        run_antechamber_for_all(charged_files)
    else:
        print("No charged .mol2 files were generated. Skipping parameter generation.")

if __name__ == "__main__":
    main()
