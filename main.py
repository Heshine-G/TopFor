import argparse
import os
from modules.processor import NonStandardAminoAcidProcessor
from modules.run_antechamber import run_antechamber_for_all

def main():
    parser = argparse.ArgumentParser(description="Full pipeline for processing non-standard amino acids")
    parser.add_argument("input", help="MOL2 file or a .txt file containing list of MOL2s")
    args = parser.parse_args()

    processor = NonStandardAminoAcidProcessor(args.input)
    charged_files = processor.process()

    if charged_files:
        run_antechamber_for_all(charged_files)
    else:
        print("No charged .mol2 files were generated. Skipping parameter generation.")

if __name__ == "__main__":
    main()
