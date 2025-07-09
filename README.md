## Non-Standard Amino Acid Parameterization Pipeline

### This pipeline automates the process of preparing non-standard amino acids for molecular simulations using AMBER. It:

1. Adds ACE/NME capping groups to the N/C terminal of the amino acid, respectively.

2. Adds hydrogens and partial charges.

3. Generates parameter files via Antechamber and related tools.

### Directory Structure:

        main.py                 
        modules/
        ├── __init__.py                  
        ├── processor.py                 
        ├── capping.py                   
        ├── run_antechamber.py          
        └── remove.py     

### Environment Setup
Follow these steps to set up everything from scratch:

#### 1. Clone the Repository
   
$ git clone https://github.com/Heshine-G/nsaa-paramgen-amber.git

$ cd nsaa-paramgen-amber

#### 2. Create the Conda Environment

Use the following command to create a working environment with all dependencies:

$ conda create -n paramgen python=3.10 numpy pandas rdkit pymol-open-source ambertools -c conda-forge

$ conda activate paramgen

$ python main.py -h


