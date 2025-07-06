# nsaa-paramgen-amber

### Non-Standard Amino Acid Parameterization Pipeline
#### This pipeline automates the process of preparing non-standard amino acids for molecular simulations using AMBER. It:

1. Adds ACE/NME capping groups.

2. Adds hydrogens and partial charges.

3. Generates parameter files via Antechamber and related tools.

main_pipeline.py                 
└── modules/
    ├── __init__.py                  
    ├── processor.py                 
    ├── capping.py                   
    ├── run_antechamber.py          
    └── remove.py     


