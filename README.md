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
   
<pre><code>$ git clone https://github.com/Heshine-G/nsaa-paramgen-amber.git
$ cd nsaa-paramgen-amber </code></pre>

#### 2. Create the Conda Environment

Use the following command to create a working environment with all dependencies:

<pre><code>$ conda create -n paramgen python=3.10 numpy pandas rdkit pymol-open-source ambertools -c conda-forge
$ conda activate paramgen
$ python main.py -h</code></pre>

### Input Requirements
You need either:

- A single .mol2 file of a non-standard amino acid (cleaned: no hydrogens, no terminal OXT).

- Or a .txt file containing a list of .mol2 file paths (one per line).

### Run the Pipeline
Assuming you have your .mol2 files in the current working directory:

#### Single file input:
<pre><code>$ python main.py DAB.mol2</code></pre>

#### Batch input with list:
<pre><code>$ python main.py mol2_list.txt</code></pre>
