## nsaa-paramgen-amber

#### nsaa-paramgen-amber is a command-line tool for generating AMBER parameters for non-standard amino acids. It automates capping, charge calculation, and parameter file creation for use with AMBER force fields.

### Features

- Single-file mode: Process one MOL2 file at a time.

- Batch mode: Process multiple MOL2 files via glob patterns or a text file list.

- Automated pipeline: Capping, Antechamber, TLeap, PREPGEN, and PARMCHK2 steps are chained automatically.

### Requirements

- Operating System: Linux, macOS, or Windows with WSL/Git Bash

- AMBER Tools: antechamber, tleap, prepgen, parmchk2 installed and in your PATH

- Conda: Miniforge or Anaconda installed

### Installation

#### 1. Clone the repository

<pre><code>git clone https://github.com/yourusername/nsaa-paramgen-amber.git
cd nsaa-paramgen-amber</code></pre>

#### 2. Create and activate the conda environment

<pre><code>conda create -n paramgen python=3.10 numpy pandas rdkit pymol-open-source ambertools -c conda-forge
conda activate paramgen</code></pre>

#### 3. Make the CLI launcher script executable

<pre><code>chmod +x topfor</code></pre>

#### 4. Add the project directory to your PATH

<pre><code>echo 'export PATH="$PATH:$(pwd)"' >> ~/.bashrc
source ~/.bashrc</code></pre>

> Note: If you clone elsewhere, adjust $(pwd) to the absolute path.

### Usage

#### Single-file mode

<pre><code>topfor -i myresidue.mol2</code></pre>

#### Batch mode via glob pattern

<pre><code>topfor -b "*.mol2"</code></pre>

#### Batch mode via text file list

1. Create a text file (e.g., list.txt) with one .mol2 path per line.

2. Run:

<pre><code>topfor -b list.txt</code></pre>

#### Flags

. -i, --input    : Path to a single .mol2 file

. -b, --batch   : Glob pattern or .txt file containing .mol2 paths

### Output

For each input file NAME.mol2, the tool creates a directory NAME/ containing:

    - NAME.ac, NAME.mol2 (charged structures)

    - NAME.lib (TLeap library)

    - NAME.prepin (PREPGEN input)

    - NAME.mc (MC file)

    - NAME.frcmod, NAME_GAFF.frcmod, NAME_FF14SB.frcmod (parameter modification files)

### Examples

# Single file
<pre><code>topfor -i AIB.mol2</code></pre>

# Batch via glob
<pre><code>topfor -b "residues/*.mol2"</code></pre>

# Batch via list
<pre><code>echo "AIB.mol2" > list.txt
echo "FGA.mol2" >> list.txt
topfor -b list.txt</code></pre>

### Troubleshooting

- bash: topfor: command not found: Ensure topfor is executable (chmod +x topfor) and your repo folder is added to PATH.

- CRLF errors: Convert topfor script to LF line endings:

<pre><code> dos2unix topfor</code></pre>
