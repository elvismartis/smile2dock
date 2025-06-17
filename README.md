# smile2dock: SMILES-to-3D Molecular Format Converter
 ![image](smile2dock.png)
 
A Python tool for converting SMILES strings—single molecules or large libraries—into multiple 3D molecular formats (PDB, PDBQT, MOL2, SDF) for molecular docking and computational chemistry. The tool also calculates key physicochemical and pharmacokinetic properties using RDKit and Open Babel.
Features
- Input: Single SMILES string or a file (one SMILES per line)
- 3D Structure Generation: Generates 3D coordinates and optimizes geometry
- Format Conversion: Exports to PDB, PDBQT, MOL2, and SDF formats
- Property Calculation: Calculates molecular weight, LogP, H-bond donors/acceptors, TPSA, rotatable bonds, and more
- Batch Processing: Handles large compound libraries efficiently
- Extensible: Easily add more descriptors or output formats

## Install the package

```
git clone git@github.com:elvismartis/smile2dock.git
```

## Requirements
- Python 3.7+
- RDKit
- Open Babel (with Python bindings: openbabel and pybel)

## Install dependencies
 ```
python3 -m venv </PATH/>cheminfo
source </PATH/>cheminfo/bin/activate
python3 -m pip install openbabel-wheel pybel
```
## Usage:

### General usage:
```
python smile2dock.py -h 
```

### single molecule conversion

```
python smile2dock.py "CCO" -o ethanol
```
### Single molecule conversion with Tanimoto similarity
```
python smile2dock.py "CCO" -o ethanol --reference "CCN"
```

### Batch Processing for molecular conversion
Given a file `molecules.smi` with one SMILES per line:
```
python smile2dock.py molecules.smi -o output_directory
```

### Batch Processing for molecular conversion with Tanimoto similarity
```
python smile2dock.py molecules.smi -o output_dir --reference "CCN"
```

## Optional Arguments:

### input
```
-n, --num_confs: Number of 3D conformers to generate (default: 10)
-o, --output: Output directory or base name
```

### Output
For each molecule, the following files are generated:
```
pdb: Protein Data Bank format
pdbqt: AutoDock Vina format
mol2: SYBYL MOL2 format
sdf: Structure Data File
```

# Known Issues
## Protonation states: 
- The net charge must be specified for the moment (*workaround: use morphitedl*)

## Incorrect aliphatic ring conformations
- It has been observed that on certain occasions, aliphatic ring conformations are incorrect. (*work in process to fix this*)

