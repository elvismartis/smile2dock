# smile2dock: SMILES-to-3D Molecular Format Converter

![image](smile2dock.png)

smile2dock converts SMILES strings (single molecules or libraries) into 3D structures and common molecular file formats (PDB, PDBQT, MOL2, SDF). It also computes several molecular descriptors with RDKit and supports protonation enumeration via Dimorphite-DL.

- `smile2dock.py` — Improved core implementation with input validation, logging, resource management and tests.

## What is new in v3

- Validates SMILES and numeric CLI args
- Centralized constants and configuration
- Timestamped, consistent logging
- Context-managed file operations and explicit RDKit object cleanup
- MMFF geometry optimization with convergence checks
- Protonation enumeration with `dimorphite_dl` (with pH validation)

## Requirements

- Python 3.8+
- RDKit
- OpenBabel (`openbabel` / `pybel`)
- dimorphite-dl
- numpy

See `requirements.txt` for recommended versions and the GNN package reference.





## Example usage

Standard single-molecule conversion (v3):

```bash
python3 smile2dock.py -i "CCO" -o ethanol
```

Protonation enumeration:

```bash
python3 smile2dock.py -i "CCO" --protonate --ph_min 7.0 --ph_max 7.4 -o ethanol_ph
```

Batch processing:

```bash
python3 smile2dock.py -i molecules.smi -o batch_output --protonate
```

## Command-line flags (summary)

- `-i/--input` : SMILES string or path to a SMILES file (required)
- `-o/--output`: Output directory or base filename
- `-n/--num_confs`: Number of conformers to generate (default: 10)
- `--protonate` : Enable protonation enumeration
- `--ph_min/--ph_max`: pH range for protonation

Run `-h/--help` on either script for full option lists.

## Output

For each molecule the tool will attempt to produce:

- `*.pdb`  — PDB format
- `*.pdbqt` — PDBQT (AutoDock Vina)
- `*.mol2` — MOL2
- `*.sdf`  — SDF

The tool also logs computed properties (MW, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, ring counts, etc.).

## CHANGELOG

See `CHANGELOG.md` for a full list of changes. Notable packaging updates in v3.1.0: `setup.py` was removed in favor of `pyproject.toml` (PEP 621).

## Citation

If you use Dimorphite-DL, please cite:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small molecules. J Cheminform 11:14. doi: 10.1186/s13321-019-0336-9.

## Contact

Open an issue or email: elvis.afmartis@gmail.com

## Warranty Statement
"smile2dock is distributed in the hope that it will be useful, but without any warranty. No author or distributor accepts responsibility to anyone for the consequences of using it or for whether it serves any particular purpose or works at all, unless he says so in writing."