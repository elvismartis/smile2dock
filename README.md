# smile2dock: SMILES-to-3D Molecular Format Converter

![image](smile2dock.png)

smile2dock converts SMILES strings (single molecules or libraries) into 3D structures and common molecular file formats (PDB, PDBQT, MOL2, SDF). It also computes several molecular descriptors with RDKit and supports protonation enumeration via Dimorphite-DL.

This repository contains two main improved entry points:

- `smile2dock_v3.py` — Improved core implementation with input validation, logging, resource management and tests.
- `smile2dock_v3_GNNimplicitsolvent.py` — Adds optional GNNImplicitSolvent-based implicit-solvent minimization on top of the v3 base.

## What is new in v3

- Validates SMILES and numeric CLI args
- Centralized constants and configuration
- Timestamped, consistent logging
- Context-managed file operations and explicit RDKit object cleanup
- MMFF geometry optimization with convergence checks
- Protonation enumeration with `dimorphite_dl` (with pH validation)
- Optional GNN-based implicit-solvent minimization support
- Unit tests in `tests/` with a small test runner

## Requirements

- Python 3.8+
- RDKit
- OpenBabel (`openbabel` / `pybel`)
- dimorphite-dl
- numpy

Optional (for GNN features):

- `GNNImplicitSolvent` (installed from GitHub)
- PyTorch, torch-geometric, OpenMM, OpenFF (required by GNNImplicitSolvent)

See `requirements.txt` for recommended versions and the GNN package reference.

## Installation

Clone the repository and create a dedicated environment (recommended: conda/mamba) for platform-sensitive packages such as RDKit and OpenMM. This project now uses PEP 621 (`pyproject.toml`) for packaging.

Quick start (macOS / zsh):

```bash
git clone git@github.com:elvismartis/smile2dock.git
cd smile2dock
# create a conda/mamba environment (example)
mamba create -n smile2dock python=3.10 -y
mamba activate smile2dock
# install platform-sensitive core deps via conda-forge
mamba install -c conda-forge rdkit openbabel -y
# create and activate a lightweight editable install for development
python3 -m pip install --upgrade pip
python3 -m pip install -e '.[dev]'
```

Notes:
- `rdkit`, `openbabel`, `openmm`, and `openff-toolkit` are best installed with conda/mamba.
- The project uses `pyproject.toml` as the single source of packaging metadata (no `setup.py`).
- To install optional GNN extras (heavy ML/MD stack) follow platform-specific steps for PyTorch and torch-geometric first, then:

```bash
# after installing platform-appropriate torch and torch-geometric
python3 -m pip install -e '.[gnn]'
```

If you do not plan to use the GNN-based minimizer, skip the `.[gnn]` install step.

## Example usage

## Reproducible Conda Environments (conda-lock)

This repository includes an `environment.yml` and a workflow that generates a deterministic `conda-lock.yml` for linux/osx/windows. Use `conda-lock` to produce and consume cross-platform lockfiles so CI and contributors can recreate identical conda environments.

Regenerate the lockfile locally (same command used by CI):

```bash
# install conda-lock (pip)
python -m pip install --upgrade pip
python -m pip install conda-lock

# generate conda-lock.yml for common platforms
conda-lock -f environment.yml --platform linux-64 --platform osx-64 --platform win-64 -o conda-lock.yml
```

Install an environment from a lockfile (platform-specific). See the `conda-lock` documentation for full options; an example flow is:

```bash
# example: install the linux-64 lock into a prefix using conda-lock's install helper
conda-lock install --name testenv --file conda-lock.yml --platform linux-64

# or create a prefix using mamba with the explicit package list generated from the lockfile
# (consult conda-lock docs for `conda-lock render` / `conda-lock export` helpers to create exact spec files)
```

Notes:
- CI includes a workflow (`.github/workflows/conda-lock.yml`) that updates `conda-lock.yml` on `main` when you push or via manual dispatch.
- Committing `conda-lock.yml` provides deterministic installs and a stable cache key for CI.


Standard single-molecule conversion (v3):

```bash
python3 smile2dock_v3.py -i "CCO" -o ethanol
```

Protonation enumeration:

```bash
python3 smile2dock_v3.py -i "CCO" --protonate --ph_min 7.0 --ph_max 7.4 -o ethanol_ph
```

GNN implicit-solvent minimization (requires `GNNImplicitSolvent`):

```bash
python3 smile2dock_v3_GNNimplicitsolvent.py -i "CCO" --use_gnn --solvent DMSO -o ethanol_gnn
```

Batch processing:

```bash
python3 smile2dock_v3.py -i molecules.smi -o batch_output --protonate
```

## Command-line flags (summary)

- `-i/--input` : SMILES string or path to a SMILES file (required)
- `-o/--output`: Output directory or base filename
- `-n/--num_confs`: Number of conformers to generate (default: 10)
- `--protonate` : Enable protonation enumeration
- `--ph_min/--ph_max`: pH range for protonation
- `--use_gnn` : Use GNN-based implicit-solvent optimization (GNN script only)
- `--solvent` : Solvent name for GNN implicit-solvent minimization

Run `-h/--help` on either script for full option lists.

## Output

For each molecule the tool will attempt to produce:

- `*.pdb`  — PDB format
- `*.pdbqt` — PDBQT (AutoDock Vina)
- `*.mol2` — MOL2
- `*.sdf`  — SDF

The tool also logs computed properties (MW, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, ring counts, etc.).

## Tests

A basic unit-test suite lives under `tests/`. The GNN tests are skipped automatically if `GNNImplicitSolvent` is not installed.

After installing the `dev` extras you can run the tests:

```bash
python3 -m pytest -q
```

If you prefer to build first, the repository supports PEP 517 builds; the project has produced wheel/sdist artifacts under `dist/` after a local build:

```bash
# build artifacts (in a venv)
python3 -m pip install --upgrade build
python3 -m build
ls -l dist
```

## Installing PyTorch and torch-geometric (guide)

The `.[gnn]` extras depend on PyTorch and `torch-geometric`, which provide platform-specific wheels (CPU / CUDA) and often require matching versions. Follow these steps for common cases.

- CPU-only (recommended for testing without CUDA):

```bash
# install a compatible CPU-only PyTorch wheel
python3 -m pip install "torch>=2.0.0" --index-url https://download.pytorch.org/whl/cpu
# then install torch-geometric following their quickstart (matching the torch version):
python3 -m pip install torch-scatter -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
python3 -m pip install torch-sparse -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
python3 -m pip install torch-cluster -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
python3 -m pip install torch-spline-conv -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
python3 -m pip install torch-geometric
```

- CUDA (GPU) installation (pick CUDA version that matches your system):

1. Install the matching PyTorch CUDA wheel, see https://pytorch.org for instructions. Example for CUDA 11.8:

```bash
python3 -m pip install "torch>=2.0.0+cu118" --index-url https://download.pytorch.org/whl/cu118
```

2. Install torch-geometric wheels that match your torch+CUDA combination. See the PyG installation guide at https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html for exact commands and wheel URLs.

Notes:
- Always consult https://pytorch.org and https://pytorch-geometric.readthedocs.io for the latest supported versions and wheel URLs.
- After installing PyTorch & torch-geometric you can install the gnn extras for the project:

```bash
# after platform-appropriate torch & torch-geometric are installed
python3 -m pip install -e '.[gnn]'
```

## CHANGELOG

See `CHANGELOG.md` for a full list of changes. Notable packaging updates in v3.1.0: `setup.py` was removed in favor of `pyproject.toml` (PEP 621) and the GNN extras Git URL was synced to the fjclark repository for `GNNImplicitSolvent`.

## Citation

If you use Dimorphite-DL, please cite:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small molecules. J Cheminform 11:14. doi: 10.1186/s13321-019-0336-9.

If you use the GNNImplicitSolvent functionality, please cite that project's authors (see their repository for citation details).

## Contact

Open an issue or email: elvis.afmartis@gmail.com

## Warranty Statement
"smile2dock is distributed in the hope that it will be useful, but without any warranty. No author or distributor accepts responsibility to anyone for the consequences of using it or for whether it serves any particular purpose or works at all, unless he says so in writing."