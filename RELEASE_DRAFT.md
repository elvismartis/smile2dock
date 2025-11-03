# Release draft — v3.1.0

Release: v3.1.0
Date: 2025-11-03

Summary
-------
- PEP 621 packaging: project now uses `pyproject.toml` as the single source of packaging metadata. `setup.py` removed.
- Added `smile2dock_v3.py` (core v3) and `smile2dock_v3_GNNimplicitsolvent.py` (optional GNN-based implicit-solvent minimizer).
- Added `requirements.txt` and optional extras `[gnn]` and `[dev]` in `pyproject.toml`.
- Unit tests and a lightweight test harness were added under `tests/`.
- Local PEP 517 builds produced wheel and sdist artifacts under `dist/`.

Artifacts
---------
- `dist/smile2dock-3.1.0-py3-none-any.whl`
- `dist/smile2dock-3.1.0.tar.gz`

Notes for release page
----------------------
- Installation instructions recommend using conda/mamba for RDKit/OpenMM/OpenFF.
- The `.[gnn]` extras depend on PyTorch and torch-geometric; include platform-specific install steps in the release notes (CPU vs CUDA wheels).
- Tests require RDKit and other chemistry toolkits; recommend running tests in a conda environment with RDKit installed.

Changelog excerpt
-----------------
See `CHANGELOG.md` for full details. Key items:

- v3.0.0: core refactor, validation, logging, MMFF improvements, protonation handling.
- v3.1.0: GNN integration, packaging updates, tests, and pyproject-based builds.

Publishing checklist
-------------------
1. Ensure CI builds and test matrix (Linux/macOS/Windows) install RDKit via conda and run unit tests.
2. Verify wheel builds on target platforms and include platform notes for optional GNN extras.
3. Tag the repository and create the GitHub release with the contents above and attach generated artifacts.
