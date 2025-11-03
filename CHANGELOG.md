# Changelog

## [Version 3.0.0] - 2025-11-03

### Added
- Type hints throughout the codebase for better IDE support and code maintainability
- Comprehensive docstrings for all functions with parameter descriptions
- Module-level constants for configuration values
- New validation function `validate_ph_range()` for pH parameters
- New function `calculate_molecular_properties()` to centralize property calculations
- New function `optimize_conformer()` for individual conformer optimization
- New function `validate_args()` for command-line argument validation
- Progress tracking for batch processing operations
- Detailed logging format with timestamps
- Conformer generation success tracking
- Optimization convergence checking
- Resource cleanup mechanisms in finally blocks
- Proper context managers for file operations

### Changed
- Improved error handling throughout the codebase
- Replaced print statements with proper logging
- Restructured main functions for better maintainability
- Enhanced SMILES validation process
- Improved memory management with explicit cleanup
- Better organization of molecular property calculations
- More robust file handling with proper error checks
- Enhanced batch processing with progress reporting
- Improved conformer generation and optimization workflow
- Better handling of protonation states
- More consistent error reporting patterns

### Fixed
- Memory leaks from RDKit mol objects
- Inconsistent error handling in file operations
- Missing validation for pH ranges
- Missing validation for numerical parameters
- Lack of cleanup in batch processing
- Incomplete error handling in SMILES processing
- Missing checks for optimization convergence
- Inconsistent output formatting

### Technical Improvements
1. Code Structure
   - Organized imports by functionality
   - Grouped related functions together
   - Added type annotations for function parameters and returns
   - Centralized configuration in module-level constants

2. Error Handling
   - Added try-except blocks for file operations
   - Added validation for all input parameters
   - Improved error messages with more context
   - Added proper cleanup in error cases

3. Performance
   - Added tracking of conformer generation success
   - Improved memory management
   - Added progress reporting for long operations
   - Enhanced batch processing efficiency

4. Documentation
   - Added detailed function docstrings
   - Improved inline comments
   - Added parameter descriptions
   - Added return value documentation

### Constants Added
- `DEFAULT_PH_MIN = 6.4`
- `DEFAULT_PH_MAX = 8.4`
- `DEFAULT_PRECISION = 1.0`
- `DEFAULT_MAX_VARIANTS = 128`
- `DEFAULT_NUM_CONFS = 10`
- `DEFAULT_MMFF_MAXITERS = 1000`
- `SUPPORTED_FORMATS = ["pdb", "mol2", "sdf", "pdbqt"]`
- Added `PROPERTY_DESCRIPTORS` dictionary for molecular properties

### New Functions
1. `validate_ph_range(ph_min: float, ph_max: float) -> bool`
   - Validates pH range parameters
   - Ensures values are between 0 and 14
   - Checks that minimum pH is less than maximum pH

2. `calculate_molecular_properties(mol: Chem.Mol) -> Dict[str, float]`
   - Centralizes property calculation logic
   - Handles calculation errors gracefully
   - Returns dictionary of properties

3. `optimize_conformer(mol: Chem.Mol, conf_id: int) -> bool`
   - Handles single conformer optimization
   - Checks for optimization convergence
   - Returns success status

4. `validate_args(args: argparse.Namespace) -> bool`
   - Validates command line arguments
   - Checks numerical parameter ranges
   - Returns validation status

### Modified Functions
1. `smiles_to_3d()`
   - Added type hints
   - Improved error handling
   - Added memory management
   - Enhanced conformer generation

2. `batch_process()`
   - Added progress tracking
   - Improved file handling
   - Added proper resource cleanup
   - Enhanced error reporting

3. `single_process()`
   - Added input validation
   - Improved error handling
   - Enhanced output formatting

4. `protonate_smiles()`
   - Added parameter validation
   - Improved error handling
   - Enhanced return value handling

### Dependencies
- No changes to external dependencies
- Same requirements as version 2:
  - RDKit
  - OpenBabel
  - dimorphite_dl

### Notes
- All changes maintain backward compatibility
- No changes to command-line interface arguments
- Output formats remain the same
- Existing scripts using the tool should continue to work as before

## [Version 3.1.0] - 2025-11-03

### Added
- `smile2dock_v3_GNNimplicitsolvent.py` — new entry point that integrates the GNNImplicitSolvent package for implicit-solvent minimization of RDKit conformers
- `requirements.txt` — consolidated dependency manifest including optional GNN installation from GitHub
- Unit tests and test runner under `tests/` to validate core behaviors and provide CI-friendly checks
- `README.md` updated with v3 and GNN usage, installation notes and examples

### Changed
- Documented optional GNN features and heavy ML/MD dependency requirements in the documentation and README
- Improved packaging and developer instructions (use of conda/mamba recommended for RDKit/OpenMM/OpenFF)

### Fixed
- Clarified behavior when `GNNImplicitSolvent` is not installed: the GNN-enabled script gracefully reports unavailability and falls back (or aborts when requested)

### Notes & Caveats
- GNNImplicitSolvent integration is optional and depends on heavy third-party toolchains (PyTorch, torch-geometric, OpenMM, OpenFF). Installing these may require a GPU-enabled environment and careful dependency management (conda recommended).
- Unit tests are lightweight and skip GNN-specific tests if `GNNImplicitSolvent` is not importable.

### Files added/updated
- `smile2dock_v3_GNNimplicitsolvent.py` (new): v3 base + optional GNN-based implicit solvation minimization
- `smile2dock_v3.py` (existing v3 core): minor refinements and validation
- `requirements.txt` (new): pinned suggestions and Git install for `GNNImplicitSolvent`
- `tests/` (new): unit tests and small test runner `tests/run_tests.py`
- `README.md` (updated): documents new features, usage, and testing

### Packaging & build

- Removed `setup.py` in favor of `pyproject.toml` (PEP 621) as the single source of packaging metadata.
- Added a `pyproject.toml` with build-system, project metadata, dependencies, and optional extras (`gnn`, `dev`).
- Synced `GNNImplicitSolvent` Git URL in `pyproject.toml` to `https://github.com/fjclark/GNNImplicitSolvent.git` to match prior configuration.
- Local PEP 517 builds were performed and produced artifacts under `dist/` (wheel and sdist).

### How to test
1. Create a dedicated environment (recommended: mamba/conda) and install core deps
2. (Optional) Install GNN stack if you plan to use the `--use_gnn` mode
3. Run the included tests:

```bash
pip install -r requirements.txt
python3 -m tests.run_tests
```

If GNN is not installed the suite will still run but skip GNN-specific checks.