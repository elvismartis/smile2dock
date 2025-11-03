#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SMILES-to-3D Molecular Format Converter v3.0 with GNNImplicitSolvent Support

Enhanced version with improved error handling, resource management,
code organization, and GNN-based implicit solvation minimization.

For Bug reports, send an E-mail to:
elvis.afmartis@gmail.com
"""

import os
import argparse
import logging
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Crippen
from openbabel import pybel
import dimorphite_dl

try:
    from GNNImplicitSolvent import minimize_mol as gnn_minimize_mol
    HAS_GNN = True
except ImportError:
    HAS_GNN = False
    logging.warning("GNNImplicitSolvent not found. GNN-based solvation will be disabled.")

# Constants for default values and configuration
DEFAULT_PH_MIN = 6.4
DEFAULT_PH_MAX = 8.4
DEFAULT_PRECISION = 1.0
DEFAULT_MAX_VARIANTS = 128
DEFAULT_NUM_CONFS = 10
DEFAULT_MMFF_MAXITERS = 1000
SUPPORTED_FORMATS = ["pdb", "mol2", "sdf", "pdbqt"]
SUPPORTED_SOLVENTS = [
    "water", "chloroform", "acetone", "acetonitrile", 
    "ethylacetate", "THF", "DCM", "ethanol", 
    "methanol", "DMSO"
]

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def ensure_output_dir(output_base: str) -> None:
    """Create the output directory if it does not exist."""
    out_dir = os.path.dirname(output_base)
    if out_dir and not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as e:
            logging.error(f"Failed to create output directory {out_dir}: {e}")
            raise

def validate_ph_range(ph_min: float, ph_max: float) -> bool:
    """Validate pH range parameters."""
    if not isinstance(ph_min, (int, float)) or not isinstance(ph_max, (int, float)):
        logging.error("pH values must be numeric")
        return False
    
    if not (0 <= ph_min <= 14 and 0 <= ph_max <= 14):
        logging.error("pH values must be between 0 and 14")
        return False
    
    if ph_min > ph_max:
        logging.error(f"ph_min ({ph_min}) must be less than or equal to ph_max ({ph_max})")
        return False
    
    return True

def is_valid_smiles(smiles: str) -> bool:
    """Check if a string is a valid SMILES using RDKit."""
    if not isinstance(smiles, str):
        return False
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def protonate_smiles(smiles: str, ph_min: float = DEFAULT_PH_MIN, 
                    ph_max: float = DEFAULT_PH_MAX, 
                    precision: float = DEFAULT_PRECISION, 
                    max_variants: int = DEFAULT_MAX_VARIANTS) -> List[str]:
    """Protonate SMILES using Dimorphite-DL for specified pH range."""
    if not validate_ph_range(ph_min, ph_max):
        return [smiles]
        
    try:
        protonated_smiles = dimorphite_dl.protonate_smiles(
            smiles, 
            ph_min=ph_min, 
            ph_max=ph_max, 
            precision=precision,
            max_variants=max_variants
        )
        return protonated_smiles
    except Exception as e:
        logging.error(f"Error protonating {smiles}: {str(e)}")
        return [smiles]  # Return original if protonation fails

def smiles_to_3d(smiles: str, output_base: str = "molecule", 
                num_confs: int = DEFAULT_NUM_CONFS, 
                optimize: bool = True,
                protonate: bool = False,
                ph_min: float = DEFAULT_PH_MIN,
                ph_max: float = DEFAULT_PH_MAX,
                use_gnn: bool = False,
                solvent: str = None) -> Tuple[Optional[Chem.Mol], Optional[Dict], Optional[List[str]]]:
    """
    Convert SMILES to multiple 3D formats and calculate properties.
    
    Args:
        smiles: Input SMILES string
        output_base: Base name for output files
        num_confs: Number of conformers to generate
        optimize: Whether to optimize the conformers
        protonate: Whether to protonate the molecule
        ph_min: Minimum pH for protonation
        ph_max: Maximum pH for protonation
        use_gnn: Whether to use GNN-based solvation for optimization
        solvent: Solvent to use for GNN-based optimization
    
    Returns:
        Tuple of (RDKit mol object, dict of properties, list of protonated variants)
    """
    if not is_valid_smiles(smiles):
        logging.error(f"Invalid SMILES: {smiles}")
        return None, None, None

    if use_gnn and not HAS_GNN:
        logging.error("GNN-based solvation requested but GNNImplicitSolvent not available")
        return None, None, None

    if use_gnn and solvent not in SUPPORTED_SOLVENTS:
        logging.error(f"Unsupported solvent for GNN-based solvation: {solvent}")
        return None, None, None

    protonated_variants = None
    mol = None

    try:
        # Protonate SMILES if requested
        if protonate:
            protonated_variants = protonate_smiles(smiles, ph_min, ph_max)
            logging.info(f"Generated {len(protonated_variants)} protonation states for pH {ph_min}-{ph_max}")
            for i, variant in enumerate(protonated_variants):
                logging.info(f"  Variant {i+1}: {variant}")
            
            if len(protonated_variants) > 1:
                logging.info(f"Using first protonation variant: {protonated_variants[0]}")
            smiles = protonated_variants[0]
        else:
            protonated_variants = [smiles]

        # RDKit: Parse and add hydrogens
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Generate 3D conformers
        try:
            AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_confs,
                randomSeed=42,
                enforceChirality=True
            )
        except Exception as e:
            logging.error(f"Conformer generation failed for {smiles}: {e}")
            return None, None, protonated_variants

        # Optimize conformers
        if optimize:
            if use_gnn and HAS_GNN:
                try:
                    # Use GNNImplicitSolvent for optimization
                    logging.info(f"Using GNN-based solvation with solvent: {solvent}")
                    mol, energies = gnn_minimize_mol(mol, solvent)
                    logging.info(f"GNN optimization completed successfully")
                except Exception as e:
                    logging.error(f"GNN optimization failed: {e}")
                    return None, None, protonated_variants
            else:
                # Use standard MMFF optimization
                for conf_id in range(mol.GetNumConformers()):
                    try:
                        res = AllChem.MMFFOptimizeMolecule(
                            mol,
                            confId=conf_id,
                            mmffVariant='MMFF94s',
                            maxIters=DEFAULT_MMFF_MAXITERS
                        )
                        if res == -1:
                            logging.warning(f"Optimization did not converge for conformer {conf_id}")
                    except Exception as e:
                        logging.warning(f"Optimization failed for conformer {conf_id}: {e}")

        # Convert to OpenBabel molecule for format export
        try:
            sdf_data = Chem.MolToMolBlock(mol)
            ob_mol = pybel.readstring("mol", sdf_data)
        except Exception as e:
            logging.error(f"OpenBabel conversion failed: {e}")
            return None, None, protonated_variants

        # Ensure output directory exists
        ensure_output_dir(output_base)

        # Write different file formats
        for fmt in SUPPORTED_FORMATS:
            output = f"{output_base}.{fmt}"
            try:
                ob_mol.write(fmt, output, overwrite=True)
                logging.info(f"Wrote {fmt.upper()} file: {output}")
            except Exception as e:
                logging.error(f"Failed to write {output}: {e}")

        # Calculate molecular properties
        try:
            properties = {
                "Molecular Weight": Descriptors.ExactMolWt(mol),
                "Crippen_LogP": Crippen.MolLogP(mol),
                "Crippen_MR": Crippen.MolMR(mol),
                "H-Bond Donors": Descriptors.NumHDonors(mol),
                "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
                "TPSA": Descriptors.TPSA(mol),
                "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
                "#Aliphatic Rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
                "#Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
                "#Heteroaromatic Rings": rdMolDescriptors.CalcNumHeterocycles(mol)
            }
        except Exception as e:
            logging.error(f"Property calculation failed: {e}")
            properties = {}

        return mol, properties, protonated_variants

    except Exception as e:
        logging.error(f"Error processing {smiles}: {e}")
        return None, None, protonated_variants
    
    finally:
        if mol is not None:
            del mol

def main():
    parser = argparse.ArgumentParser(
        description='SMILES to 3D converter with protonation, property calculation, and GNN-based solvation'
    )

    # Core arguments
    parser.add_argument('-i', '--input', required=True,
                      help='SMILES string or input file')
    parser.add_argument('-o', '--output', default="output",
                      help='Output directory or base name')
    parser.add_argument('-n', '--num_confs', type=int, default=DEFAULT_NUM_CONFS,
                      help=f'Number of conformers to generate (default: {DEFAULT_NUM_CONFS})')

    # Protonation arguments
    parser.add_argument('--protonate', action='store_true',
                      help='Enable protonation using Dimorphite-DL')
    parser.add_argument('--ph_min', type=float, default=DEFAULT_PH_MIN,
                      help=f'Minimum pH for protonation (default: {DEFAULT_PH_MIN})')
    parser.add_argument('--ph_max', type=float, default=DEFAULT_PH_MAX,
                      help=f'Maximum pH for protonation (default: {DEFAULT_PH_MAX})')

    # GNN solvation arguments
    parser.add_argument('--use_gnn', action='store_true',
                      help='Use GNN-based implicit solvation for optimization')
    parser.add_argument('--solvent', choices=SUPPORTED_SOLVENTS,
                      help='Solvent to use for GNN-based optimization')

    args = parser.parse_args()

    # Validate arguments
    if args.use_gnn and not args.solvent:
        parser.error("--solvent is required when using --use_gnn")

    if args.protonate and not validate_ph_range(args.ph_min, args.ph_max):
        exit(1)

    try:
        smiles = args.input
        if os.path.isfile(args.input):
            with open(args.input) as f:
                smiles = f.read().strip()

        mol, props, variants = smiles_to_3d(
            smiles,
            args.output,
            args.num_confs,
            optimize=True,
            protonate=args.protonate,
            ph_min=args.ph_min,
            ph_max=args.ph_max,
            use_gnn=args.use_gnn,
            solvent=args.solvent
        )

        if props:
            logging.info("\nMolecular Properties:")
            for k, v in props.items():
                logging.info(f"{k}: {v:.2f}")

            if args.protonate and variants:
                logging.info(f"\nProtonation states at pH {args.ph_min}-{args.ph_max}:")
                for i, variant in enumerate(variants, 1):
                    logging.info(f"  {i}. {variant}")

    except Exception as e:
        logging.error(f"Fatal error: {e}")
        exit(1)

if __name__ == "__main__":
    main()