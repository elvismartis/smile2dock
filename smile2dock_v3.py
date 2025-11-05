#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SMILES-to-3D Molecular Format Converter v3.0

Enhanced version with improved error handling, resource management,
and code organization.

For Bug reports, send an E-mail to:
elvis.afmartis@gmail.com
"""

import os
import argparse
import logging
from typing import Dict, List, Optional, Tuple, Union
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors, Crippen
from rdkit.Chem.Fingerprints import FingerprintMols
from openbabel import pybel
import dimorphite_dl

# Constants for default values and configuration
DEFAULT_PH_MIN = 6.4
DEFAULT_PH_MAX = 8.4
DEFAULT_PRECISION = 1.0
DEFAULT_MAX_VARIANTS = 128
DEFAULT_NUM_CONFS = 10
DEFAULT_MMFF_MAXITERS = 1000
SUPPORTED_FORMATS = ["pdb", "mol2", "sdf", "pdbqt"]
PROPERTY_DESCRIPTORS = {
    "Molecular Weight": Descriptors.ExactMolWt,
    "Crippen_LogP": Crippen.MolLogP,
    "Crippen_MR": Crippen.MolMR,
    "H-Bond Donors": Descriptors.NumHDonors,
    "H-Bond Acceptors": Descriptors.NumHAcceptors,
    "TPSA": Descriptors.TPSA,
    "Rotatable Bonds": Descriptors.NumRotatableBonds,
    "NUM of Aliphatic Rings": rdMolDescriptors.CalcNumAliphaticRings,
    "NUM of Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings,
    "NUM of Heteroaromatic Rings": rdMolDescriptors.CalcNumHeterocycles
}

# Configure logging with more detailed format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def ensure_output_dir(output_base: str) -> None:
    """Create the output directory if it does not exist.
    
    Args:
        output_base: Base path for output files
    """
    out_dir = os.path.dirname(output_base)
    if out_dir and not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as e:
            logging.error(f"Failed to create output directory {out_dir}: {e}")
            raise

def validate_ph_range(ph_min: float, ph_max: float) -> bool:
    """Validate pH range parameters.
    
    Args:
        ph_min: Minimum pH value
        ph_max: Maximum pH value
        
    Returns:
        True if valid, False otherwise
    """
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
    """Check if a string is a valid SMILES using RDKit.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        True if valid SMILES, False otherwise
    """
    if not isinstance(smiles, str):
        return False
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def protonate_smiles(smiles: str, ph_min: float = DEFAULT_PH_MIN, 
                    ph_max: float = DEFAULT_PH_MAX, 
                    precision: float = DEFAULT_PRECISION, 
                    max_variants: int = DEFAULT_MAX_VARIANTS) -> List[str]:
    """Protonate SMILES using Dimorphite-DL for specified pH range.
    
    Args:
        smiles: Input SMILES string
        ph_min: Minimum pH value
        ph_max: Maximum pH value
        precision: pH precision factor
        max_variants: Maximum number of protonation variants
        
    Returns:
        List of protonated SMILES strings
    """
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

def calculate_molecular_properties(mol: Chem.Mol) -> Dict[str, float]:
    """Calculate molecular properties using RDKit descriptors.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dictionary of property names and their values
    """
    properties = {}
    for name, descriptor in PROPERTY_DESCRIPTORS.items():
        try:
            properties[name] = descriptor(mol)
        except Exception as e:
            logging.warning(f"Failed to calculate {name}: {e}")
            properties[name] = float('nan')
    return properties

def optimize_conformer(mol: Chem.Mol, conf_id: int) -> bool:
    """Optimize a single conformer using MMFF94s.
    
    Args:
        mol: RDKit molecule object
        conf_id: Conformer ID to optimize
        
    Returns:
        True if optimization succeeded, False otherwise
    """
    try:
        res = AllChem.MMFFOptimizeMolecule(
            mol,
            confId=conf_id,
            mmffVariant='MMFF94s',
            maxIters=DEFAULT_MMFF_MAXITERS
        )
        if res == -1:
            logging.warning(f"Optimization did not converge for conformer {conf_id}")
            return False
        return True
    except Exception as e:
        logging.warning(f"Optimization failed for conformer {conf_id}: {e}")
        return False

def smiles_to_3d(smiles: str, output_base: str = "molecule", 
                num_confs: int = DEFAULT_NUM_CONFS, 
                optimize: bool = True, 
                protonate: bool = False, 
                ph_min: float = DEFAULT_PH_MIN, 
                ph_max: float = DEFAULT_PH_MAX) -> Tuple[Optional[Chem.Mol], Optional[Dict], Optional[List[str]]]:
    """Convert SMILES to multiple 3D formats and calculate properties.
    
    Args:
        smiles: Input SMILES string
        output_base: Base name for output files
        num_confs: Number of conformers to generate
        optimize: Whether to optimize the conformers
        protonate: Whether to protonate the molecule
        ph_min: Minimum pH for protonation
        ph_max: Maximum pH for protonation
        
    Returns:
        Tuple of (RDKit mol object, dict of properties, list of protonated variants)
    """
    if not is_valid_smiles(smiles):
        logging.error(f"Invalid SMILES: {smiles}")
        return None, None, None

    mol = None
    protonated_variants = None
    try:
        # Protonate SMILES if requested
        if protonate:
            protonated_variants = protonate_smiles(smiles, ph_min, ph_max)
            logging.info(f"Generated {len(protonated_variants)} protonation states for pH {ph_min}-{ph_max}")
            for i, variant in enumerate(protonated_variants):
                logging.info(f"  Variant {i+1}: {variant}")
            
            # Use the first variant for 3D generation
            if len(protonated_variants) > 1:
                logging.info(f"Using first protonation variant: {protonated_variants[0]}")
            smiles = protonated_variants[0]
        else:
            protonated_variants = [smiles]

        # RDKit: Parse and add hydrogens
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Generate 3D conformers
        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=num_confs,
            randomSeed=42,
            enforceChirality=True,
            returnFailedIds=True
        )
        
        if len(conf_ids) == 0:
            logging.error("Failed to generate any conformers")
            return None, None, protonated_variants

        logging.info(f"Generated {len(conf_ids)} conformers")

        if optimize:
            optimization_results = [optimize_conformer(mol, conf_id) for conf_id in conf_ids]
            successful_opts = sum(1 for result in optimization_results if result)
            logging.info(f"Successfully optimized {successful_opts}/{len(conf_ids)} conformers")

        # Convert to OpenBabel molecule for format export
        sdf_data = Chem.MolToMolBlock(mol)
        ob_mol = pybel.readstring("mol", sdf_data)

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

        # Calculate properties
        properties = calculate_molecular_properties(mol)
        
        return mol, properties, protonated_variants

    except Exception as e:
        logging.error(f"Error processing {smiles}: {e}")
        return None, None, protonated_variants
    
    finally:
        # Cleanup to free memory
        if mol is not None:
            del mol

def calculate_tanimoto_similarity(mol1: Chem.Mol, mol2: Chem.Mol, 
                                fp_type: str = "morgan", 
                                radius: int = 2, 
                                n_bits: int = 2048) -> Optional[float]:
    """Calculate Tanimoto similarity between two molecules.
    
    Args:
        mol1: First molecule
        mol2: Second molecule
        fp_type: Fingerprint type ('morgan' or 'rdkit')
        radius: Morgan fingerprint radius
        n_bits: Number of bits in fingerprint
        
    Returns:
        Tanimoto similarity score or None if calculation fails
    """
    if mol1 is None or mol2 is None:
        logging.warning("One or both molecules are None, cannot calculate similarity")
        return None
        
    try:
        if fp_type == "morgan":
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=n_bits)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=n_bits)
        else:
            fp1 = FingerprintMols.FingerprintMol(mol1)
            fp2 = FingerprintMols.FingerprintMol(mol2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except Exception as e:
        logging.error(f"Error calculating similarity: {e}")
        return None

def batch_process(input_file: str, output_dir: str, 
                 reference_smiles: Optional[str] = None, 
                 fp_type: str = "morgan", 
                 radius: int = 2, 
                 n_bits: int = 2048, 
                 protonate: bool = False, 
                 ph_min: float = DEFAULT_PH_MIN, 
                 ph_max: float = DEFAULT_PH_MAX, 
                 num_confs: int = DEFAULT_NUM_CONFS) -> None:
    """Process a file of SMILES strings with optional protonation and similarity calculation.
    
    Args:
        input_file: Path to input file containing SMILES strings
        output_dir: Directory for output files
        reference_smiles: Optional reference SMILES for similarity calculation
        fp_type: Fingerprint type for similarity calculation
        radius: Morgan fingerprint radius
        n_bits: Number of bits in fingerprint
        protonate: Whether to protonate molecules
        ph_min: Minimum pH for protonation
        ph_max: Maximum pH for protonation
        num_confs: Number of conformers to generate
    """
    if not os.path.exists(input_file):
        logging.error(f"Input file not found: {input_file}")
        return

    ref_mol = None
    if reference_smiles:
        if not is_valid_smiles(reference_smiles):
            logging.error(f"Invalid reference SMILES: {reference_smiles}")
        else:
            ref_mol = Chem.MolFromSmiles(reference_smiles)

    # Create output directory
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        logging.error(f"Failed to create output directory {output_dir}: {e}")
        return

    # Process SMILES with proper resource management
    protonation_file = None
    try:
        if protonate:
            protonation_file_path = os.path.join(output_dir, "protonation_states.txt")
            protonation_file = open(protonation_file_path, 'w')
            protonation_file.write("Original_SMILES\tProtonated_SMILES\tpH_Range\n")

        total_lines = sum(1 for _ in open(input_file))
        with open(input_file) as f:
            for idx, line in enumerate(f, 1):
                smiles = line.strip()
                if not smiles:
                    continue

                logging.info(f"Processing molecule {idx}/{total_lines}")
                base_name = os.path.join(output_dir, f"mol_{idx}")
                
                mol, props, protonated_variants = smiles_to_3d(
                    smiles, base_name, num_confs, protonate=protonate, ph_min=ph_min, ph_max=ph_max
                )
                
                if props:
                    logging.info(f"Properties for molecule {idx}:")
                    for k, v in props.items():
                        logging.info(f"  {k}: {v:.2f}")

                    # Log protonation states
                    if protonate and protonated_variants and protonation_file:
                        for variant in protonated_variants:
                            protonation_file.write(f"{smiles}\t{variant}\t{ph_min}-{ph_max}\n")

                    # Calculate similarity if reference provided
                    if ref_mol and mol:
                        similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                        if similarity is not None:
                            logging.info(f"Tanimoto similarity to reference: {similarity:.4f}")

                if mol:  # Cleanup
                    del mol

    except Exception as e:
        logging.error(f"Error in batch processing: {e}")
    finally:
        if protonation_file:
            protonation_file.close()

def single_process(smiles: str, output_base: str, 
                  num_confs: int = DEFAULT_NUM_CONFS, 
                  reference_smiles: Optional[str] = None, 
                  fp_type: str = "morgan", 
                  radius: int = 2, 
                  n_bits: int = 2048, 
                  protonate: bool = False, 
                  ph_min: float = DEFAULT_PH_MIN, 
                  ph_max: float = DEFAULT_PH_MAX) -> None:
    """Process a single SMILES string with optional protonation and similarity calculation.
    
    Args:
        smiles: Input SMILES string
        output_base: Base name for output files
        num_confs: Number of conformers to generate
        reference_smiles: Optional reference SMILES for similarity calculation
        fp_type: Fingerprint type for similarity calculation
        radius: Morgan fingerprint radius
        n_bits: Number of bits in fingerprint
        protonate: Whether to protonate molecule
        ph_min: Minimum pH for protonation
        ph_max: Maximum pH for protonation
    """
    if not is_valid_smiles(smiles):
        logging.error(f"Invalid input SMILES: {smiles}")
        return

    try:
        mol, props, protonated_variants = smiles_to_3d(
            smiles, output_base, num_confs, protonate=protonate, ph_min=ph_min, ph_max=ph_max
        )

        if props:
            logging.info("Generated files: %s", [f"{output_base}.{fmt}" for fmt in SUPPORTED_FORMATS])
            logging.info("\nMolecular Properties:")
            for k, v in props.items():
                logging.info(f"{k}: {v:.2f}")

            # Display protonation states
            if protonate and protonated_variants:
                logging.info(f"\nProtonation states at pH {ph_min}-{ph_max}:")
                for i, variant in enumerate(protonated_variants, 1):
                    logging.info(f"  {i}. {variant}")

            # Calculate similarity if reference provided
            if reference_smiles:
                if not is_valid_smiles(reference_smiles):
                    logging.error(f"Invalid reference SMILES: {reference_smiles}")
                else:
                    ref_mol = Chem.MolFromSmiles(reference_smiles)
                    similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                    if similarity is not None:
                        logging.info(f"Tanimoto similarity to reference: {similarity:.4f}")

        if mol:  # Cleanup
            del mol
            
    except Exception as e:
        logging.error(f"Error in single processing: {e}")

def validate_args(args: argparse.Namespace) -> bool:
    """Validate command line arguments.
    
    Args:
        args: Parsed command line arguments
        
    Returns:
        True if arguments are valid, False otherwise
    """
    if args.num_confs < 1:
        logging.error("Number of conformers must be positive")
        return False
        
    if args.protonate and not validate_ph_range(args.ph_min, args.ph_max):
        return False
        
    if args.radius < 0:
        logging.error("Fingerprint radius must be non-negative")
        return False
        
    if args.bits < 1:
        logging.error("Number of fingerprint bits must be positive")
        return False
        
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='SMILES to 3D converter with protonation, property calculation, and similarity analysis'
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
    parser.add_argument('--precision', type=float, default=DEFAULT_PRECISION,
                      help=f'pKa precision factor (default: {DEFAULT_PRECISION})')
    parser.add_argument('--max_variants', type=int, default=DEFAULT_MAX_VARIANTS,
                      help=f'Maximum protonation variants (default: {DEFAULT_MAX_VARIANTS})')

    # Similarity arguments
    parser.add_argument('--reference',
                      help='Reference SMILES for Tanimoto similarity')
    parser.add_argument('--fp_type', choices=["morgan", "rdkit"], default="morgan",
                      help='Fingerprint type (default: morgan)')
    parser.add_argument('--radius', type=int, default=2,
                      help='Morgan fingerprint radius (default: 2)')
    parser.add_argument('--bits', type=int, default=2048,
                      help='Fingerprint bit size (default: 2048)')

    args = parser.parse_args()

    if not validate_args(args):
        exit(1)

    try:
        if os.path.isfile(args.input):
            os.makedirs(args.output, exist_ok=True)
            batch_process(
                args.input, args.output, args.reference, args.fp_type, args.radius, args.bits,
                args.protonate, args.ph_min, args.ph_max, args.num_confs
            )
        else:
            single_process(
                args.input, args.output, args.num_confs, args.reference, args.fp_type, args.radius, args.bits,
                args.protonate, args.ph_min, args.ph_max
            )
    except Exception as e:
        logging.error(f"Fatal error: {e}")
        exit(1)