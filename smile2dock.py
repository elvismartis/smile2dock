#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SMILES-to-3D Molecular Format Converter

For Bug reports, send an E-mail to:
elvis.afmartis@gmail.com
"""

import os
import argparse

# Handle missing dependencies gracefully
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors, Crippen
    from rdkit.Chem.Fingerprints import FingerprintMols
except ImportError as e:
    print("Error: RDKit must be installed. Install with: conda install -c conda-forge rdkit")
    raise

try:
    from openbabel import pybel
except ImportError as e:
    print("Error: Open Babel with Python bindings must be installed. Install with: conda install -c conda-forge openbabel")
    raise

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s"
)

def ensure_output_dir(output_base):
    """Ensure the output directory exists for a given base filename."""
    outdir = os.path.dirname(output_base)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

def smiles_to_3d(smiles, output_base="molecule", num_confs=10, optimize=True):
    """Convert SMILES to multiple 3D formats and calculate properties"""
    try:
        # RDKit: Parse and add hydrogens
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.error(f"Invalid SMILES: {smiles}")
            return None, None
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
            return None, None

        if optimize:
            for conf_id in range(mol.GetNumConformers()):
                try:
                    AllChem.MMFFOptimizeMolecule(
                        mol,
                        confId=conf_id,
                        mmffVariant='MMFF94s',
                        maxIters=1000
                    )
                except Exception as e:
                    logging.warning(f"Optimization failed for conformer {conf_id} of {smiles}: {e}")

        # Convert to OpenBabel molecule for format export
        try:
            sdf_data = Chem.MolToMolBlock(mol)
            ob_mol = pybel.readstring("mol", sdf_data)
        except Exception as e:
            logging.error(f"OpenBabel conversion failed: {e}")
            return None, None

        # Ensure output directory exists
        ensure_output_dir(output_base)

        # Write different file formats
        for fmt in ["pdb", "mol2", "sdf", "pdbqt"]:
            output = f"{output_base}.{fmt}"
            try:
                ob_mol.write(fmt, output, overwrite=True)
            except Exception as e:
                logging.error(f"Failed to write {output}: {e}")

        # Calculate properties (corrected descriptors)
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

        return mol, properties

    except Exception as e:
        logging.error(f"Error processing {smiles}: {e}")
        return None, None

def calculate_tanimoto_similarity(mol1, mol2, fp_type="morgan", radius=2, n_bits=2048):
    """Calculate Tanimoto similarity between two molecules"""
    if mol1 is None or mol2 is None:
        logging.warning("One or both molecules are None, cannot calculate similarity.")
        return None
    if fp_type == "morgan":
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=n_bits)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=n_bits)
    else:
        fp1 = FingerprintMols.FingerprintMol(mol1)
        fp2 = FingerprintMols.FingerprintMol(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def batch_process(input_file, output_dir, reference_smiles=None, fp_type="morgan", radius=2, n_bits=2048):
    """Process a file of SMILES strings, optionally computing similarity to a reference"""
    ref_mol = None
    if reference_smiles:
        ref_mol = Chem.MolFromSmiles(reference_smiles)
        if ref_mol is None:
            logging.error(f"Invalid reference SMILES: {reference_smiles}")
            ref_mol = None

    with open(input_file) as f:
        for idx, line in enumerate(f):
            smiles = line.strip()
            if not smiles:
                continue
            base_name = os.path.join(output_dir, f"mol_{idx+1}")
            mol, props = smiles_to_3d(smiles, base_name)
            if props:
                logging.info(f"Processed {smiles}")
                logging.info("Properties:")
                for k, v in props.items():
                    if isinstance(v, float):
                        logging.info(f"  {k}: {v:.2f}")
                    else:
                        logging.info(f"  {k}: {v}")
                if ref_mol:
                    similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                    if similarity is not None:
                        logging.info(f"Tanimoto similarity to reference: {similarity:.4f}")
                    else:
                        logging.warning("Similarity could not be calculated.")

def single_process(smiles, output_base, num_confs, reference_smiles=None, fp_type="morgan", radius=2, n_bits=2048):
    """Process a single SMILES string, optionally computing similarity to a reference"""
    mol, props = smiles_to_3d(smiles, output_base, num_confs)
    if props:
        print("Generated files:", [f"{output_base}.{fmt}" for fmt in ["pdb", "mol2", "sdf", "pdbqt"]])
        print("\nMolecular Properties:")
        for k, v in props.items():
            if isinstance(v, float):
                print(f"{k}: {v:.2f}")
            else:
                print(f"{k}: {v}")
        if reference_smiles:
            ref_mol = Chem.MolFromSmiles(reference_smiles)
            if ref_mol:
                similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                if similarity is not None:
                    print(f"Tanimoto similarity to reference: {similarity:.4f}")
                else:
                    print("Similarity could not be calculated.")
            else:
                print(f"Invalid reference SMILES: {reference_smiles}")

def is_valid_smiles(smiles):
    """Check if a string is a valid SMILES using RDKit."""
    return Chem.MolFromSmiles(smiles) is not None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SMILES to 3D converter with property and similarity calculation')
    parser.add_argument('-i', '--input', required=True, help='SMILES string or input file')
    parser.add_argument('-o', '--output', default="output", help='Output directory or base name')
    parser.add_argument('-n', '--num_confs', type=int, default=10, help='Number of conformers to generate')
    parser.add_argument('--reference', help='Reference SMILES for Tanimoto similarity')
    parser.add_argument('--fp_type', choices=["morgan", "rdkit"], default="rdkit", help='Fingerprint type for similarity')
    parser.add_argument('--radius', type=int, default=2, help='Morgan fingerprint radius (default: 2)')
    parser.add_argument('--bits', type=int, default=2048, help='Fingerprint bit size (default: 2048)')

    args = parser.parse_args()

    if os.path.isfile(args.input):
        os.makedirs(args.output, exist_ok=True)
        batch_process(args.input, args.output, args.reference, args.fp_type, args.radius, args.bits)
    else:
        if not is_valid_smiles(args.input):
            logging.error("Input is not a valid file or SMILES string. Please provide a valid SMILES or input file.")
            exit(1)
        ensure_output_dir(args.output)
        single_process(args.input, args.output, args.num_confs, args.reference, args.fp_type, args.radius, args.bits)
