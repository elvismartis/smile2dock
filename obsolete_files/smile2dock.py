import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors, Crippen
from rdkit.Chem.Fingerprints import FingerprintMols
from openbabel import pybel

import dimorphite_dl

def protonate_smiles(smiles, ph_min=6.4, ph_max=8.4, precision=1.0, max_variants=128):
    """Protonate SMILES using Dimorphite-DL for specified pH range"""
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
        print(f"Error protonating {smiles}: {str(e)}")
        return [smiles]  # Return original if protonation fails

def smiles_to_3d(smiles, output_base="molecule", num_confs=10, optimize=True, protonate=False, ph_min=6.4, ph_max=8.4):
    """Convert SMILES to multiple 3D formats and calculate properties"""
    try:
        # Protonate SMILES if requested
        if protonate:
            protonated_variants = protonate_smiles(smiles, ph_min, ph_max)
            print(f"Generated {len(protonated_variants)} protonation states for pH {ph_min}-{ph_max}")
            for i, variant in enumerate(protonated_variants):
                print(f"  Variant {i+1}: {variant}")
            
            # Use the first variant for 3D generation (or let user choose)
            if len(protonated_variants) > 1:
                print(f"Using first protonation variant: {protonated_variants[0]}")
            smiles = protonated_variants[0]
        
        # RDKit: Parse and add hydrogens
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            return None, None, None
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformers
        AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, randomSeed=42)
        if optimize:
            for conf_id in range(mol.GetNumConformers()):
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)

        # Convert to OpenBabel molecule for format export
        sdf_data = Chem.MolToMolBlock(mol)
        ob_mol = pybel.readstring("mol", sdf_data)
        
        # Write different file formats
        for fmt in ["pdb", "mol2", "sdf", "pdbqt"]:
            output = f"{output_base}.{fmt}"
            ob_mol.write(fmt, output, overwrite=True)
        
        # Calculate properties
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
            "#Heteroaromatic Rings": rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
        }

        # Return protonated variants if generated
        protonated_variants = protonate_smiles(smiles, ph_min, ph_max) if protonate else None
        return mol, properties, protonated_variants

    except Exception as e:
        print(f"Error processing {smiles}: {str(e)}")
        return None, None, None

def calculate_tanimoto_similarity(mol1, mol2, fp_type="morgan", radius=2, n_bits=2048):
    """Calculate Tanimoto similarity between two molecules"""
    if mol1 is None or mol2 is None:
        return None
    if fp_type == "morgan":
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=n_bits)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=n_bits)
    else:
        fp1 = FingerprintMols.FingerprintMol(mol1)
        fp2 = FingerprintMols.FingerprintMol(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def batch_process(input_file, output_dir, reference_smiles=None, fp_type="morgan", radius=2, n_bits=2048, 
                 protonate=False, ph_min=6.4, ph_max=8.4, num_confs=10):
    """Process a file of SMILES strings with optional protonation and similarity calculation"""
    if reference_smiles:
        ref_mol = Chem.MolFromSmiles(reference_smiles)
        if ref_mol is None:
            print(f"Invalid reference SMILES: {reference_smiles}")
            ref_mol = None
    else:
        ref_mol = None

    # Create protonation output file if requested
    protonation_file = None
    if protonate:
        protonation_file = open(os.path.join(output_dir, "protonation_states.txt"), 'w')
        protonation_file.write(f"Original_SMILES\tProtonated_SMILES\tpH_Range\n")

    with open(input_file) as f:
        for idx, line in enumerate(f):
            smiles = line.strip()
            if smiles:
                base_name = os.path.join(output_dir, f"mol_{idx+1}")
                mol, props, protonated_variants = smiles_to_3d(
                    smiles, base_name, num_confs, protonate=protonate, ph_min=ph_min, ph_max=ph_max
                )
                
                if props:
                    print(f"Processed {smiles}")
                    print("Properties:", props)
                    
                    # Log protonation states
                    if protonate and protonated_variants and protonation_file:
                        for variant in protonated_variants:
                            protonation_file.write(f"{smiles}\t{variant}\t{ph_min}-{ph_max}\n")
                    
                    # Calculate similarity if reference provided
                    if ref_mol:
                        similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                        if similarity:
                            print(f"Tanimoto similarity to reference: {similarity:.4f}")

    if protonation_file:
        protonation_file.close()
        print(f"Protonation states saved to: {os.path.join(output_dir, 'protonation_states.txt')}")

def single_process(smiles, output_base, num_confs, reference_smiles=None, fp_type="morgan", radius=2, n_bits=2048,
                  protonate=False, ph_min=6.4, ph_max=8.4):
    """Process a single SMILES string with optional protonation and similarity calculation"""
    mol, props, protonated_variants = smiles_to_3d(
        smiles, output_base, num_confs, protonate=protonate, ph_min=ph_min, ph_max=ph_max
    )
    
    if props:
        print("Generated files:", [f"{output_base}.{fmt}" for fmt in ["pdb", "mol2", "sdf", "pdbqt"]])
        print("\nMolecular Properties:")
        for k, v in props.items():
            print(f"{k}: {v:.2f}")
        
        # Display protonation states
        if protonate and protonated_variants:
            print(f"\nProtonation states at pH {ph_min}-{ph_max}:")
            for i, variant in enumerate(protonated_variants, 1):
                print(f"  {i}. {variant}")
        
        # Calculate similarity if reference provided
        if reference_smiles:
            ref_mol = Chem.MolFromSmiles(reference_smiles)
            if ref_mol:
                similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                if similarity:
                    print(f"Tanimoto similarity to reference: {similarity:.4f}")
            else:
                print(f"Invalid reference SMILES: {reference_smiles}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='SMILES to 3D converter with protonation, property calculation, and similarity analysis'
    )
    
    # Core arguments
    parser.add_argument('-i', '--input', help='SMILES string or input file')
    parser.add_argument('-o', '--output', default="output", help='Output directory or base name')
    parser.add_argument('-n', '--num_confs', type=int, default=10, help='Number of conformers to generate')
    
    # Protonation arguments
    parser.add_argument('--protonate', action='store_true', help='Enable protonation using Dimorphite-DL')
    parser.add_argument('--ph_min', type=float, default=6.4, help='Minimum pH for protonation (default: 6.4)')
    parser.add_argument('--ph_max', type=float, default=8.4, help='Maximum pH for protonation (default: 8.4)')
    parser.add_argument('--precision', type=float, default=1.0, help='pKa precision factor (default: 1.0)')
    parser.add_argument('--max_variants', type=int, default=128, help='Maximum protonation variants (default: 128)')
    
    # Similarity arguments
    parser.add_argument('--reference', help='Reference SMILES for Tanimoto similarity')
    parser.add_argument('--fp_type', choices=["morgan", "rdkit"], default="morgan", help='Fingerprint type')
    parser.add_argument('--radius', type=int, default=2, help='Morgan fingerprint radius (default: 2)')
    parser.add_argument('--bits', type=int, default=2048, help='Fingerprint bit size (default: 2048)')

    args = parser.parse_args()

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
