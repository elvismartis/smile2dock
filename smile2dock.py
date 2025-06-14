import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from openbabel import pybel

def smiles_to_3d(smiles, output_base="molecule", num_confs=10, optimize=True):
    """Convert SMILES to multiple 3D formats and calculate properties"""
    try:
        # RDKit: Parse and add hydrogens
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            return None, None
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
            "Molecular Weight": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "H-Bond Donors": Descriptors.NumHDonors(mol),
            "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
            "TPSA": Descriptors.TPSA(mol),
            "Rotatable Bonds": Descriptors.NumRotatableBonds(mol)
        }
        return mol, properties

    except Exception as e:
        print(f"Error processing {smiles}: {str(e)}")
        return None, None

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

def batch_process(input_file, output_dir, reference_smiles=None, fp_type="morgan", radius=2, n_bits=2048):
    """Process a file of SMILES strings, optionally computing similarity to a reference"""
    if reference_smiles:
        ref_mol = Chem.MolFromSmiles(reference_smiles)
        if ref_mol is None:
            print(f"Invalid reference SMILES: {reference_smiles}")
            ref_mol = None
    else:
        ref_mol = None

    with open(input_file) as f:
        for idx, line in enumerate(f):
            smiles = line.strip()
            if smiles:
                base_name = os.path.join(output_dir, f"mol_{idx+1}")
                mol, props = smiles_to_3d(smiles, base_name)
                if props:
                    print(f"Processed {smiles}")
                    print("Properties:", props)
                    if ref_mol:
                        similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                        print(f"Tanimoto similarity to reference: {similarity:.4f}")

def single_process(smiles, output_base, num_confs, reference_smiles=None, fp_type="morgan", radius=2, n_bits=2048):
    """Process a single SMILES string, optionally computing similarity to a reference"""
    mol, props = smiles_to_3d(smiles, output_base, num_confs)
    if props:
        print("Generated files:", [f"{output_base}.{fmt}" for fmt in ["pdb", "mol2", "sdf", "pdbqt"]])
        print("\nMolecular Properties:")
        for k, v in props.items():
            print(f"{k}: {v:.2f}")
        if reference_smiles:
            ref_mol = Chem.MolFromSmiles(reference_smiles)
            if ref_mol:
                similarity = calculate_tanimoto_similarity(mol, ref_mol, fp_type, radius, n_bits)
                print(f"Tanimoto similarity to reference: {similarity:.4f}")
            else:
                print(f"Invalid reference SMILES: {reference_smiles}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SMILES to 3D converter with property and similarity calculation')
    parser.add_argument('input', help='SMILES string or input file')
    parser.add_argument('-o', '--output', default="output", help='Output directory or base name')
    parser.add_argument('-n', '--num_confs', type=int, default=10, help='Number of conformers to generate')
    parser.add_argument('--reference', help='Reference SMILES for Tanimoto similarity')
    parser.add_argument('--fp_type', choices=["morgan", "rdkit"], default="morgan", help='Fingerprint type for similarity')
    parser.add_argument('--radius', type=int, default=2, help='Morgan fingerprint radius (default: 2)')
    parser.add_argument('--bits', type=int, default=2048, help='Fingerprint bit size (default: 2048)')

    args = parser.parse_args()

    if os.path.isfile(args.input):
        os.makedirs(args.output, exist_ok=True)
        batch_process(args.input, args.output, args.reference, args.fp_type, args.radius, args.bits)
    else:
        single_process(args.input, args.output, args.num_confs, args.reference, args.fp_type, args.radius, args.bits)
