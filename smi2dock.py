from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from openbabel import pybel
import argparse
import os
from io import StringIO

def smiles_to_3d(smiles, output_base="molecule", num_confs=10, optimize=True):
    """Convert SMILES to multiple 3D formats and calculate properties"""
    try:
        # RDKit processing
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformers
        AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, randomSeed=42)
        
        if optimize:
            for conf_id in range(mol.GetNumConformers()):
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)

        # Convert to OpenBabel format
        sdf_data = Chem.MolToMolBlock(mol)
        ob_mol = pybel.readstring("mol", sdf_data)
        
        # Write different file formats for each conformer
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
        
        return properties

    except Exception as e:
        print(f"Error processing {smiles}: {str(e)}")
        return None

def batch_process(input_file, output_dir):
    """Process a file of SMILES strings"""
    with open(input_file) as f:
        for line in f:
            smiles = line.strip()
            if smiles:
                base_name = os.path.join(output_dir, smiles[:50])
                props = smiles_to_3d(smiles, base_name)
                if props:
                    print(f"Processed {smiles}")
                    print("Properties:", props)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SMILES to 3D converter')
    parser.add_argument('input', help='SMILES string or input file')
    parser.add_argument('-o', '--output', default="output", 
                      help='Output directory/base name')
    parser.add_argument('-n', '--num_confs', type=int, default=10,
                      help='Number of conformers to generate')
    
    args = parser.parse_args()
    
    if os.path.isfile(args.input):
        os.makedirs(args.output, exist_ok=True)
        batch_process(args.input, args.output)
    else:
        props = smiles_to_3d(args.input, args.output, args.num_confs)
        if props:
            print("Generated files:", 
                 [f"{args.output}.{fmt}" for fmt in ["pdb", "mol2", "sdf", "pdbqt"]])
            print("\nMolecular Properties:")
            for k, v in props.items():
                print(f"{k}: {v:.2f}")
