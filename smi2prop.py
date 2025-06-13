# Example property calculation extension
from rdkit.Chem import rdMolDescriptors

def calculate_additional_properties(mol):
    return {
        "QED": rdMolDescriptors._qed.QED(mol),
        "SAScore": rdMolDescriptors.CalcSAScore(mol),
        "RingCount": rdMolDescriptors.CalcNumRings(mol),
        "AMR": rdMolDescriptors.CalcAMR(mol)
    }

