from rdkit import Chem
import dimorphite_dl
print(" Adapted from  Cheminform 11:14. doi:10.1186/s13321-019-0336-9")
print("git@github.com:UnixJunkie/dimorphite_dl.git")
# Using the dimorphite_dl.run() function, you can run Dimorphite-DL exactly as
# you would from the command line. Here's an example:
dimorphite_dl.run(
   smiles="",
   min_ph=7.4,
   max_ph=7.4,
   output_file="output.smi"
)
print("Output of first test saved to output.smi...")

# Using the dimorphite_dl.run_with_mol_list() function, you can also pass a
# list of RDKit Mol objects. The first argument is always the list.
# This part will pushed to main program smile2dock.py in the upcoming release.
smiles = ["C[C@](F)(Br)CC(O)=O", "CCCCCN"]
mols = [Chem.MolFromSmiles(s) for s in smiles]
for i, mol in enumerate(mols):
    mol.SetProp("msg","Orig SMILES: " + smiles[i])

protonated_mols = dimorphite_dl.run_with_mol_list(
    mols,
    min_ph=7.4,
    max_ph=7.4,
)
print([Chem.MolToSmiles(m) for m in protonated_mols])

# Note that properties are preserved.
print([m.GetProp("msg") for m in protonated_mols])

