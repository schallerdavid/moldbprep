from moldbprep.standardize import standardize_mol
import pytest
from rdkit import Chem


@pytest.mark.parametrize("input_smiles, output_smiles", [
    ("[Cl-].CNC1CCN(CCCN)CC1", "C[NH2+]C1CC[NH+](CCC[NH3+])CC1"),
    ("CC[O-]", "CCO"),
    ("CNC(N)=N", "CNC(N)=[NH2+]"),
    ("CCC(N)=N", "CCC(N)=[NH2+]"),
    ("CCC(O)=O", "CCC([O-])=O"),
    ("CCP(O)(O)=O" ,"CCP([O-])([O-])=O"),
    ("CS(=O)(=O)NC=O", "CS(=O)(=O)[N-]C=O"),
    ("CC1=NNN=N1", "CC1=N[N-]N=N1"),
    ("CC(=O)C=CO", "CC(=O)C=C[O-]")
])
def test_standardize_mol(input_smiles, output_smiles):
    assert standardize_mol(Chem.MolFromSmiles(input_smiles)) == Chem.MolToSmiles(Chem.MolFromSmiles(output_smiles))
