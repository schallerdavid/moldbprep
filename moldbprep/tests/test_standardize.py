from moldbprep.standardize import largest_fragment, protonate_mol
import pytest
from rdkit import Chem


@pytest.mark.parametrize("input_smiles, output_smiles", [
    ("[Cl-].[NH3+]CCC1=CC=CC=C1", "[NH3+]CCC1=CC=CC=C1"),
    ("[O-]C=O.[NH3+]CCC1=CC=CC=C1", "[NH3+]CCC1=CC=CC=C1")
])
def test_largest_fragment(input_smiles, output_smiles):
    assert Chem.MolToSmiles(largest_fragment(Chem.MolFromSmiles(input_smiles))) == \
           Chem.MolToSmiles(Chem.MolFromSmiles(output_smiles))


@pytest.mark.parametrize("input_smiles, output_smiles", [
    ("CNC1CCN(CCCN)CC1", "C[NH2+]C1CC[NH+](CCC[NH3+])CC1"),
    ("CC[O-]", "CCO"),
    ("CNC(N)=N", "CNC(N)=[NH2+]"),
    ("CCC(N)=N", "CCC(N)=[NH2+]"),
    ("CCC(O)=O", "CCC([O-])=O"),
    ("CCP(O)(O)=O", "CCP([O-])([O-])=O"),
    ("CS(=O)(=O)NC=O", "CS(=O)(=O)[N-]C=O"),
    ("CC1=NNN=N1", "CC1=N[N-]N=N1"),
    ("CC(=O)C=CO", "CC(=O)C=C[O-]")
])
def test_protonate_mol(input_smiles, output_smiles):
    assert Chem.MolToSmiles(protonate_mol(Chem.MolFromSmiles(input_smiles))) == \
           Chem.MolToSmiles(Chem.MolFromSmiles(output_smiles))
