from moldbprep.standardize import protonate_mol, enumerate_stereo_isomers, merge_ids
import pandas as pd
import pytest
from rdkit import Chem


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


def test_enumerate_stereo_isomers():
    assert sorted([Chem.MolToSmiles(mol) for mol in enumerate_stereo_isomers(Chem.MolFromSmiles('FC(Cl)Br'), 8)]) == \
        sorted(['F[C@H](Cl)Br', 'F[C@@H](Cl)Br'])
    assert len(enumerate_stereo_isomers(Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1'), 4)) == 4


def test_merge_ids():
    assert merge_ids(pd.DataFrame([['C', '1', ''],
                                   ['C', '2', ''],
                                   ['A', '3', ''],
                                   ['C', '', '1'],
                                   ['B', '', '2']], columns=['smiles', 'DB1', 'DB2']), ['DB1', 'DB2']).equals(
           pd.DataFrame([['A', '3', ''],
                         ['B', '', '2'],
                         ['C', '1,2', '1']], columns=['smiles', 'DB1', 'DB2']))
