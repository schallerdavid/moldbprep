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
    ("CC(=O)C=CO", "CC(=O)C=C[O-]"),
    ("c1(C2CCCCC2)[nH]nnn1", "c1(C2CCCCC2)nnn[n-]1"),
    ("c1(C2CCCCC2)n[nH]nn1", "c1(C2CCCCC2)nn[n-]n1"),
    ("c1(C2CCCCC2)nnn[n-]1", "c1(C2CCCCC2)nnn[n-]1"),
    ("c1(C2CCCCC2)nn[n-]n1", "c1(C2CCCCC2)nn[n-]n1"),
    ("O=C1C(O)=C(O)C([C@H](O)CO)O1", "O=C1C(O)=C([O-])C([C@H](O)CO)O1"),
    ("O=C1C(O)=C([O-])C([C@H](O)CO)O1", "O=C1C(O)=C([O-])C([C@H](O)CO)O1"),
    ("O=C(O)C1CCCCC1", "O=C([O-])C1CCCCC1"),
    ("S(=O)(=O)(NC(=O)NCCCC)c1ccc(C)cc1", "S(=O)(=O)([N-]C(=O)NCCCC)c1ccc(C)cc1"),
    ("S(=O)(=O)(Nc1noc(C)c1)c1ccc(N)cc1", "S(=O)(=O)(Nc1noc(C)c1)c1ccc(N)cc1"),
    ("P(=O)(O)(O)[C@@H]1[C@H](C)O1", "P(=O)([O-])([O-])[C@@H]1[C@H](C)O1"),
    ("S(=O)(=O)(O)c1ccc(C)cc1", "S(=O)(=O)([O-])c1ccc(C)cc1"),
    ("CNC(=O)Oc1ccc2N(C)[C@H]3N(C)CC[C@@]3(C)c2c1", "CNC(=O)Oc1ccc2N(C)[C@H]3[NH+](C)CC[C@@]3(C)c2c1"),
    ("S(=O)(=O)(NC(=N)N)c1ccc(N)cc1", "S(=O)(=O)(NC(=N)N)c1ccc(N)cc1"),
    ("S(CCN/C(=N/C#N)/NC)Cc1c(C)nc[nH]1", "S(CCN/C(=N/C#N)/NC)Cc1c(C)nc[nH]1"),
    ("Clc1cc2nccc(N[C@H](CCCN(CC)CC)C)c2cc1", "Clc1cc2nccc(N[C@H](CCC[NH+](CC)CC)C)c2cc1"),
    ("O=C(O)[C@@H](N)Cc1nc[nH]c1", "O=C([O-])[C@@H]([NH3+])Cc1nc[nH]c1"),
    ("FC1(F)C(N)CCCC1", "FC1(F)C(N)CCCC1"),
    ("[C@H](CN1C[C@@H](C)O[C@@H](C)C1)(Cc1ccc(C(CC)(C)C)cc1)C",
     "[C@H](CN1C[C@@H](C)O[C@@H](C)C1)(Cc1ccc(C(CC)(C)C)cc1)C"),
    ("O[C@@H](CNC)c1cc(O)c(O)cc1", "O[C@@H](C[NH2+]C)c1cc(O)c(O)cc1"),
    ("O(C)c1c(OC)cc(Cc2c(N)nc(N)nc2)cc1OC", "O(C)c1c(OC)cc(Cc2c(N)nc(N)nc2)cc1OC"),
    ("O=C(NN)c1ccncc1", "O=C(NN)c1ccncc1")
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
