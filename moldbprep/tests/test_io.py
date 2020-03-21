from moldbprep.io import count_sdf_mols, sdf_properties, time_to_text, sdf_text
import pytest
import os
from rdkit import Chem


def test_count_sdf_mols():
    assert count_sdf_mols(os.path.join(os.getcwd(), "moldbprep", "data", "db1.sdf")) == 3


def test_sdf_properties():
    assert sdf_properties(os.path.join(os.getcwd(), "moldbprep", "data", "db1.sdf")) == ['ID', 'vendor']


@pytest.mark.parametrize("time, text", [
    (5, '5 s'),
    (61, '1.0 min'),
    (3601, '1.0 h'),
    (86401, '1.0 d'),
    (1209601, '1.0 weeks'),
    (62899253, 'years'),
])
def test_time_to_text(time, text):
    assert time_to_text(time) == text


def test_sdf_text():
    assert sdf_text(Chem.MolFromSmiles('CCC'), {'db1': '1', 'db2': ''}) == \
           '\n' \
           '     RDKit          2D\n' \
           '\n' \
           '  3  2  0  0  0  0  0  0  0  0999 V2000\n' \
           '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
           '    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
           '    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n' \
           '  1  2  1  0\n' \
           '  2  3  1  0\n' \
           'M  END\n' \
           '>  <db1>\n' \
           '1\n' \
           '\n' \
           '>  <db2>\n' \
           '\n' \
           '\n' \
           '$$$$\n'
