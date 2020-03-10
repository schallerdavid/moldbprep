from moldbprep.io import count_sdf_mols, sdf_properties
import os


def test_count_sdf_mols():
    assert count_sdf_mols(os.path.join(os.getcwd(), os.pardir, "data", "sample.sdf")) == 3

def test_sdf_properties():
    assert sdf_properties(os.path.join(os.getcwd(), os.pardir, "data", "sample.sdf")) == ['ID', 'vendor']
