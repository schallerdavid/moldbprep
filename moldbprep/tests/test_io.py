from moldbprep.io import count_sdf_mols, sdf_properties


def test_count_sdf_mols():
    assert count_sdf_mols("../data/sample.sdf") == 3

def test_sdf_properties():
    assert sdf_properties("../data/sample.sdf") == ['ID', 'vendor']
