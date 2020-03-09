from moldbprep.io import count_sdf_mols


def test_count_sdf_mols():
    assert count_sdf_mols("../data/sample.sdf") == 3
