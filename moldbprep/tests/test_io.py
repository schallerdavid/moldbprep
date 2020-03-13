from moldbprep.io import count_sdf_mols, sdf_properties, time_to_text
import pytest
import os


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
