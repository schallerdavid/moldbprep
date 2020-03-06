"""
Unit and regression test for the moldbprep package.
"""

# Import package, test suite, and other packages as needed
import moldbprep
import pytest
import sys

def test_moldbprep_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "moldbprep" in sys.modules
