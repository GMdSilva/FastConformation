"""
Unit and regression test for the decaf_e_dev package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import decaf_e_dev


def test_decaf_e_dev_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "decaf_e_dev" in sys.modules
