import numpy as np

import sys
import os

# Add the parent directory (containing differintP) to the import path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from differintP.functions import *  # type: ignore

# Define constants to be used in tests.
poch_first_argument = 1
poch_second_argument = 5
poch_true_answer = 120

def test_pochhammer():
    assert poch(poch_first_argument, poch_second_argument) == poch_true_answer
    assert poch(-1, 3) == 0
    # assert poch(-1.5, 0.5) == np.inf   # commented as in original
    assert np.round(poch(1j, 1), 3) == 0.000 + 1.000j
    assert poch(-10, 2) == 90

def test_ML_cosh_root():
    xs = np.arange(0.1, 10, 0.1)
    assert np.all(
        np.abs(
            MittagLeffler(2, 1, xs, ignore_special_cases=True) - np.cosh(np.sqrt(xs))
        ) <= 1e-3
    )

def test_ML_exp():
    xs = np.arange(0.1, 10, 0.1)
    assert np.all(
        np.abs(
            MittagLeffler(1, 1, xs, ignore_special_cases=True) - np.exp(xs)
        ) <= 1e-3
    )

def test_ML_geometric():
    xs = np.arange(0.05, 0.95, 0.05)
    assert np.all(
        np.abs(
            MittagLeffler(0, 1, xs, ignore_special_cases=True) - 1 / (1 - xs)
        ) <= 1e-1
    )

# Optional: Run doctest if called directly
if __name__ == "__main__":
    import doctest
    doctest.testmod()
