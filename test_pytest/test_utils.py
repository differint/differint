import numpy as np
import pytest

import sys
import os

# Add the parent directory (containing differintP) to the import path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from differintP.utils import *  # type: ignore
from .__init__ import test_N

stepsize = 1 / (test_N - 1)

# Testing if callable functions and arrays of function values will work.
checked_function1, test_stepsize1 = functionCheck(
    lambda x: 2 * np.exp(3 * x) * x - x**2 + x - 5, 0, 1, test_N
)
checked_function2, test_stepsize2 = functionCheck(np.ones(test_N), 0, 1, test_N)

def test_isInteger():
    assert isInteger(1)
    assert isInteger(1.0)
    assert isInteger(1 + 0j)
    assert not isInteger(1.1)
    assert not isInteger(1.1 + 0j)
    assert not isInteger(1 + 1j)

def test_isPositiveInteger():
    assert isPositiveInteger(1)
    assert not isPositiveInteger(1.1)
    assert not isPositiveInteger(-1)

def test_functionCheck():
    assert len(checked_function1) == test_N
    assert len(checked_function2) == test_N

    # Make sure it treats defined functions and arrays of function values the same.
    assert len(checked_function1) == len(checked_function2)
    assert test_stepsize1 == stepsize
    assert test_stepsize2 == stepsize
    assert test_stepsize1 == test_stepsize2

def test_checkValues():
    with pytest.raises(AssertionError):
        checkValues(0.1, 0, 1, 1.1) # type: ignore
    with pytest.raises(AssertionError):
        checkValues(0.1, 1j, 2, 100) # type: ignore
    with pytest.raises(AssertionError):
        checkValues(0.1, 1, 2j, 100) # type: ignore
    with pytest.raises(AssertionError):
        checkValues(0.1, 0, 1, -100)
    with pytest.raises(AssertionError):
        checkValues(1 + 1j, 0, 1, 100) # type: ignore
    # These should not raise
    checkValues(0.5, 0, 1, 100, support_complex_alpha=True)
    checkValues(1 + 1j, 0, 1, 100, support_complex_alpha=True) # type: ignore
    alpha_vals = np.array([0.1, 0.2])
    domain_vals = np.array([0.1, 1, 2.0, -1])
    num_vals = np.array([1.0, 100.0])
    [
        [
            [
                [
                    checkValues(alpha, domain_start, domain_end, num_points)
                    for alpha in alpha_vals
                ]
                for domain_start in domain_vals
            ]
            for domain_end in domain_vals
        ]
        for num_points in num_vals
    ]
