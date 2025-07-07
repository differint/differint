import numpy as np

import sys
import os

# Add the parent directory (containing differintP) to the import path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from differintP.core import Gamma
from differintP.solvers import PCsolver
from differintP.functions import MittagLeffler

PC_x_power = np.linspace(0, 1, 100) ** 5.5

# Get FODE function for solving.
PC_func_power = lambda x, y: 1 / 24 * Gamma(5 + 1.5) * x**4 + x ** (8 + 2 * 1.5) - y**2
PC_func_ML = lambda x, y: y

def test_PC_solution_three_halves():
    result = PCsolver([0, 0], 1.5, PC_func_power, 0, 1, 100)
    assert np.all(np.abs(result - PC_x_power) <= 1e-2)

def test_PC_solution_ML():
    xs = np.linspace(0, 1, 100)
    ML_alpha = MittagLeffler(5.5, 1, xs**5.5)
    result = PCsolver([1, 0, 0, 0, 0, 0], 5.5, PC_func_ML)
    assert np.all(np.abs(result - ML_alpha) <= 1e-2)

def test_PC_solution_linear():
    xs = np.linspace(0, 1, 100)
    result = PCsolver([1, 1], 1.5, lambda x, y: y - x - 1)
    assert np.all(np.abs(result - (xs + 1)) <= 1e-2)

# Optional: Run doctest (only if you want this to run on manual invocation)
if __name__ == "__main__":
    import doctest
    doctest.testmod()
