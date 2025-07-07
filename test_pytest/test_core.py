import numpy as np

import sys
import os

# Add the parent directory (containing differintP) to the import path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from differintP.core import *  # type: ignore
from differintP.functions import MittagLeffler
from .__init__ import test_N

# Define constants to be used in tests.
size_coefficient_array = 20
sqrtpi2 = 0.88622692545275794
truevaluepoly = 0.94031597258
truevaluepoly_caputo = 1.50450555  # 8 / (3 * np.sqrt(np.pi))
truevaluepoly_caputo_higher = 2 / Gamma(1.5)

# Get SQRT results for checking accuracy.
GL_r = GL(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
GL_result = GL_r[-1]
GL_length = len(GL_r)

GLI_r = GLI(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
GLI_result = GLI_r[-1]
GLI_length = len(GLI_r)

RL_r = RL(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
RL_result = RL_r[-1]
RL_length = len(RL_r)

# --- Exponential function ---
alpha_exp = 0.5
groundtruth_exp = MittagLeffler(1, 1 - 0.5, np.array([1]))[0]

# --- Sine function ---
alpha_sin = 0.5
groundtruth_sin = np.sin(1 + alpha_sin * np.pi / 2)

#######################
# GLpoint
#######################

def test_GLpoint_sqrt_accuracy():
    assert abs(GLpoint(0.5, lambda x: x**0.5, 0.0, 1.0, 1024) - sqrtpi2) <= 1e-3

def test_GLpoint_accuracy_polynomial():
    assert abs(GLpoint(0.5, lambda x: x**2 - 1, 0.0, 1.0, 1024) - truevaluepoly) <= 1e-3

def test_GLpoint_accuracy_exp():
    val = GLpoint(alpha_exp, np.exp, 0, 1, 1024)
    assert abs(val - groundtruth_exp) < 1e-3

#######################
# GL
#######################

def test_GL_accuracy_sqrt():
    assert abs(GL_result - sqrtpi2) <= 1e-4

def test_GL_accuracy_polynomial():
    assert abs(GL(0.5, lambda x: x**2 - 1, 0.0, 1.0, 1024)[-1] - truevaluepoly) <= 1e-3

def test_GL_accuracy_exp():
    val = GL(alpha_exp, np.exp, 0, 1, 1024)[-1]
    assert abs(val - groundtruth_exp) < 1e-3

#######################
# GLI
#######################

def test_GLI_accuracy_sqrt():
    assert abs(GLI_result - sqrtpi2) <= 1e-4

def test_GLI_accuracy_polynomial():
    assert abs(GLI(0.5, lambda x: x**2 - 1, 0.0, 1.0, 1024)[-1] - truevaluepoly) <= 6e-3 # low accuracy

def test_GLI_accuracy_exp():
    val = GLI(alpha_exp, np.exp, 0, 1, 1024)[-1]
    assert abs(val - groundtruth_exp) < 5e-3 # low accuracy

#######################
# RLpoint
#######################

def test_RLpoint_sqrt_accuracy():
    assert abs(RLpoint(0.5, lambda x: x**0.5, 0.0, 1.0, 1024) - sqrtpi2) <= 1e-3

def test_RLpoint_accuracy_polynomial():
    assert abs(RLpoint(0.5, lambda x: x**2 - 1, 0.0, 1.0, 1024) - truevaluepoly) <= 1e-2

def test_RLpoint_accuracy_exp():
    val = RLpoint(alpha_exp, np.exp, 0, 1, 1024)
    assert abs(val - groundtruth_exp) < 1e-3

#######################
# RL
#######################

def test_RL_accuracy_sqrt():
    assert abs(RL_result - sqrtpi2) <= 1e-4

def test_RL_accuracy_polynomial():
    assert abs(RL(0.5, lambda x: x**2 - 1, 0.0, 1.0, 1024)[-1] - truevaluepoly) <= 1e-3

def test_RL_accuracy_exp():
    val = RL(alpha_exp, np.exp, 0, 1, 1024)[-1]
    assert abs(val - groundtruth_exp) < 1e-3

#######################
# Caputo 1p
#######################

def test_CaputoL1point_accuracy_sqrt():
    assert abs(CaputoL1point(0.5, lambda x: x**0.5, 0, 1.0, 1024) - sqrtpi2) <= 1e-2

def test_CaputoL1point_accuracy_polynomial():
    assert abs(CaputoL1point(0.5, lambda x: x**2 - 1, 0, 1.0, 1024) - truevaluepoly_caputo) <= 1e-3

def test_CaputoL1point_accuracy_exp():
    val = CaputoL1point(alpha_exp, np.exp, 0, 1, 1024)
    assert abs(val - groundtruth_exp) < 0.6 # really bad accuracy

#######################
# Caputo 2p
#######################

def test_CaputoL2point_accuracy_polynomial():
    assert abs(CaputoL2point(1.5, lambda x: x**2, 0, 1.0, 1024) - truevaluepoly_caputo_higher) <= 1e-1

#######################
# Caputo 2pC
#######################

def test_CaputoL2Cpoint_accuracy_polynomial_higher():
    assert abs(CaputoL2Cpoint(1.5, lambda x: x**2, 0, 1.0, 1024) - truevaluepoly_caputo_higher) <= 1e-1

def test_CaputoL2Cpoint_accuracy_polynomial():
    assert abs(CaputoL2Cpoint(0.5, lambda x: x**2, 0, 1.0, 1024) - truevaluepoly_caputo) <= 1e-3

def test_CaputoL2Cpoint_accuracy_exp():
    val = CaputoL2Cpoint(alpha_exp, np.exp, 0, 1, 1024)
    assert abs(val - groundtruth_exp) < 2 # unusable

#######################
# General algorithm output checks
#######################

def test_GL_result_length():
    assert GL_length == test_N

def test_GLI_result_length():
    assert GLI_length == test_N

def test_RL_result_length():
    assert RL_length == test_N

def test_RL_matrix_shape():
    assert np.shape(RLmatrix(0.4, test_N)) == (test_N, test_N)

def test_GL_binomial_coefficient_array_size():
    assert len(GLcoeffs(0.5, size_coefficient_array)) - 1 == size_coefficient_array

# Optional: Run doctest if called directly
if __name__ == "__main__":
    import doctest
    doctest.testmod()
