import unittest
import numpy as np
import differint

# Define constants to be used in tests.
poch_first_argument = 1
poch_second_argument = 5
poch_true_answer = 120
size_matrix = 15
size_coefficient_array = 20

test_N = 512

sqrtpi2 = 0.88622692545275794
truevaluepoly = 0.94031597258

INTER = GLIinterpolat(1)

stepsize = 1/(test_N-1)
checked_function, test_stepsize = functionCheck(test_func,0,1,test_N)

# Get results for checking accuracy.
GL_r = GL(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
GL_result = GL_r[-1]

GLI_r = GLI(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
GLI_result = GLI_r[-1]
GLI_length = len(GLI_r)

RL_r = RL(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
RL_result = RL_r[-1]
RL_length = len(RL_r)

class HelperTestCases(unittest.TestCase):
    """ Tests for helper functions. """
    
    def test_isInteger(self):
        with self.assertRaises(AssertionError):
            isInteger(1.1)
    
    def test_pochhammer(self):
        self.assertEqual(poch(poch_first_argument, poch_second_argument), poch_true_answer)
        
    def test_functionCheck(self):
        self.assertEqual(len(checked_function), test_N)
        self.assertEqual(test_stepsize, stepsize)
        
    def test_GL_binomial_coefficient_array_size(self):
        self.assertEqual(len(GLcoeffs(0.5,size_coefficient_array))-1,size_coefficient_array)
        
class TestInterpolantCoefficients(unittest.TestCase):
    """ Test the correctness of the interpolant coefficients. """
    def check_coefficients(self):
        self.assertEqual(INTER.prv, -0.125)
        self.assertEqual(INTER.crr, 0.75)
        self.assertEqual(INTER.nxt, 0.375)
    
class TestAlgorithms(unittest.TestCase):
    """ Tests for correct size of algorithm results. """
    
    def test_GLI_result_length(self):
        self.assertEqual(GLI_length,test_N)
    
    def test_RL_result_length(self):
        self.assertEqual(RL_length, test_N)
    
    """ Tests for algorithm accuracy. """
    
    def test_GLpoint_sqrt_accuracy(self):
        self.assertTrue(abs(GLpoint(0.5,lambda x: x**0.5,0.,1.,1024)-sqrtpi2) <= 1e-3)
    
    def test_GLpoint_accuracy_polynomial(self):
        self.assertTrue(abs(GLpoint(0.5,lambda x: x**2-1,0.,1.,1024)-truevaluepoly) <= 1e-3)
        
    def test_GL_accuracy_sqrt(self):
        self.assertTrue(abs(GL_result - sqrtpi2) <= 1e-4)
        
    def test_GLI_accuracy_sqrt(self):
        self.assertTrue(abs(GLI_result - sqrtpi2) <= 1e-4)
        
    def test_RLpoint_sqrt_accuracy(self):
        self.assertTrue(abs(RLpoint(0.5,lambda x: x**0.5,0.,1.,1024)-sqrtpi2) <= 1e-3)
        
    def test_RLpoint_accuracy_polynomial(self):
        self.assertTrue(abs(RLpoint(0.5,lambda x: x**2-1,0.,1.,1024)-truevaluepoly) <= 1e-2)
        
    def test_RL_accuracy_sqrt(self):
        self.assertTrue(abs(RL_result - sqrtpi2) <= 1e-4)
        
if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    
    # Ensure all docstring examples work.
    import doctest
    doctest.testmod()