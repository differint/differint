import unittest
import numpy as np

# Import from sibling directory.
from differint.differint import *

# Define constants to be used in tests.
poch_first_argument = 1
poch_second_argument = 5
poch_true_answer = 120
size_coefficient_array = 20
test_N = 512
sqrtpi2 = 0.88622692545275794
truevaluepoly = 0.94031597258
truevaluepoly_caputo = 1.50450555 # 8 / (3 * np.sqrt(np.pi))
truevaluepoly_caputo_higher = 2 / Gamma(1.5)
PC_x_power = np.linspace(0, 1, 100) ** 5.5

INTER = GLIinterpolat(1)

stepsize = 1/(test_N-1)

# Testing if callable functions and arrays of function values will work.
checked_function1, test_stepsize1 = functionCheck(lambda x: 2*np.exp(3*x)*x - x**2 + x - 5,0,1,test_N)
checked_function2, test_stepsize2 = functionCheck(np.ones(test_N),0,1,test_N)

# Get results for checking accuracy.
GL_r = GL(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
GL_result = GL_r[-1]

GLI_r = GLI(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
GLI_result = GLI_r[-1]
GLI_length = len(GLI_r)

RL_r = RL(0.5, lambda x: np.sqrt(x), 0, 1, test_N)
RL_result = RL_r[-1]
RL_length = len(RL_r)

# Get FODE function for solving.
PC_func_power = lambda x, y : 1/24 * Gamma(5 + 1.5) * x**4 + x**(8 + 2 * 1.5) - y**2
PC_func_ML = lambda x,y : y

class HelperTestCases(unittest.TestCase):
    """ Tests for helper functions. """
    
    def test_isInteger(self):
        self.assertTrue(isInteger(1))
        self.assertTrue(isInteger(1.0))
        self.assertTrue(isInteger(1+0j))
        self.assertFalse(isInteger(1.1))
        self.assertFalse(isInteger(1.1+0j))
        self.assertFalse(isInteger(1+1j))

    def test_isPositiveInteger(self):
        self.assertTrue(isPositiveInteger(1))
        self.assertFalse(isPositiveInteger(1.1))
        self.assertFalse(isPositiveInteger(-1))
    
    def test_pochhammer(self):
        self.assertEqual(poch(poch_first_argument, poch_second_argument), poch_true_answer)
        self.assertEqual(poch(-1, 3), 0)
        self.assertEqual(poch(-1.5, 0.5), np.inf)
        self.assertEqual(np.round(poch(1j, 1), 3), 0.000+1.000j)
        self.assertEqual(poch(-10, 2), 90)
        
    def test_functionCheck(self):
        self.assertEqual(len(checked_function1), test_N)
        self.assertEqual(len(checked_function2), test_N)
        
        # Make sure it treats defined functions and arrays of function values the same.
        self.assertEqual(len(checked_function1), len(checked_function2))
        self.assertEqual(test_stepsize1, stepsize)
        self.assertEqual(test_stepsize2, stepsize)
        self.assertEqual(test_stepsize1, test_stepsize2)
        
    def test_GL_binomial_coefficient_array_size(self):
        self.assertEqual(len(GLcoeffs(0.5,size_coefficient_array))-1,size_coefficient_array)
        
    def test_checkValues(self):
        with self.assertRaises(AssertionError):
            checkValues(0.1, 0, 1, 1.1)
        with self.assertRaises(AssertionError):
            checkValues(0.1, 1j, 2, 100)
        with self.assertRaises(AssertionError):
            checkValues(0.1, 1, 2j, 100)
        with self.assertRaises(AssertionError):
            checkValues(0.1, 0, 1, -100)
        with self.assertRaises(AssertionError):
            checkValues(1+1j, 0, 1, 100)
        checkValues(0.5, 0, 1, 100, support_complex_alpha=True)
        checkValues(1+1j, 0, 1, 100, support_complex_alpha=True)
        alpha_vals = np.array([0.1, 0.2])
        domain_vals = np.array([0.1, 1, 2.0, -1])
        num_vals = np.array([1., 100.0])
        [[[[checkValues(alpha, domain_start, domain_end, num_points) for alpha in alpha_vals] 
                                                                     for domain_start in domain_vals] 
                                                                     for domain_end in domain_vals]
                                                                     for num_points in num_vals]
            
    """ Unit tests for gamma function. """
    
    def testFiveFactorial(self):
        self.assertEqual(Gamma(6),120)
        
    def testNegativePoles(self):
        self.assertEqual(Gamma(-2),np.inf)
        
    def testRealValue(self):
        self.assertEqual(np.round(Gamma(1.25), 12), 0.906402477055)

    def testComplexValue(self):
        self.assertEqual(np.round(Gamma(1j), 4), -0.1549-0.498j)

    """ Unit tests for Mittag-Leffler function. """

    def test_ML_cosh_root(self):
        xs = np.arange(10, 0.1)
        self.assertTrue((np.abs(MittagLeffler(2, 1, xs, ignore_special_cases=True)\
                                        - np.cosh(np.sqrt(xs))) <= 1e-3).all())

    def test_ML_exp(self):
        xs = np.arange(10, 0.1)
        self.assertTrue((np.abs(MittagLeffler(1, 1, xs, ignore_special_cases=True)\
                                        - np.exp(xs)) <= 1e-3).all())

    def test_ML_geometric(self):
        xs = np.arange(1, 0.05)
        self.assertTrue((np.abs(MittagLeffler(0, 1, xs, ignore_special_cases=True)\
                                        - 1 / (1 - xs)) <= 1e-3).all())
        
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
        
    def test_RL_matrix_shape(self):
        self.assertTrue(np.shape(RLmatrix(0.4, test_N)) == (test_N, test_N))
    
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

    def test_CaputoL1point_accuracy_sqrt(self):
        self.assertTrue(abs(CaputoL1point(0.5, lambda x: x**0.5, 0, 1., 1024)-sqrtpi2) <= 1e-2)

    def test_CaputoL1point_accuracy_polynomial(self):
        self.assertTrue(abs(CaputoL1point(0.5, lambda x: x**2-1, 0, 1., 1024)-truevaluepoly_caputo) <= 1e-3)
    
    def test_CaputoL2point_accuracy_polynomial(self):
        self.assertTrue(abs(CaputoL2point(1.5, lambda x: x**2, 0, 1., 1024)-truevaluepoly_caputo_higher) <= 1e-1)

    def test_CaputoL2Cpoint_accuracy_polynomial_higher(self):
        self.assertTrue(abs(CaputoL2Cpoint(1.5, lambda x: x**2, 0, 1., 1024)-truevaluepoly_caputo_higher) <= 1e-1)

    def test_CaputoL2Cpoint_accuracy_polynomial(self):
        self.assertTrue(abs(CaputoL2Cpoint(0.5, lambda x: x**2, 0, 1., 1024)-truevaluepoly_caputo) <= 1e-3)

class TestSolvers(unittest.TestCase):
    """ Tests for the correct solution to the equations. """
    def test_PC_solution_three_halves(self):
        self.assertTrue((np.abs(PCsolver([0, 0], 1.5, PC_func_power, 0, 1, 100)-PC_x_power) <= 1e-2).all())

    def test_PC_solution_ML(self):
        xs = np.linspace(0, 1, 100)
        ML_alpha = MittagLeffler(5.5, 1, xs ** 5.5)
        self.assertTrue((np.abs(PCsolver([1, 0, 0, 0, 0, 0], 5.5, PC_func_ML)-ML_alpha) <= 1e-2).all())

    def test_PC_solution_linear(self):
        xs = np.linspace(0, 1, 100)
        self.assertTrue((np.abs(PCsolver([1, 1], 1.5, lambda x,y : y-x-1)-(xs+1)) <= 1e-2).all())

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    
    # Ensure all docstring examples work.
    import doctest
    doctest.testmod()
