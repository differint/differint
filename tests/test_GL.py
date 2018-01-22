import unittest
import numpy as np
import differint

class GLTestCase(unittest.TestCase):
    """ Tests for algorithm accuracy. """
    
    def test_GLpoint_sqrt_accuracy(self):
        self.assertTrue(abs(GLpoint(0.5,lambda x: x**0.5,0.,1.,100)-np.sqrt(np.pi)/2) <= 1e-3)
    
    def test_GLpoint_accuracy_polynomial(self):
        self.assertTrue(abs(GLpoint(0.5,lambda x: x**2-1,0.,1.,100)-0.94031597258) <= 1e-3)
        
    def test_RLpoint_sqrt_accuracy(self):
        self.assertTrue(abs(RLtrap(0.5,lambda x: x**0.5,0.,1.,100)-np.sqrt(np.pi)/2) <= 1e-3)
        
    def test_RLpoint_accuracy_polynomial(self):
        self.assertTrue(abs(RLtrap(0.5,lambda x: x**2-1,0.,1.,100)-0.94031597258) <= 1e-3)
        
if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    
    # Ensure all docstring examples work.
    import doctest
    doctest.testmod()
