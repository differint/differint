import unittest
import numpy as np
import differint

class GLTestCase(unittest.TestCase):
    """Test for GL"""
    
    def test_GL_sqrt_x(self):
        self.assertTrue(abs(GL(0.5,f,0.,1.,100)-np.sqrt(np.pi)/2) <= 1e-3)
        
    def test_GL_endpoints(self):
        self.assertRaises(ValueError, GL, 0.5, f, 1., 0., 100)
        
    def test_GL_arg_is_func_handle(self):
        self.assertRaises(TypeError, GL, 0.5, "O hai Mark", 0., 1., 100)
        
if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
