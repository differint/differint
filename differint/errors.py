# File for computation of errors for respective algorithms
import numpy as np

from differint.differint.differint import *

def CaputoL1Error(alpha, f_name, domain_start=0, domain_end=1, num_points=100):
	''' Calculate a bound on the error from applying the CaputoL1 algorithm on a function.
		inf_norm(true - calculated) <= C * inf_norm(f'') * h ** (2 - alpha)

	see Ming Li et al (2011). A numerical evaluation and regularization of Caputo fractional
	    derivatives. Journal of Physics: Conference Series.
	'''
	# Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

	norm_f = np.max(np.abs(np.diff(f_values, 2) / step_size ** 2))
	return (alpha + 1) / Gamma(1 - alpha) * norm_f * step_size ** (2 - alpha)