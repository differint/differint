# File for computation of errors for respective algorithms
import numpy as np

from differint.differint.differint import *

def CaputoL1Error(alpha, f_name, domain_start=0, domain_end=1, num_points=100):
    ''' Calculate a bound on the error from applying the CaputoL1 algorithm on a function.
        inf_norm(true - calculated) <= C * inf_norm(f'') * h ** (2 - alpha)
        where C = (alpha + 1)/Gamma(1 - alpha) * k^(2 - alpha)

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

def RLError(alpha, f_name, domain_start=0, domain_end=1, num_points=100):
    ''' Calculate a bound on the error from applying the RL Differintegral to a function. The 
        formula is accurate for a trapizoid approximation for numerical integration.
        |E| <= C_alpha * max(|f''|) * h^2 * b^alpha
        where b is the right endpoint, C_alpha is some constant dependant on alpha, and h is the spacing between points.

   	see Karniadakis, G.E.. (2019). Handbook of Fractional Calculus with Applications
    	Volume 3: Numerical Methods. De Gruyter.
    '''
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    norm_f = np.max(np.abs(np.diff(f_values, 2) / step_size ** 2))
    norm_f = np.nan_to_num(norm_f, nan=0)
    # some const "that depends only on alpha"... but what is it??
    C_alpha = 2 + alpha # guess
    return C_alpha * norm_f * domain_end ** alpha * step_size ** 2