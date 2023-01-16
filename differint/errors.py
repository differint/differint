# File for computation of errors for respective algorithms
import numpy as np

from differint.differint.differint import *

def find_local_lipschitz_const()

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
    ''' Calculate a bound on the error from applying the RL Differintegral to a function. The error
        calculation comes from the error on the trapizoid approximation for numerical integration,
        |E| <= (b-a)^3/(12n^2) * max(|f''|)
        where b is the right endpoint, a is the left endpoint, and n is the number of points.
    '''
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    norm_f = np.max(np.abs(np.diff(f_values, 2) / step_size ** 2))
    norm_f = np.nan_to_num(norm_f, nan=0)
    #print(norm_f)
    #print((domain_end - domain_start) ** 3 / (12 * num_points ** 2))
    #return (domain_end - domain_start) ** 3 / (12 * num_points ** 2) * norm_f
    # some const "that depends only on alpha"... but what is it??
    C_alpha = 1
    return C_alpha * norm_f * domain_end ** alpha * step_size ** 2