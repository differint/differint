import numpy as np
import math

def GL(alpha, f_name, domain_start = 0., domain_end = 1., num_points = 100):
    """Computes the GL fractional derivative of a function at a point.
       
       Parameters
       ==========
        alpha : float
            The order of the differintegral to be computed.
        f_name : function handle, lambda function, list, or 1d-array of 
                 function values
            This is the function that is to be differintegrated.
        domain_start : float
            The left-endpoint of the function domain. Default value is 0.
        domain_end : float
            The right-endpoint of the function domain; the point at which the 
            differintegral is being evaluated. Default value is 1.
        num_points : integer
            The number of points in the domain. Default value is 100.
            
        Examples:
        >>> DF_sqrt = GL(-0.5, lambda x: np.sqrt(x))
        >>> DF_sqrt = GL(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> DF_sqrt = GL(0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    if domain_start > domain_end:
        msg = 'Domain start {} is greater than the domain end {}. Would you like to switch them? (y/n)'.format(str(domain_start),str(domain_end))
        print(msg)
        change_domain = input()
        if change_domain == 'y':
            domain_start, domain_end = domain_end, domain_start
        else:
            raise ValueError('Right endpoint must be greater than the left endpoint.')
    
        # Define the function domain and obtain function values.
    if hasattr(f_name, '__call__'):
        x = np.linspace(domain_start, domain_end, num_points + 1)
        f_values = list(map(lambda t: f_name(t), x))
    else:
        num_points = np.size(f_name)
        f_values = f_name
        
    # Calculate the GL differintegral, avoiding the explicit calculation of
    # the gamma function.
    GL_previous = f_values[1]
    for index in range(2,num_points):
        GL_current = GL_previous*(num_points-alpha-index-1)/(num_points-index) + f_values[index]
        GL_previous = GL_current
        
    return GL_current*(num_points/(domain_end - domain_start))**alpha

def RLcoeffs(index_k, index_j, alpha):
    """Calculates coefficients for the RL differintegral operator."""
    
    if index_j == 0:
        return ((index_k-1)**(1-alpha)-(index_k+alpha-1)*index_k**-alpha)
    elif index_j == index_k:
        return 1
    else:
        return ((index_k-index_j+1)**(1-alpha)+(index_k-index_j-1)**(1-alpha)-2*(index_k-index_j)**(1-alpha))
    
def RLtrap(alpha, f_name, domain_start = None, domain_end = None, num_points = None):
    """Calculate the RL differintegral with the trapezoid rule.
    
    Parameters
       ==========
        alpha : float
            The order of the differintegral to be computed.
        f_name : function handle, lambda function, list, or 1d-array of 
                 function values
            This is the function that is to be differintegrated.
        domain_start : float
            The left-endpoint of the function domain. Default value is 0.
        domain_end : float
            The right-endpoint of the function domain; the point at which the 
            differintegral is being evaluated. Default value is 1.
        num_points : integer
            The number of points in the domain. Default value is 100.
            
        Examples:
        >>> RL_sqrt = RLtrap(0.5, lambda x: np.sqrt(x))
        >>> RL_sqrt = RLtrap(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> RL_sqrt = RLtrap(-0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    if domain_start > domain_end:
        msg = 'Domain start {} is greater than the domain end {}. Would you like to switch them? (y/n)'.format(str(domain_start),str(domain_end))
        print(msg)
        change_domain = input()
        if change_domain == 'y':
            domain_start, domain_end = domain_end, domain_start
        else:
            raise ValueError('Right endpoint must be greater than the left endpoint.')
    
    # Define the function domain and obtain function values.
    if hasattr(f_name, '__call__'):
        x = np.linspace(domain_start, domain_end, num_points + 1)
        f_values = list(map(lambda t: f_name(t), x))
        step_size = x[1] - x[0]
    else:
        num_points = np.size(f_name) - 1
        f_values = f_name
        step_size = (domain_end - domain_start)/num_points
    
    C = 1/math.gamma(2-alpha)
    
    RL = 0
    for index_j in range(num_points+1):
        coeff = RLcoeffs(num_points, index_j, alpha)
        RL += coeff*f_values[index_j]
        
    RL *= C*step_size**-alpha
    return RL