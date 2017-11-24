import numpy as np
import math

def checkvalues(alpha, f_name, domain_start, domain_end, num_points):
    """ Checks for valid inputs."""
    if domain_end < domain_start:
        raise ValueError('Right endpoint must be larger than the left endpoint.')
        
    if not hasattr(f_name, '__call__'):
        msg = 'Input arg {} is not callable! Please input a function handle.'.format(str(f_name))
        raise TypeError(msg)
    
    # Set the default domain left-endpoint and number of domain points.
    if domain_start == None:
        domain_start = 0.
        
    if domain_end == None:
        domain_end = 1.
    
    if num_points == None:
        num_points = 100
    num_points = math.ceil(num_points)    
        
    return alpha, f_name, domain_start, domain_end, num_points

def GL(alpha, f_name, domain_start = None, domain_end = None, num_points = None):
    """Computes the GL fractional derivative of a function at a point.
       
       Parameters
       ==========
        alpha : float
            The order of the differintegral to be computed.
        f_name : function handle
            Either a lambda function or function handle. This is the function 
            that is to be differintegrated.
        domain_start : float
            The left-endpoint of the function domain.
        domain_end : float
            The right-endpoint of the function domain; the point at which the 
            differintegral is being evaluated.
        num_points : integer
            The number of points in the domain.
            
        Examples:
        >>> DF_sqrt = GL(-0.5, lambda x: np.sqrt(x))
        >>> DF_sqrt = GL(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> DF_sqrt = GL(0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    alpha, f_name, domain_start, domain_end, num_points = checkvalues(alpha, f_name, domain_start, domain_end, num_points)
    
    # Define the function domain.
    x = np.linspace(domain_start, domain_end, num_points + 1)
        
    # Calculate the GL differintegral, avoiding the explicit calculation of
    # the gamma function.
    GL_previous = f_name(x[1])
    for index in range(1,num_points):
        GL_current = GL_previous*(num_points-alpha-index-1)/(num_points-index) + f_name(x[index+1])
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
        f_name : function handle
            Either a lambda function or function handle. This is the function 
            that is to be differintegrated.
        domain_start : float
            The left-endpoint of the function domain.
        domain_end : float
            The right-endpoint of the function domain; the point at which the 
            differintegral is being evaluated.
        num_points : integer
            The number of points in the domain.
            
        Examples:
        >>> RL_sqrt = RLtrap(0.5, lambda x: np.sqrt(x))
        >>> RL_sqrt = RLtrap(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> RL_sqrt = RLtrap(-0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    alpha, f_name, domain_start, domain_end, num_points = checkvalues(alpha, f_name, domain_start, domain_end, num_points)
    
    C = 1/math.gamma(2-alpha)
    x = np.linspace(domain_start, domain_end, num_points+1)
    step_size = x[1] - x[0]
    
    RL = 0
    for index_j in range(num_points+1):
        coeff = RLcoeffs(num_points, index_j, alpha)
        RL += coeff*f_name(x[index_j])
        
    RL *= C*step_size**-alpha
    return RL