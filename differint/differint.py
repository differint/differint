import numpy as np
import math

class GLIinterpolat:
    """ Class for computing interpolation of function values for the 
        improved GL algorithm. 
        
        Using a class here helps avoid type flexibility for these constants.
        """
    
    def __init__(self,alpha):
        # Determine coefficients for quadratic interpolation.
        self.nxt = alpha*(2+alpha)/8
        self.crr = (4-alpha*alpha)/4
        self.prv = alpha*(alpha-2)/8

def isInteger(n):
    assert (n >= 0 and (type(n) is type(0))), "n must be a positive integer or zero: %r" % n

def checkValues(alpha, domain_start, domain_end, num_points):
    """ Type checking for valid inputs. """
    
    assert type(num_points) is type(1), "num_points is not an integer: %r" % num_points
    
    assert type(domain_start) is type(0.0) \
        or type(domain_start) is type(0), "domain_start must be integer or float: %r" % domain_start
        
    assert type(domain_end) is type(0.0) \
        or type(domain_end) is type(0), "domain_end must be integer or float: %r" % domain_end
        
    # Currently there is no support for complex orders (17 Jan 2017).
    assert type(alpha) is not type(1+1j), "alpha must be real: %r" % alpha
    
    return   

def functionCheck(f_name, domain_start, domain_end, num_points):
    """ Check if function is callable and assign function values. """
    
    # Define the function domain and obtain function values.
    if hasattr(f_name, '__call__'):
        # If f_name is callable, call it and save to a list.
        x = np.linspace(domain_start, domain_end, num_points)
        f_values = list(map(lambda t: f_name(t), x)) 
        step_size = x[1] - x[0]
    else:
        num_points = np.size(f_name)
        f_values = f_name
        step_size = (domain_end - domain_start)/(num_points-1)
    return f_values, step_size

def poch(a,n):
    """ Returns the Pochhammer symbol (a)_n. """
    
    # First, check if 'a' is a real number (this is currently only working for reals).
    assert type(a) is not type(1+1j), "a must be real: %r" % a
    isInteger(n)
    
    # Compute the Pochhammer symbol.
    if n == 0:
        return 1.0
    else:
        poch = 1
        for j in range(n):
            poch *= (a + j)
        return poch
    
def GLcoeffs(alpha,n):
    """ Computes the GL coefficient array of size n. 
    
        These coefficients can be used for both the GL 
        and the improved GL algorithm.
    """ 
    
    # Validate input.
    isInteger(n)
    
    # Get generalized binomial coefficients.
    GL_filter = np.zeros(n+1,)
    GL_filter[0] = 1
    
    for i in range(n):
        GL_filter[i+1] = GL_filter[i]*(-alpha + i)/(i+1)
    
    return GL_filter

def GLpoint(alpha, f_name, domain_start = 0., domain_end = 1., num_points = 100):
    """ Computes the GL fractional derivative of a function at a point.
       
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
        >>> DF_sqrt = GLpoint(-0.5, lambda x: np.sqrt(x))
        >>> DF_sqrt = GLpoint(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> DF_sqrt = GLpoint(0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
        
    # Calculate the GL differintegral, avoiding the explicit calculation of
    # the gamma function.
    GL_previous = f_values[1]
    for index in range(2,num_points):
        GL_current = GL_previous*(num_points-alpha-index-1)/(num_points-index) + f_values[index]
        GL_previous = GL_current
        
    return GL_current*(num_points/(domain_end - domain_start))**alpha

def GL(alpha, f_name, domain_start = 0.0, domain_end = 1.0, num_points = 100):
    """ Computes the GL fractional derivative of a function for an entire array
        of function values.
        
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
        >>> DF_poly = GL(-0.5, lambda x: x**2 - 1)
        >>> DF_sqrt = GL(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> DF_sqrt = GL(0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
        
    # Get the convolution filter.
    b_coeffs = GLcoeffs(alpha, num_points)
    
    # Real Fourier transforms for convolution filter and array of function values.
    B = np.fft.rfft(b_coeffs)
    F = np.fft.rfft(f_values)
    
    result = np.fft.irfft(F*B)*step_size**-alpha
    
    return result
    

def GLI(alpha, f_name, domain_start = 0.0, domain_end = 1.0, num_points = 100):
    """ Computes the 'improved' GL fractional derivative of a function for an 
        entire array of function values. The 'improved' definition uses the 
        3-point Lagrange interpolation found in:
            
            Oldham, K. & Spanier, J. (1974). The Fractional Calculus: Theory
                and Applications of Differentiation and Integration to Arbitrary 
                Order. Academic Press, Inc.
        
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
        >>> DF_poly = GLI(-0.5, lambda x: x**2 - 1)
        >>> DF_sqrt = GLI(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> DF_sqrt = GLI(0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Get interpolating values.
    IN = GLIinterpolat(0.5)
    I = [IN.prv,IN.crr,IN.nxt]
    
    # Get array of generalized binomial coefficients.
    b_coeffs = GLcoeffs(0.5,num_points)
    
    # Calculate the improved GL differintegral using convolution.
    GLI = np.zeros(num_points)
    for i in range(3,num_points):
        F = f_name[:i]
        L = len(F)
        B = b_coeffs[:(L-2)]
        G = np.convolve(F,B,'valid')
        GLI[i] = sum(G*I)
        
    return GLI*step_size**-alpha

def RLcoeffs(index_k, index_j, alpha):
    """Calculates coefficients for the RL differintegral operator.
    
    see Baleanu, D., Diethelm, K., Scalas, E., and Trujillo, J.J. (2012). Fractional
        Calculus: Models and Numerical Methods. World Scientific.
    """
    
    if index_j == 0:
        return ((index_k-1)**(1-alpha)-(index_k+alpha-1)*index_k**-alpha)
    elif index_j == index_k:
        return 1
    else:
        return ((index_k-index_j+1)**(1-alpha)+(index_k-index_j-1)**(1-alpha)-2*(index_k-index_j)**(1-alpha))
    
def RLpoint(alpha, f_name, domain_start = None, domain_end = None, num_points = None):
    """Calculate the RL differintegral at a point with the trapezoid rule.
    
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
        >>> RL_sqrt = RLpoint(0.5, lambda x: np.sqrt(x))
        >>> RL_sqrt = RLpoint(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> RL_sqrt = RLpoint(-0.5, f, 0., 1., 100)    # Where f is defined elsewhere.
    """
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    C = 1/math.gamma(2-alpha)
    
    RL = 0
    for index_j in range(num_points+1):
        coeff = RLcoeffs(num_points, index_j, alpha)
        RL += coeff*f_values[index_j]
        
    RL *= C*step_size**-alpha
    return RL

def RLmatrix(alpha, N):
    """ Define the coefficient matrix for the RL algorithm. """
    
    coeffMatrix = np.zeros((N,N))
    for i in range(N):
        for j in range(i):
            coeffMatrix[i,j] = RLcoeffs(i,j,alpha)
    
    # Place 1 on the main diagonal.
    np.fill_diagonal(coeffMatrix,1)
    return coeffMatrix/math.gamma(2-alpha)

def RL(alpha, f_name, domain_start = 0.0, domain_end = 1.0, num_points = 100):
    """ Calculate the RL algorithm using a trapezoid rule over 
        an array of function values.
        
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
            
    Output
    ======
        RL : float 1d-array
            Each element of the array is the RL differintegral evaluated at the 
            corresponding function array index.
    
    Examples:
        >>> RL_sqrt = RL(0.5, lambda x: np.sqrt(x))
        >>> RL_sqrt = RL(0.5, lambda x: np.sqrt(x), 0., 1., 100)
        >>> RL_poly = RL(-0.5, lambda x: x**2 - 1, 0., 1., 100)
    """
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Calculate the RL differintegral.
    D = RLmatrix(alpha, num_points)
    RL = step_size**-alpha*np.dot(D, f_values)
    return RL