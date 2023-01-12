from __future__ import print_function

import numpy as np

def isInteger(n):
    if n.imag:
        return False
    if float(n.real).is_integer():
        return True
    else:
        return False

def isPositiveInteger(n):
    return isInteger(n) and n > 0

def checkValues(alpha, domain_start, domain_end, num_points, support_complex_alpha=False):
    """ Type checking for valid inputs. """
    
    assert isPositiveInteger(num_points), "num_points is not an integer: %r" % num_points
    
    assert isinstance(domain_start, (int, np.integer, float, np.floating)),\
                     "domain_start must be integer or float: %r" % domain_start
        
    assert isinstance(domain_end, (int, np.integer, float, np.floating)),\
                     "domain_end must be integer or float: %r" % domain_end

    if not support_complex_alpha:
        assert not isinstance(alpha, complex), "Complex alpha not supported for this algorithm."
    
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
    """ Returns the Pochhammer symbol (a)_n. a can be any complex or real number 
        except the negative integers and 0. n can be any nonnegative real.
    """
    if isPositiveInteger(n):
        # Compute the Pochhammer symbol.
        n = int(n)
        if n == 0:
            return 1.0
        else:
            poch = 1
            for j in range(n):
                poch *= (a + j)
            return poch

    # if a and a + n are both nonpositive integers, we can use another formula...
    # see here https://www.mathworks.com/help/symbolic/sym.pochhammer.html
    if isPositiveInteger(-1 * a) and isPositiveInteger(-1 * a - n):
        sign = -1 if np.abs(n % 2) == 1 else 1
        return sign * Gamma(1 - a) / Gamma(1 - a - n)
    return Gamma(a + n) / Gamma(a)
    
def Gamma(z):
    """ Paul Godfrey's Gamma function implementation valid for z complex.
        This is converted from Godfrey's Gamma.m Matlab file available at
        https://www.mathworks.com/matlabcentral/fileexchange/3572-gamma.
        15 significant digits of accuracy for real z and 13
        significant digits for other values.
    """
    if not (type(z) == type(1+1j)):
        if isPositiveInteger(-1 * z):
            return np.inf
        from math import gamma
        return gamma(z)

    siz = np.size(z)
    zz = z
    f = np.zeros(2,)
        
    # Find negative real parts of z and make them positive.
    if type(z) == 'complex':
        Z = [z.real,z.imag]
        if Z[0] < 0:
            Z[0]  = -Z[0]
            z = np.asarray(Z)
            z = z.astype(complex)
    
    g = 607/128.
    
    c = [0.99999999999999709182,\
          57.156235665862923517,\
         -59.597960355475491248,\
          14.136097974741747174,\
        -0.49191381609762019978,\
        .33994649984811888699e-4,\
        .46523628927048575665e-4,\
       -.98374475304879564677e-4,\
        .15808870322491248884e-3,\
       -.21026444172410488319e-3,\
        .21743961811521264320e-3,\
       -.16431810653676389022e-3,\
        .84418223983852743293e-4,\
       -.26190838401581408670e-4,\
        .36899182659531622704e-5]
    
    if z == 0 or z == 1:
        return 1.
    
    # Adjust for negative poles.
    if (np.round(zz) == zz) and (zz.imag == 0) and (zz.real <= 0):
        return np.inf
        
    z = z - 1
    zh = z + 0.5
    zgh = zh + g
    
    # Trick for avoiding floating-point overflow above z = 141.
    zp = zgh**(zh*0.5)
    ss = 0.
    
    for pp in range(len(c)-1,0,-1):
        ss += c[pp]/(z+pp)
        
    sq2pi =  2.5066282746310005024157652848110;
    f = (sq2pi*(c[0]+ss))*((zp*np.exp(-zgh))*zp)
    
    # Adjust for negative real parts.
    #if zz.real < 0:
    #    F = [f.real,f.imag]
    #    F[0] = -np.pi/(zz.real*F[0]*np.sin(np.pi*zz.real))
    #    f = np.asarray(F)
    #    f = f.astype(complex)
    
    if type(zz) == 'complex':
        return f.astype(complex)
    elif isPositiveInteger(zz):
        f = np.round(f)
        return f.astype(int)
    else:
        return f

def Beta(x,y):
    """ Beta function using Godfrey's Gamma function. """
    
    return Gamma(x)*Gamma(y)/Gamma(x+y)
    
def MittagLeffler(a, b, x, num_terms=50, *, ignore_special_cases=False):
    ''' Calculate the Mittag-Leffler function by checking for special cases, and trying to 
        reduce the parameters. If neither of those work, it just brute forces it.
        
        Parameters
       ==========
        a : float
            The first parameter of the Mittag-Leffler function.
        b : float
            The second parameter of the Mittag-Leffler function
        x : float or 1D-array of floats
            The value or values to be evaluated at.
        num_terms : int
            The number of terms to calculate in the sum. Ignored if 
            a special case can be used instead. Default value is 100.
        ignore_special_cases : bool
            Don't use the special cases, use the series definition.
            Probably only useful for testing. Default value is False.
    '''
    # check for quick special cases
    if not ignore_special_cases:
        if a == 0:
            if (np.abs(x) < 1).all():
                return 1 / Gamma(b) * 1 / (1 - x)
            return x * np.inf
        elif a == 0.5 and b == 1:
            # requires calculation of the complementary error function
            pass
        elif a == 1 and b == 1:
            return np.exp(x)
        elif a == 2 and b == 1:
            return np.cosh(np.sqrt(x))
        elif a == 1 and b == 2:
            return (np.exp(x) - 1) / x
        elif a == 2 and b == 2:
            return np.sinh(np.sqrt(x)) / np.sqrt(x)
    # manually calculate with series definition
    exponents = np.arange(num_terms)
    exp_vals = np.array([x]).T ** exponents
    gamma_vals = np.array([Gamma(exponent * a + b) for exponent in exponents])
    return np.sum(exp_vals / gamma_vals, axis=1)


def GLcoeffs(alpha,n):
    """ Computes the GL coefficient array of size n. 
    
        These coefficients can be used for both the GL 
        and the improved GL algorithm.
    """ 
    
    # Validate input.
    isPositiveInteger(n)
    
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
        >>> DF_poly = GLpoint(-0.5, lambda x: 3*x**2 - 9*x + 2)
        >>> DF_sqrt = GLpoint(0.5, lambda x: np.sqrt(x), 0., 1., 100)
    """
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
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
    """
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
       
    # Get the convolution filter.
    b_coeffs = GLcoeffs(alpha, num_points-1)
    
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
        >>> GLI_poly = GLI(-0.5, lambda x: x**2 - 1)
        >>> GLI_sqrt = GLI(0.5, lambda x: np.sqrt(x), 0., 1., 100)
    """
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    # Get interpolating values.
    IN = GLIinterpolat(0.5)
    I = [IN.prv,IN.crr,IN.nxt]
    
    # Get array of generalized binomial coefficients.
    b_coeffs = GLcoeffs(0.5,num_points)
    
    # Calculate the improved GL differintegral using convolution.
    GLI = np.zeros(num_points)
    for i in range(3,num_points):
        F = f_values[:i]
        L = len(F)
        B = b_coeffs[:(L-2)]
        G = np.convolve(F,B,'valid')
        GLI[i] = sum(G*I)
        
    return GLI*step_size**-alpha

def CRONE(alpha, f_name):
    """Calculates the GL derivative approximation using the CRONE operator.
    
    
    
        see Mathieu, B., Melchior, P., Oustaloup, A., and Ceyral, Ch. (2003). Fractional
            differentiation for edge detection. Signal Processing, 83, pp. 2421 -- 2432.
    
    """
    class Error(Exception):
        pass

    class InputError(Error):
        def __init__(self, expr, msg):
            self.expr = expr
            self.msg = msg

    def CRONEfilter(siz, alpha):
        """Creates CRONE convolution filter."""
                
        if (siz % 2) != 0:
            w = siz
            stop = int((siz-1)/2)
            print(stop)
        else:
            w = siz + 1
            stop = int(siz/2)
            
        D = GLcoeffs(alpha, stop)
        D1 = D
        D = np.flip(D, axis = 0)
        
        np.append(D,0)
        np.append(D,-D1)
        
        return D
    
    if len(np.shape(f_name)) > 1:
        [rows,cols] = np.shape(f_name)
        imgx = np.zeros((rows,cols))
        imgy = np.zeros((rows,cols))
                
        # Define the CRONE operators with the correct sizes.
        CRONEx = CRONEfilter(cols, alpha) # cols is the width of the matrix
        CRONEy = CRONEfilter(rows, alpha) # rows is the height of the matrix
        
        for i in range(rows):
            imgx[i,:] = np.convolve(f_name[i,:], CRONEx, mode = 'same')
            
        for j in range(cols):
            imgy[:,j] = np.convolve(f_name[:,j], CRONEy, mode = 'same')
            
        return imgx, imgy
    
    elif len(np.shape(f_name)) == 1:
        w = len(f_name)
        CRONEx = CRONEfilter(w, alpha) # w is the length of the array
        
        imgx = np.convolve(f_name, CRONEx, mode = 'same')
        
        return imgx
    
    else:
        raise InputError(f_name, 'f_name must have dimension <= 2')

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
    
def RLpoint(alpha, f_name, domain_start = 0.0, domain_end = 1.0, num_points = 100):
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
        >>> RL_poly = RLpoint(0.5, lambda x: x**2 - 4*x - 1, 0., 1., 100)
    """
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    C = 1/Gamma(2-alpha)
    
    RL = 0
    for index_j in range(num_points):
        coeff = RLcoeffs(num_points-1, index_j, alpha)
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
    return coeffMatrix/Gamma(2-alpha)

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
        >>> RL_poly = RL(0.5, lambda x: x**2 - 1, 0., 1., 100)
    """
    
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)
    
    # Calculate the RL differintegral.
    D = RLmatrix(alpha, num_points)
    RL = step_size**-alpha*np.dot(D, f_values)
    return RL

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

def CaputoL1point(alpha, f_name, domain_start=0, domain_end=1, num_points=100):
    ''' Calculate the Caputo derivative of a function at a point using the L1 method.

    see Karniadakis, G.E.. (2019). Handbook of Fractional Calculus with Applications
    Volume 3: Numerical Methods. De Gruyter.

    Parameters
    ==========
        alpha : float
            The order of the differintegral to be computed. Must be in (0, 1)
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
        L1 : float
            The Caputo L1 integral evaluated at the corresponding point.
    '''

    if alpha <= 0 or alpha >= 1:
        raise ValueError('Alpha must be in (0, 1) for this method.')

    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    f_values = np.array(f_values)
    j_values = np.arange(0, num_points-1)
    coefficients = (j_values + 1) ** (1 - alpha) - (j_values) ** (1 - alpha)
    f_differences = f_values[1:] - f_values[:-1]
    f_differences = f_differences[::-1]
    L1 = 1 / Gamma(2 - alpha) * np.sum(np.multiply(coefficients * step_size**(-alpha), f_differences))
    
    return L1

def CaputoL2point(alpha, f_name, domain_start=0, domain_end=1, num_points=100):
    ''' Calculate the Caputo derivative of a function at a point using the L2 method.
        A note: this method requires evaluation of the point f(domain_end + step size),
        and currently will only work if `f_name` is a callable function.

    see Karniadakis, G.E.. (2019). Handbook of Fractional Calculus with Applications
    Volume 3: Numerical Methods. De Gruyter.

    Parameters
    ==========
        alpha : float
            The order of the differintegral to be computed. Must be in (1, 2).
        f_name : function handle or lambda function
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
        L2 : float
            The Caputo L2 integral evaluated at the corresponding point.
    '''
    if alpha <= 1 or alpha >= 2:
        raise ValueError('Alpha must be in (1, 2) for this method.')
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    def b_coes(alpha, j):
        return (j + 1) ** (2 - alpha) - j ** (2 - alpha)

    # start with the point outside of the domain
    L2 = b_coes(alpha, 0) * (f_values[num_points - 2] + f_name(num_points * step_size) - 2 * f_values[num_points - 1]) #f_name(num_points * step_size)
    for k in range(1, num_points - 2):
        L2 += b_coes(alpha, k) * (f_values[num_points - 2 - k] + f_values[num_points - k] - 2 * f_values[num_points - k - 1])
    return L2 * step_size ** (-1 * alpha) / Gamma(3 - alpha)


def CaputoL2Cpoint(alpha, f_name, domain_start=0, domain_end=1, num_points=100):
    ''' Calculate the Caputo derivative of a function at a point using the L2C method.
        A note: this method requires evaluation of the points f(domain_end + step size)
        and f(-step_size), and currently will only work if `f_name` is a callable 
        function.

    see Karniadakis, G.E.. (2019). Handbook of Fractional Calculus with Applications
    Volume 3: Numerical Methods. De Gruyter.

    Parameters
    ==========
        alpha : float
            The order of the differintegral to be computed. Must be in (0, 2).
        f_name : function handle or lambda function
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
        L2C : float
            The Caputo L2C integral evaluated at the corresponding point.
    '''
    if alpha <= 0 or alpha >= 2:
        raise ValueError('Alpha must be in (0, 1) or (1, 2) for this method.')

    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    def b_coes(alpha, j):
        return (j + 1) ** (2 - alpha) - j ** (2 - alpha)

    # start with the points outside of the domain
    L2C = b_coes(alpha, 0) * (f_values[num_points - 3] - f_values[num_points - 2] - f_values[num_points - 1] + f_name(num_points * step_size)) #f_name(num_points * step_size)
    L2C += b_coes(alpha, num_points - 2) * (f_name(-1 * step_size) + f_values[2] - f_values[1] - f_values[0])
    for k in range(1, num_points - 2):
        L2C += b_coes(alpha, k) * (f_values[num_points - 3 - k] - f_values[num_points - k - 2] - f_values[num_points - k - 1] + f_values[num_points - k]) 
    L2C *= step_size ** (-1 * alpha) / Gamma(3 - alpha) * 0.5

    return L2C

def CaputoFromRLpoint(alpha, f_name, domain_start=0, domain_end=1, num_points=100):
    ''' Calculate the Caputo derivative of a function at a point using the conversion
        formula from the RL differintegrals. DOESN'T CURRENTLY WORK.

    see Du, R., Yan, Y. and Liang, Z., (2019). A high-order scheme to
        approximate the caputo fractional derivative and its application
        to solve the fractional diffusion wave equation, Journal of
        Computational Physics, 376, pp. 1312-1330

    Parameters
    ==========
        alpha : float
            The order of the differintegral to be computed. Must be in (1, 2).
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
        C : float
            The Caputo integral evaluated at the corresponding point.
    '''
    if alpha <= 1 or alpha >= 2:
        raise ValueError('Alpha must be in (1, 2) for this method.')

    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start
    
    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    C = 0
    C -= f_values[0] * domain_end ** (-1 * alpha) / Gamma(1 - alpha)
    C -= (f_values[1] - f_values[0]) / step_size * domain_end ** (1 - alpha) / Gamma(2 - alpha)
    C += RLpoint(alpha - 2, f_name, domain_start, float(domain_end + step_size), num_points) / step_size ** 2
    C -= 2 * RLpoint(alpha - 2, f_name, domain_start, float(domain_end), num_points) / step_size ** 2
    C -= RLpoint(alpha - 2, f_name, domain_start, float(domain_end - step_size), num_points) / step_size ** 2
    return C

def PCcoeffs(alpha, j, n):
    if 1 < alpha:
        if j == 0:
            return (n+1)**alpha * (alpha - n) + n**alpha * (2 * n - alpha - 1) - (n - 1)**(alpha + 1)
        elif j == n:
            return 2 ** (alpha + 1) - alpha - 3
        return (n - j + 2) ** (alpha + 1) + 3 * (n - j) ** (alpha + 1) - 3 * (n - j + 1) ** (alpha + 1) - (n - j - 1) ** (alpha + 1)

def PCsolver(initial_values, alpha, f_name, domain_start=0, domain_end=1, num_points=100):
    """ Solve an equation of the form D[y(x)]=f(x, y(x)) using the predictor-corrector
        method, modified to be compatible with fractional derivatives.

    see Deng, W. (2007) Short memory principle and a predictor–corrector approach for 
        fractional differential equations. Journal of Computational and Applied 
        Mathematics.

    test examples from
        Baskonus, H.M., Bulut, H. (2015) On the numerical solutions of some fractional
        ordinary differential equations by fractional Adams-Bashforth-Moulton method.
        De Gruyter.
        Weilbeer, M. (2005) Efficient Numerical Methods for Fractional Differential
        Equations and their Analytical Background. 
        
    Parameters
    ==========
        initial_values : float 1d-array
            A list of initial values for the IVP. There should be as many IVs
            as ceil(alpha).
        alpha : float
            The order of the differintegral in the equation to be computed.
        f_name : function handle or lambda function
            This is the function on the right side of the equation, and should
            accept two variables; first the independant variable, and second
            the equation to be solved.
        domain_start : float
            The left-endpoint of the function domain. Default value is 0.
        domain_end : float
            The right-endpoint of the function domain; the point at which the 
            differintegral is being evaluated. Default value is 1.
        num_points : integer
            The number of points in the domain. Default value is 100.
            
    Output
    ======
        y_correction : float 1d-array
            The calculated solution to the IVP at each of the points 
            between the left and right endpoint.

    Examples:
        >>> f_name = lambda x, y : y - x - 1
        >>> initial_values = [1, 1]
        >>> y_solved = PCsolver(initial_values, 1.5, f_name)
        >>> theoretical = np.linspace(0, 1, 100) + 1
        >>> np.allclose(y_solved, theoretical)
        True
    """
    x_points = np.linspace(domain_start, domain_end, num_points)
    step_size = x_points[1] - x_points[0]
    y_correction = np.zeros(num_points, dtype='complex_')
    y_prediction = np.zeros(num_points, dtype='complex_')
    
    y_prediction[0] = initial_values[0]            
    y_correction[0] = initial_values[0]
    for x_index in range(num_points - 1):
        initial_value_contribution = 0
        if 1 < alpha and alpha < 2:
            initial_value_contribution = initial_values[1] * step_size
        elif 2 < alpha:
            for k in range(1, int(np.ceil(alpha))):
                initial_value_contribution += initial_values[k] / np.math.factorial(k) * (x_points[x_index + 1] ** k - x_points[x_index] ** k) 
        elif alpha < 1:
            raise ValueError('Not yet supported!')
        y_prediction[x_index + 1] += initial_value_contribution
        y_prediction[x_index + 1] += y_correction[x_index]
        y_prediction[x_index + 1] += step_size ** alpha / Gamma(alpha + 1) * f_name(x_points[x_index], y_correction[x_index])
        subsum = 0
        for j in range(x_index + 1):
            subsum += PCcoeffs(alpha, j, x_index) * f_name(x_points[j], y_correction[j])
        y_prediction[x_index + 1] += step_size ** alpha / Gamma(alpha + 2) * subsum

        y_correction[x_index + 1] += initial_value_contribution
        y_correction[x_index + 1] += y_correction[x_index]
        y_correction[x_index + 1] += step_size ** alpha / Gamma(alpha + 2) * alpha * f_name(x_points[x_index], y_correction[x_index])
        y_correction[x_index + 1] += step_size ** alpha / Gamma(alpha + 2) * f_name(x_points[x_index + 1], y_prediction[x_index + 1])
        y_correction[x_index + 1] += step_size ** alpha / Gamma(alpha + 2) * subsum
    
    return y_correction
