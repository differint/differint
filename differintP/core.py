from __future__ import print_function

from typing import Callable

import numpy as np

from scipy.special import gamma as Gamma
from numba import njit

#  CuPy dependency for GPU-accelerated GL_gpu
from .gpu_utils import cupy_manager

from .utils import checkValues, functionCheck





#########################################################################################
########################################## GL ###########################################
#########################################################################################


def GLcoeffs(alpha: float, n: int) -> np.ndarray:
    """Vectorized GL coefficient computation"""
    """ Computes the GL coefficient array of size n.

        These coefficients can be used for both the GL
        and the improved GL algorithm.
    """
    if n == 0:
        return np.array([1.0])

    # Preallocate factors array
    factors = np.ones(n + 1)

    # Compute the multiplicative factors for positions 1 to n
    numerators = -alpha + np.arange(n)
    denominators = np.arange(1, n + 1)
    factors[1:] = numerators / denominators

    # Compute cumulative product
    return np.cumprod(factors)


def GL(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> np.ndarray:
    """Optimized GL fractional derivative using precomputation"""
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
    # Domain handling
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Generate points and get function values
    x = np.linspace(domain_start, domain_end, num_points)
    step_size = x[1] - x[0]

    if callable(f_name):
        f_values = f_name(x)
    else:
        f_values = np.asarray(f_name)
        if len(f_values) != num_points:
            raise ValueError("Function array length doesn't match num_points")

    # Precompute coefficients (vectorized)
    b_coeffs = GLcoeffs(alpha, num_points - 1)

    # FFT convolution
    B = np.fft.rfft(b_coeffs, n=num_points)
    F = np.fft.rfft(f_values)
    result = np.fft.irfft(F * B, n=num_points)[:num_points] * step_size**-alpha

    return result


@njit
def _GLpoint_loop(alpha: float, f_values: np.ndarray, step: float) -> float:
    k = f_values.shape[0] - 1
    acc = 0.0
    c_val = 1.0
    for j in range(k + 1):
        acc += c_val * f_values[k - j]
        if j < k:
            c_val *= (-alpha + j) / (j + 1)
    return step ** (-alpha) * acc


def GLpoint(
    alpha: float,
    f_name: Callable[[np.ndarray], np.ndarray] | list[float] | np.ndarray,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """
    Efficient, robust single-point Grünwald-Letnikov fractional derivative
    using a direct recurrence (C++-style) in a JIT-compiled kernel.

    Parameters
    ----------
    alpha : float
        The order of the fractional derivative.
    f_name : Callable[[np.ndarray], np.ndarray] or Sequence[float] or np.ndarray
        The function to differentiate (callable or array-like).
    domain_start : float, optional
        The starting value of the domain (default is 0.0).
    domain_end : float, optional
        The ending value of the domain (default is 1.0).
    num_points : int, optional
        Number of discretization points (default is 100).

    Returns
    -------
    float
        The Grünwald-Letnikov fractional derivative at the endpoint.
    """
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    x = np.linspace(domain_start, domain_end, num_points)
    if callable(f_name):
        f_values = f_name(x)
    else:
        f_values = np.asarray(f_name)
        if len(f_values) != num_points:
            raise ValueError("Function array length doesn't match num_points")

    step = (domain_end - domain_start) / (num_points - 1)

    return _GLpoint_loop(alpha, f_values, step)


#########################################################################################
######################################## GL Gpu #########################################
#########################################################################################


def _gpu_GLcoeffs(
    alpha: float,
    n: int,
) -> np.ndarray:
    """GPU version of GL coefficients"""
    if cupy_manager.HAS_CUPY:
        k = cupy_manager.cp.arange(n + 1)  # type: ignore
        factors = cupy_manager.cp.ones(n + 1)  # type: ignore
        if n > 0:
            numerators = -alpha + cupy_manager.cp.arange(n)  # type: ignore
            denominators = cupy_manager.cp.arange(1, n + 1)  # type: ignore
            factors[1:] = numerators / denominators
        return cupy_manager.cp.cumprod(factors)  # type: ignore
    else:
        raise RuntimeError("CuPy not available")


def GL_gpu(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> np.ndarray:
    """GPU-accelerated GL fractional derivative"""

    """
    GPU-accelerated GL fractional derivative (requires CuPy)

    Note: This function requires CuPy with CUDA support. Install with:
        pip install cupy-cuda11x   or   cupy-cuda12x
        the code is only tested with cupy-cuda12x
    """

    if not cupy_manager.HAS_CUPY:
        raise RuntimeError(
            "CuPy not available. Please install CuPy with CUDA support to use this function.\n"
            "Install via: pip install cupy-cuda11x (replace 11x with your CUDA version)"
        )

    # Convert to GPU arrays
    x = cupy_manager.cp.linspace(domain_start, domain_end, num_points)  # type: ignore
    if callable(f_name):
        f_values = cupy_manager.cp.asarray(f_name(cupy_manager.cp.asnumpy(x)))  # type: ignore
    else:
        f_values = cupy_manager.cp.asarray(f_name)  # type: ignore

    # GPU-accelerated computation
    b_coeffs = _gpu_GLcoeffs(alpha, num_points - 1)
    B = cupy_manager.cp.fft.rfft(b_coeffs)  # type: ignore
    F = cupy_manager.cp.fft.rfft(f_values)  # type: ignore
    result = cupy_manager.cp.fft.irfft(F * B, n=num_points)[:num_points] * (  # type: ignore
        (x[1] - x[0]) ** -alpha
    )

    return cupy_manager.cp.asnumpy(result)  # Convert back to CPU # type: ignore


#########################################################################################
##################################### GLI - Crone #######################################
#########################################################################################

@njit
def _GLI_core(alpha: float, f_values: np.ndarray, b_coeffs: np.ndarray) -> np.ndarray:
    """
    Developer function: Efficiently computes the improved Grünwald-Letnikov (GLI) fractional derivative core.

    Given function values and GL coefficients (as 1D NumPy arrays), this routine:
    - Applies a 3-point quadratic (Lagrange) interpolation in a moving window.
    - For each point, performs a local convolution with a flipped GL coefficient array.
    - Combines results with fixed interpolation weights for high-order accuracy.

    Intended for internal use with Numba JIT compilation for speed.
    Inputs must be 1D NumPy arrays of matching length.
    """
    num_points = f_values.shape[0]
    GLI = np.zeros(num_points)
    prv = alpha * (alpha - 2) / 8
    crr = (4 - alpha * alpha) / 4
    nxt = alpha * (2 + alpha) / 8
    I = np.array([prv, crr, nxt])
    for i in range(3, num_points):
        L = i
        F = f_values[:L]
        B = b_coeffs[: (L - 2)]
        G = np.zeros(3)
        B_flip = B[::-1]
        for m in range(3):
            s = 0.0
            for n in range(len(B)):
                s += F[m + n] * B_flip[n]
            G[m] = s
        GLI[i] = np.sum(G * I)
    return GLI


def GLI(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> np.ndarray:
    """
    Computes the 'improved' Grünwald-Letnikov (GL) fractional derivative of a function over an entire domain,
    using the quadratic 3-point Lagrange interpolation method described by Oldham & Spanier (1974).

    This implementation applies a three-point moving-window interpolation to the function values,
    followed by a specialized convolution with Grünwald-Letnikov coefficients at each evaluation point.
    The result is a higher-order, more accurate fractional derivative estimate compared to the standard GL method,
    especially for smooth functions.

    For performance, a Numba-accelerated helper is used internally.

    Parameters
    ----------
    alpha : float
        Order of the fractional differintegral to compute.
    f_name : Callable, list, or 1d-array
        Function, lambda, or sequence of function values to be differintegrated.
        If callable, it will be evaluated at `num_points` evenly spaced points.
    domain_start : float, optional
        Start of the domain. Default is 0.0.
    domain_end : float, optional
        End of the domain. Default is 1.0.
    num_points : int, optional
        Number of points to evaluate in the domain. Default is 100.

    Returns
    -------
    np.ndarray
        Array of 'improved' GL fractional derivative values at each grid point.

    Notes
    -----
    This method implements the "improved" Grünwald-Letnikov definition by first applying
    quadratic interpolation to each interior function value:
        interpolated[i] = a * f[i-1] + b * f[i] + c * f[i+1]
    where the coefficients (a, b, c) depend on `alpha` and are chosen to reduce discretization error.

    The resulting interpolated array is then convolved with the GL binomial coefficients,
    with the result scaled by the step size raised to `-alpha`. For performance, the function
    chooses between direct and FFT convolution based on array size.

    References
    ----------
    Oldham, K. & Spanier, J. (1974). The Fractional Calculus: Theory and Applications of Differentiation and Integration to Arbitrary Order. Academic Press.

    Examples
    --------
    >>> GLI_poly = GLI(-0.5, lambda x: x**2 - 1)
    >>> GLI_sqrt = GLI(0.5, lambda x: np.sqrt(x), 0., 1., 100)
    """
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Evaluate function on grid
    x = np.linspace(domain_start, domain_end, num_points)
    step_size = (domain_end - domain_start) / (num_points - 1)
    if callable(f_name):
        f_values = np.asarray(f_name(x))
    else:
        f_values = np.asarray(f_name)
        if len(f_values) != num_points:
            raise ValueError("Function array length doesn't match num_points")

    b_coeffs = GLcoeffs(alpha, num_points)
    GLI_vals = _GLI_core(alpha, f_values, b_coeffs)
    return GLI_vals * step_size**-alpha



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
            stop = int((siz - 1) / 2)
            print(stop)
        else:
            w = siz + 1
            stop = int(siz / 2)

        D = GLcoeffs(alpha, stop)
        D1 = D
        D = np.flip(D, axis=0)

        np.append(D, 0)
        np.append(D, -D1)

        return D

    if len(np.shape(f_name)) > 1:
        [rows, cols] = np.shape(f_name)
        imgx = np.zeros((rows, cols))
        imgy = np.zeros((rows, cols))

        # Define the CRONE operators with the correct sizes.
        CRONEx = CRONEfilter(cols, alpha)  # cols is the width of the matrix
        CRONEy = CRONEfilter(rows, alpha)  # rows is the height of the matrix

        for i in range(rows):
            imgx[i, :] = np.convolve(f_name[i, :], CRONEx, mode="same")

        for j in range(cols):
            imgy[:, j] = np.convolve(f_name[:, j], CRONEy, mode="same")

        return imgx, imgy

    elif len(np.shape(f_name)) == 1:
        w = len(f_name)
        CRONEx = CRONEfilter(w, alpha)  # w is the length of the array

        imgx = np.convolve(f_name, CRONEx, mode="same")

        return imgx

    else:
        raise InputError(f_name, "f_name must have dimension <= 2")


#########################################################################################
######################################### RLmatrix ######################################
#########################################################################################


def RLmatrix(alpha, N):
    """Vectorized RL coefficient matrix generation"""
    # Precompute all required powers
    k = np.arange(N)
    v = np.zeros(N + 2)  # +2 to avoid index issues
    v[1:] = np.power(np.arange(1, N + 2), 1 - alpha)

    # Initialize coefficient matrix
    coeffMatrix = np.zeros((N, N))

    # Set diagonal to 1
    np.fill_diagonal(coeffMatrix, 1)

    # First column (j=0)
    i_vals = np.arange(1, N)
    coeffMatrix[i_vals, 0] = v[i_vals - 1] - (i_vals + alpha - 1) * np.power(
        i_vals, -alpha
    )

    # Main coefficients using vectorization
    for k_val in range(1, N - 1):
        rows = np.arange(k_val + 1, N)
        cols = rows - k_val
        coeffMatrix[rows, cols] = v[k_val + 1] + v[k_val - 1] - 2 * v[k_val]

    # Normalize with Gamma function
    return coeffMatrix / Gamma(2 - alpha)

#########################################################################################
####################################### RLcoeffs ########################################
#########################################################################################

def RLcoeffs(index_k, index_j, alpha):
    """Calculates coefficients for the RL differintegral operator.

    see Baleanu, D., Diethelm, K., Scalas, E., and Trujillo, J.J. (2012). Fractional
        Calculus: Models and Numerical Methods. World Scientific.
    """
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
    if index_j == 0:
        return (index_k - 1) ** (1 - alpha) - (index_k + alpha - 1) * index_k**-alpha
    elif index_j == index_k:
        return 1
    else:
        return (
            (index_k - index_j + 1) ** (1 - alpha)
            + (index_k - index_j - 1) ** (1 - alpha)
            - 2 * (index_k - index_j) ** (1 - alpha)
        )

#########################################################################################
######################################## RLpoint ########################################
#########################################################################################

def RLpoint(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """Optimized RL fractional derivative calculation 8x - 60x speed"""
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

    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Generate evaluation points
    x = np.linspace(domain_start, domain_end, num_points)
    step_size = x[1] - x[0]

    # Get function values (optimized)
    if callable(f_name):
        f_values = f_name(x)
    else:
        f_values = np.asarray(f_name)
        if len(f_values) != num_points:
            raise ValueError("Function array length doesn't match num_points")

    # Precompute all coefficients in vectorized form
    k = num_points - 1  # Fixed evaluation index (endpoint)
    j = np.arange(num_points)

    # Initialize coefficient array
    coeffs = np.zeros(num_points)

    # Case 1: j == 0
    mask_j0 = j == 0
    if k > 0:  # Only compute if k > 0
        coeffs[mask_j0] = (k - 1) ** (1 - alpha) - (k + alpha - 1) * k**-alpha

    # Case 2: j == k
    mask_jk = j == k
    coeffs[mask_jk] = 1

    # Case 3: All other indices
    mask_other = ~mask_j0 & ~mask_jk
    d = k - j[mask_other]
    coeffs[mask_other] = (
        (d + 1) ** (1 - alpha) + (d - 1) ** (1 - alpha) - 2 * d ** (1 - alpha)
    )

    # Final calculation
    C = 1 / Gamma(2 - alpha)
    return C * step_size**-alpha * np.dot(coeffs, f_values)

#########################################################################################
########################################### RL ##########################################
#########################################################################################

def RL(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> np.ndarray:
    """Optimized RL fractional derivative calculation 14x speed"""
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
    # Domain validation and flipping
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Generate evaluation points
    x = np.linspace(domain_start, domain_end, num_points)
    step_size = x[1] - x[0]

    # Get function values (optimized)
    if callable(f_name):
        f_values = f_name(x)
    else:
        f_values = np.asarray(f_name)
        if len(f_values) != num_points:
            raise ValueError("Function array length doesn't match num_points")

    # Compute RL differintegral
    D = RLmatrix(alpha, num_points)
    result = step_size**-alpha * (D @ f_values)
    return result


#########################################################################################
######################################### Caputo ########################################
#########################################################################################


def CaputoL1point(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """Calculate the Caputo derivative of a function at a point using the L1 method.

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
    """

    if alpha <= 0 or alpha >= 1:
        raise ValueError("Alpha must be in (0, 1) for this method.")

    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    f_values = np.array(f_values)
    j_values = np.arange(0, num_points - 1)
    coefficients = (j_values + 1) ** (1 - alpha) - (j_values) ** (1 - alpha)
    f_differences = f_values[1:] - f_values[:-1]
    f_differences = f_differences[::-1]
    L1 = (
        1
        / Gamma(2 - alpha)
        * np.sum(np.multiply(coefficients * step_size ** (-alpha), f_differences))
    )

    return L1


def CaputoL2point(
    alpha: float,
    f_name: Callable,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """Calculate the Caputo derivative of a function at a point using the L2 method.
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
    """
    if alpha <= 1 or alpha >= 2:
        raise ValueError("Alpha must be in (1, 2) for this method.")
    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    def b_coes(alpha, j):
        return (j + 1) ** (2 - alpha) - j ** (2 - alpha)

    # start with the point outside of the domain
    L2 = b_coes(alpha, 0) * (
        f_values[num_points - 2]
        + f_name(num_points * step_size)
        - 2 * f_values[num_points - 1]
    )  # f_name(num_points * step_size)
    for k in range(1, num_points - 1):
        L2 += b_coes(alpha, k) * (
            f_values[num_points - 2 - k]
            + f_values[num_points - k]
            - 2 * f_values[num_points - k - 1]
        )
    return L2 * step_size ** (-1 * alpha) / Gamma(3 - alpha)


def CaputoL2Cpoint(
    alpha: float,
    f_name: Callable,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """Calculate the Caputo derivative of a function at a point using the L2C method.
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
    """
    if alpha <= 0 or alpha >= 2:
        raise ValueError("Alpha must be in (0, 1) or (1, 2) for this method.")

    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    def b_coes(alpha, j):
        return (j + 1) ** (2 - alpha) - j ** (2 - alpha)

    # start with the points outside of the domain
    L2C = b_coes(alpha, 0) * (
        f_values[num_points - 3]
        - f_values[num_points - 2]
        - f_values[num_points - 1]
        + f_name(num_points * step_size)
    )  # f_name(num_points * step_size)
    L2C += b_coes(alpha, num_points - 2) * (
        f_name(-1 * step_size) + f_values[2] - f_values[1] - f_values[0]
    )
    for k in range(1, num_points - 2):
        L2C += b_coes(alpha, k) * (
            f_values[num_points - 3 - k]
            - f_values[num_points - k - 2]
            - f_values[num_points - k - 1]
            + f_values[num_points - k]
        )
    L2C *= step_size ** (-1 * alpha) / Gamma(3 - alpha) * 0.5

    return L2C


def CaputoFromRLpoint(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """Calculate the Caputo derivative of a function at a point using the conversion
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
    """
    if alpha <= 1 or alpha >= 2:
        raise ValueError("Alpha must be in (1, 2) for this method.")

    # Flip the domain limits if they are in the wrong order.
    if domain_start > domain_end:
        domain_start, domain_end = domain_end, domain_start

    # Check inputs.
    checkValues(alpha, domain_start, domain_end, num_points)
    f_values, step_size = functionCheck(f_name, domain_start, domain_end, num_points)

    C = 0
    C -= f_values[0] * domain_end ** (-1 * alpha) / Gamma(1 - alpha)
    C -= (
        (f_values[1] - f_values[0])
        / step_size
        * domain_end ** (1 - alpha)
        / Gamma(2 - alpha)
    )
    C += (
        RLpoint(
            alpha - 2, f_name, domain_start, float(domain_end + step_size), num_points
        )
        / step_size**2
    )
    C -= (
        2
        * RLpoint(alpha - 2, f_name, domain_start, float(domain_end), num_points)
        / step_size**2
    )
    C -= (
        RLpoint(
            alpha - 2, f_name, domain_start, float(domain_end - step_size), num_points
        )
        / step_size**2
    )
    return C


