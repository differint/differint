from typing import Callable

import numpy as np

from .utils import functionCheck, checkValues
from .core import GL


def GLpoint_via_GL(
    alpha: float,
    f_name: Callable[[np.ndarray], np.ndarray] | list[float] | np.ndarray,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """
    Efficiently computes the Grünwald-Letnikov fractional derivative at the endpoint
    by evaluating the full array via the optimized GL function and returning the last value.

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
    values = GL(alpha, f_name, domain_start, domain_end, num_points)
    return float(values[-1])


def GLpoint_direct(
    alpha: float,
    f_name: Callable | np.ndarray | list,
    domain_start: float = 0.0,
    domain_end: float = 1.0,
    num_points: int = 100,
) -> float:
    """Efficiently computes the Grünwald-Letnikov fractional derivative
    of a function at the right endpoint of the domain using a direct single-pass
    recurrence relation.

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
    f_values, _ = functionCheck(f_name, domain_start, domain_end, num_points)

    # Calculate the GL differintegral, avoiding the explicit calculation of
    # the gamma function.
    GL_previous = f_values[1]
    for index in range(2, num_points):
        GL_current = (
            GL_previous * (num_points - alpha - index - 1) / (num_points - index)
            + f_values[index]
        )
        GL_previous = GL_current

    return GL_current * (num_points / (domain_end - domain_start)) ** alpha  # type: ignore
