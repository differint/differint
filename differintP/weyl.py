from typing import Callable, cast, Union, List

import numpy as np

def Weyl(
    alpha: float,
    f_name: Union[np.ndarray, List[float], Callable],
    domain_start: float = 0.0,
    domain_end: float = 2 * np.pi,
    num_points: int = 100,
) -> np.ndarray:
    """
    Weyl fractional derivative (periodic, Fourier-based).

    Numerically computes the Weyl (right-sided, periodic) fractional derivative of order `alpha`
    for a function on a uniform grid, using the FFT. This method is fast and accurate for
    periodic functions on [domain_start, domain_end].

    References:
    - Samko, Kilbas, Marichev, *Fractional Integrals and Derivatives* (see Weyl derivative, Ch. 7)

    Parameters
    ----------
    alpha : float
        Order of the derivative.
    f_name : callable or array-like
        Function or array of values to differentiate.
    domain_start, domain_end : float
        Interval (should cover one period for periodic functions).
    num_points : int
        Number of grid points.

    Returns
    -------
    df : np.ndarray
        Array of Weyl fractional derivative values at grid points.
    """

    if callable(f_name):
        # f_name = cast(Callable, f_name) # type checking
        # If f_name is callable, call it and save to a list.
        x = np.linspace(domain_start, domain_end, num_points, endpoint=False)
        f_values = f_name(x)
    else:
        f_name = cast(np.ndarray | list[float], f_name)
        num_points = np.size(f_name)
        f_values = f_name

    # Compute FFT
    fhat = np.fft.fft(f_values) # type: ignore
    L = domain_end - domain_start
    k = np.fft.fftfreq(num_points, d=L / num_points) * 2 * np.pi  # Frequency in radians

    # Fractional derivative in Fourier domain
    multiplier = np.zeros_like(k, dtype=complex)
    multiplier[1:] = (1j * k[1:]) ** alpha  # k=0 stays zero

    fhat_new = fhat * multiplier
    df = np.fft.ifft(fhat_new)
    return df.real if np.all(np.isreal(f_values)) else df # type: ignore


def Riesz(
    alpha: float,
    f_name: Union[np.ndarray, List[float], Callable],
    domain_start: float = 0.0,
    domain_end: float = 2 * np.pi,
    num_points: int = 100,
) -> np.ndarray:
    """
    Riesz fractional derivative (symmetric, Fourier-based).

    Numerically computes the Riesz fractional derivative of order `alpha`
    for a function on a uniform grid using the FFT. This operator is
    symmetric (unlike Weyl) and real for real input.

    Parameters
    ----------
    alpha : float
        Order of the derivative.
    f_name : callable or array-like
        Function or array of values to differentiate.
    domain_start, domain_end : float
        Interval (should cover one period for periodic functions).
    num_points : int
        Number of grid points.

    Returns
    -------
    df : np.ndarray
        Array of Riesz fractional derivative values at grid points.
    """
    if callable(f_name):
        x = np.linspace(domain_start, domain_end, num_points, endpoint=False)
        f_values = f_name(x)
    else:
        f_values = np.asarray(f_name)
        num_points = len(f_values)

    # FFT
    fhat = np.fft.fft(f_values)
    L = domain_end - domain_start
    k = np.fft.fftfreq(num_points, d=L / num_points) * 2 * np.pi  # radians

    # Riesz symbol: -|k|^alpha (symmetric)
    multiplier = -np.abs(k) ** alpha

    # Apply in Fourier domain
    fhat_new = fhat * multiplier
    df = np.fft.ifft(fhat_new)
    return df.real if np.all(np.isreal(f_values)) else df


