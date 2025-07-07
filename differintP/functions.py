import numpy as np

from scipy.special import gamma

from .utils import isPositiveInteger



def poch(a, n):
    """Returns the Pochhammer symbol (a)_n. a can be any complex or real number
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
                poch *= a + j
            return poch

    # if a and a + n are both nonpositive integers, we can use another formula...
    # see here https://www.mathworks.com/help/symbolic/sym.pochhammer.html
    if isPositiveInteger(-1 * a) and isPositiveInteger(-1 * a - n):
        sign = -1 if np.abs(n % 2) == 1 else 1
        return sign * gamma(1 - a) / gamma(1 - a - n)
    return gamma(a + n) / gamma(a)


def Beta(
    x: int | float | np.ndarray | complex,
    y: int | float | np.ndarray | complex,
) -> int | float | np.ndarray | complex:
    """Beta function using Scipy Gamma function."""

    return gamma(x) * gamma(y) / gamma(x + y)


def MittagLeffler(
    a: float,
    b: float,
    x: np.ndarray,
    num_terms: int = 50,
    ignore_special_cases: bool = False,
) -> np.ndarray:
    """Calculate the Mittag-Leffler function by checking for special cases, and trying to
     reduce the parameters. If neither of those work, it just brute forces it.

     Parameters
    ==========
     a : float
         The first parameter of the Mittag-Leffler function.
     b : float
         The second parameter of the Mittag-Leffler function
     x : 1D-array of floats (can be len = 1)
         The value or values to be evaluated at.
     num_terms : int
         The number of terms to calculate in the sum. Ignored if
         a special case can be used instead. Default value is 100.
     ignore_special_cases : bool
         Don't use the special cases, use the series definition.
         Probably only useful for testing. Default value is False.
    """
    # check for quick special cases
    if not ignore_special_cases:
        if a == 0:
            if (np.abs(x) < 1).all():
                return 1 / gamma(b) * 1 / (1 - x)
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
    gamma_vals = np.array([gamma(exponent * a + b) for exponent in exponents])
    return np.sum(exp_vals / gamma_vals, axis=1)