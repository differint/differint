from typing import Callable, cast

import numpy as np


def isInteger(n) -> bool:
    if n.imag:
        return False
    if float(n.real).is_integer():
        return True
    else:
        return False


def isPositiveInteger(n) -> bool:
    return isInteger(n) and n > 0


def checkValues(
    alpha: float,
    domain_start: int | float,
    domain_end: int | float,
    num_points: int,
    support_complex_alpha: bool = False,
) -> bool | None:
    """Type checking for valid inputs."""

    assert isPositiveInteger(num_points), (
        "num_points is not an integer: %r" % num_points
    )

    assert isinstance(domain_start, (int, np.integer, float, np.floating)), (
        "domain_start must be integer or float: %r" % domain_start
    )

    assert isinstance(domain_end, (int, np.integer, float, np.floating)), (
        "domain_end must be integer or float: %r" % domain_end
    )

    if not support_complex_alpha:
        assert not isinstance(
            alpha, complex
        ), "Complex alpha not supported for this algorithm."

    return


def functionCheck(
    f_name: Callable | list | np.ndarray,
    domain_start: float | int,
    domain_end: float | int,
    num_points: int,
):
    """Check if function is callable and assign function values."""

    # Define the function domain and obtain function values.
    if hasattr(f_name, "__call__"):
        f_name = cast(Callable, f_name)
        # If f_name is callable, call it and save to a list.
        x = np.linspace(domain_start, domain_end, num_points)
        f_values = list(map(lambda t: f_name(t), x))
        step_size = x[1] - x[0]
    else:
        f_name = cast(np.ndarray | list[float], f_name)
        num_points = np.size(f_name)
        f_values = f_name
        step_size = (domain_end - domain_start) / (num_points - 1)
    return f_values, step_size
