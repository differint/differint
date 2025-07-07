from __future__ import print_function


import numpy as np

from scipy.special import gamma as Gamma



#########################################################################################
######################################## PCsolver #######################################
#########################################################################################


def PCcoeffs(alpha, j, n):
    if 1 < alpha:
        if j == 0:
            return (
                (n + 1) ** alpha * (alpha - n)
                + n**alpha * (2 * n - alpha - 1)
                - (n - 1) ** (alpha + 1)
            )
        elif j == n:
            return 2 ** (alpha + 1) - alpha - 3
        return (
            (n - j + 2) ** (alpha + 1)
            + 3 * (n - j) ** (alpha + 1)
            - 3 * (n - j + 1) ** (alpha + 1)
            - (n - j - 1) ** (alpha + 1)
        )


def PCsolver(
    initial_values, alpha, f_name, domain_start=0, domain_end=1, num_points=100
):
    """Solve an equation of the form D[y(x)]=f(x, y(x)) using the predictor-corrector
        method, modified to be compatible with fractional derivatives.

    see Deng, W. (2007) Short memory principle and a predictor-corrector approach for
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
    from scipy.special import factorial

    x_points = np.linspace(domain_start, domain_end, num_points)
    step_size = x_points[1] - x_points[0]
    y_correction = np.zeros(num_points, dtype="complex")
    y_prediction = np.zeros(num_points, dtype="complex")

    y_prediction[0] = initial_values[0]
    y_correction[0] = initial_values[0]
    for x_index in range(num_points - 1):
        initial_value_contribution = 0
        if 1 < alpha and alpha < 2:
            initial_value_contribution = initial_values[1] * step_size
        elif 2 < alpha:
            for k in range(1, int(np.ceil(alpha))):
                initial_value_contribution += (
                    initial_values[k]
                    / factorial(k)
                    * (x_points[x_index + 1] ** k - x_points[x_index] ** k)
                )
        elif alpha < 1:
            raise ValueError("Not yet supported!")
        y_prediction[x_index + 1] += initial_value_contribution
        y_prediction[x_index + 1] += y_correction[x_index]
        y_prediction[x_index + 1] += (
            step_size**alpha
            / Gamma(alpha + 1)
            * f_name(x_points[x_index], y_correction[x_index])
        )
        subsum = 0
        for j in range(x_index + 1):
            subsum += PCcoeffs(alpha, j, x_index) * f_name(x_points[j], y_correction[j])
        y_prediction[x_index + 1] += step_size**alpha / Gamma(alpha + 2) * subsum

        y_correction[x_index + 1] += initial_value_contribution
        y_correction[x_index + 1] += y_correction[x_index]
        y_correction[x_index + 1] += (
            step_size**alpha
            / Gamma(alpha + 2)
            * alpha
            * f_name(x_points[x_index], y_correction[x_index])
        )
        y_correction[x_index + 1] += (
            step_size**alpha
            / Gamma(alpha + 2)
            * f_name(x_points[x_index + 1], y_prediction[x_index + 1])
        )
        y_correction[x_index + 1] += step_size**alpha / Gamma(alpha + 2) * subsum

    return y_correction
