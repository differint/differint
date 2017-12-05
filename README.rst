Differint
---------

.. role:: latex(raw)
   :format: latex

This package is used for numerically calculating fractional derivatives and integrals (differintegrals). Various methods are used, which correspond to varying definitions of the differintegral.

To use the Grunwald-Letnikov differintegral of order 0.5 with a user-defined function 'func', do:

.. code-block:: python

    import differint
    fracderiv_f = differint.GL(0.5, func)

If the function data is specified as a list or array, the user can calculate the GL differintegral with:

.. code-block:: python

   fracderiv_fdata = differint.GL(0.5, f_values)

You can also specify the endpoints of your domain (default is [0,1]) and the number of points to use (default is 100).

.. code-block:: python 

    import differint
    fracderiv_f = differint.GL(0.5, func, 0., 1., 100)

The Riemann-Liouville differintegral is also available, which calculates the differintegral using a trapezoid rule. The RL differintegral can be calculated similarly to the GL differintegral as follows.

.. code-block:: python

   fracderiv_f_RL = differint.RLtrap(0.5, func)
   # Or, if the data is a list or array, do:
   fracderiv_f_RL = differint.RLtrap(0.5, f_values)
