Differint
---------

This package is used for numerically calculating fractional derivatives and integrals (differintegrals). Various methods are used, which correspond to varying definitions of the differintegral::

To use the Grunwald_Letnikov differintegral of order 0.5 with a user-defined function 'func', do::

.. code-block:: python
    import differint
    fracderiv_f = differint.GL(0.5, func)

You can also specify the endpoints of your domain (default is [0,1]) and the number of points to use (default is 100)::
.. code-block:: python 
    fracderiv_f = differint.GL(0.5, func, 0., 1., 100)
