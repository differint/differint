Differint
---------

This package is used for numerically calculating fractional derivatives and integrals (differintegrals). Various methods are used, which correspond to varying definitions of the differintegral.

To use the Gr\"uwald_Letnikov differintegral of order 0.5 with a user-defined function $f(x)$, do:

  >>> import differint
  >>> fracderiv_f = differint.GL(0.5, f)
