differint
---------

This package is used for numerically calculating fractional derivatives
and integrals (differintegrals). Options for varying definitions of the
differintegral are available, including the Grunwald-Letnikov (GL), the
‘improved’ Grunwald-Letnikov (GLI), the Riemann-Liouville (RL), and the
Caputo (coming soon!). Through the API, you can compute differintegrals
at a point or over an array of function values.

Motivation
----------

There is little in the way of readily available, easy-to-use code for
numerical fractional calculus. What is currently available are functions
that are generally either smart parts of a much larger package, or only
offer one numerical algorithm. The *differint* package offers a variety
of algorithms for computing differintegrals and several auxiliary
functions relating to generalized binomial coefficients.

Installation
------------

This project requires Python 3+ and NumPy to run.

Installation from the Python Packaging index
(https://pypi.python.org/pypi) is simple using pip.

.. code:: python

    pip install differint

Example Usage
-------------

Taking a fractional derivative is easy with the *differint* package.
Let’s take the 1/2 derivative of the square root function on the
interval [0,1], using the Riemann-Liouville definition of the fractional
derivative.

.. code:: python

    import numpy as np
    import differint as df

    def f(x):
       return x**0.5

    DF = df.RL(0.5, f)
    print(DF)

You can also specify the endpoints of the domain and the number of
points used as follows.

.. code:: python

    DF = df.RL(0.5, f, 0, 1, 128)

Tests
-----

All tests can be run with nose from the command line. Setup will
automatically install nose if it is not present on your machine.

.. code:: python

    python setup.py tests

Alternatively, you can run the test script directly.

.. code:: python

    cd <file_path>/differint/tests/
    python test.py

API Reference
-------------

In this section we cover the usage of the various functions within the
*differint* package.

+---------------------------------------------------+-------------------+
| Main Function                                     | Usage             |
+===================================================+===================+
| GLpoint                                           | Computes the GL   |
|                                                   | differintegral at |
|                                                   | a point           |
+---------------------------------------------------+-------------------+
| GL                                                | Computes the GL   |
|                                                   | differintegral    |
|                                                   | over an entire    |
|                                                   | array of function |
|                                                   | values using the  |
|                                                   | Fast Fourier      |
|                                                   | Transform         |
+---------------------------------------------------+-------------------+
| GLI                                               | Computes the      |
|                                                   | improved GL       |
|                                                   | differintegral    |
|                                                   | over an entire    |
|                                                   | array of function |
|                                                   | values            |
+---------------------------------------------------+-------------------+
| RLpoint                                           | Computes the RL   |
|                                                   | differintegral at |
|                                                   | a point           |
+---------------------------------------------------+-------------------+
| RL                                                | Computes the RL   |
|                                                   | differintegral    |
|                                                   | over an entire    |
|                                                   | array of function |
|                                                   | values using      |
|                                                   | matrix methods    |
+---------------------------------------------------+-------------------+

+-------------------------------------------------------+--------------+
| Auxiliary Function                                    | Usage        |
+=======================================================+==============+
| isInteger                                             | Determine if |
|                                                       | a number is  |
|                                                       | an integer   |
+-------------------------------------------------------+--------------+
| checkValues                                           | Used to      |
|                                                       | check for    |
|                                                       | valid        |
|                                                       | algorithm    |
|                                                       | input types  |
+-------------------------------------------------------+--------------+
| GLIinterpolat                                         | Define       |
|                                                       | interpolatin |
|                                                       | g            |
|                                                       | coefficients |
|                                                       | for the      |
|                                                       | improved GL  |
|                                                       | algorithm    |
+-------------------------------------------------------+--------------+
| functionCheck                                         | Determines   |
|                                                       | if algorithm |
|                                                       | function     |
|                                                       | input is     |
|                                                       | callable or  |
|                                                       | an array of  |
|                                                       | numbers      |
+-------------------------------------------------------+--------------+
| test_func                                             | Testing      |
|                                                       | function for |
|                                                       | docstring    |
|                                                       | examples     |
+-------------------------------------------------------+--------------+
| poch                                                  | Computes the |
|                                                       | Pochhammer   |
|                                                       | symbol       |
+-------------------------------------------------------+--------------+
| GLcoeffs                                              | Determines   |
|                                                       | the          |
|                                                       | convolution  |
|                                                       | filter       |
|                                                       | composed of  |
|                                                       | generalized  |
|                                                       | binomial     |
|                                                       | coefficients |
|                                                       | used in the  |
|                                                       | GL algorithm |
+-------------------------------------------------------+--------------+
| RLcoeffs                                              | Calculates   |
|                                                       | the          |
|                                                       | coefficients |
|                                                       | used in the  |
|                                                       | RLpoint and  |
|                                                       | RL           |
|                                                       | algorithms   |
+-------------------------------------------------------+--------------+
