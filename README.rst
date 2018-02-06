differint
---------

This package is used for numerically calculating fractional derivatives
and integrals (differintegrals). Options for varying definitions of the
differintegral are available, including the Grunwald-Letnikov (GL), the
'improved' Grunwald-Letnikov (GLI), the Riemann-Liouville (RL), and the
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

Included Files
--------------

+------------+--------------+
| Core File  | Description  |
+============+==============+
| differint/ | Contains     |
| differint. | algorithms   |
| py         | for          |
|            | fractional   |
|            | differentiat |
|            | ion          |
|            | and          |
|            | integration. |
+------------+--------------+
| tests/test | Testing      |
| .py        | suite        |
|            | containing   |
|            | all unit     |
|            | tests.       |
+------------+--------------+

Both of the above files have corresponding ``__init__.py`` files.

+-------------+--------------+
| Setup File  | Description  |
+=============+==============+
| .gitignore  | List of      |
|             | files to     |
|             | ignore       |
|             | during       |
|             | ``git``      |
|             | push/pull    |
|             | requests.    |
+-------------+--------------+
| CONTRIBUTIN | Instructions |
| G.md        | for          |
|             | potential    |
|             | contributors |
|             | to the       |
|             | *differint*  |
|             | project.     |
+-------------+--------------+
| LICENSE     | MIT license  |
|             | agreement.   |
+-------------+--------------+
| MANIFEST.in | Selects the  |
|             | README file  |
|             | for          |
|             | uploading to |
|             | PyPi.        |
+-------------+--------------+
| README.md   | This README  |
|             | file.        |
+-------------+--------------+
| README.rst  | This README  |
|             | file in      |
|             | ReStructured |
|             | Text         |
|             | format.      |
+-------------+--------------+
| **init**.py | ``__init__`` |
|             | file for     |
|             | overall      |
|             | package.     |
+-------------+--------------+
| changelog.t | List of      |
| xt          | updates to   |
|             | package.     |
+-------------+--------------+
| setup.py    | Script for   |
|             | downloading  |
|             | package from |
|             | ``pip``.     |
+-------------+--------------+

Example Usage
-------------

Taking a fractional derivative is easy with the *differint* package.
Let's take the 1/2 derivative of the square root function on the
interval [0,1], using the Riemann-Liouville definition of the fractional
derivative.

.. code:: python

    import numpy as np
    import differint.differint as df

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

+----------------+--------+
| Main Function  | Usage  |
+================+========+
| GLpoint        | Comput |
|                | es     |
|                | the GL |
|                | differ |
|                | integr |
|                | al     |
|                | at a   |
|                | point  |
+----------------+--------+
| GL             | Comput |
|                | es     |
|                | the GL |
|                | differ |
|                | integr |
|                | al     |
|                | over   |
|                | an     |
|                | entire |
|                | array  |
|                | of     |
|                | functi |
|                | on     |
|                | values |
|                | using  |
|                | the    |
|                | Fast   |
|                | Fourie |
|                | r      |
|                | Transf |
|                | orm    |
+----------------+--------+
| GLI            | Comput |
|                | es     |
|                | the    |
|                | improv |
|                | ed     |
|                | GL     |
|                | differ |
|                | integr |
|                | al     |
|                | over   |
|                | an     |
|                | entire |
|                | array  |
|                | of     |
|                | functi |
|                | on     |
|                | values |
+----------------+--------+
| RLpoint        | Comput |
|                | es     |
|                | the RL |
|                | differ |
|                | integr |
|                | al     |
|                | at a   |
|                | point  |
+----------------+--------+
| RL             | Comput |
|                | es     |
|                | the RL |
|                | differ |
|                | integr |
|                | al     |
|                | over   |
|                | an     |
|                | entire |
|                | array  |
|                | of     |
|                | functi |
|                | on     |
|                | values |
|                | using  |
|                | matrix |
|                | method |
|                | s      |
+----------------+--------+

+---------------------+--------+
| Auxiliary Function  | Usage  |
+=====================+========+
| isInteger           | Determ |
|                     | ine    |
|                     | if a   |
|                     | number |
|                     | is an  |
|                     | intege |
|                     | r      |
+---------------------+--------+
| checkValues         | Used   |
|                     | to     |
|                     | check  |
|                     | for    |
|                     | valid  |
|                     | algori |
|                     | thm    |
|                     | input  |
|                     | types  |
+---------------------+--------+
| GLIinterpolat       | Define |
|                     | interp |
|                     | olatin |
|                     | g      |
|                     | coeffi |
|                     | cients |
|                     | for    |
|                     | the    |
|                     | improv |
|                     | ed     |
|                     | GL     |
|                     | algori |
|                     | thm    |
+---------------------+--------+
| functionCheck       | Determ |
|                     | ines   |
|                     | if     |
|                     | algori |
|                     | thm    |
|                     | functi |
|                     | on     |
|                     | input  |
|                     | is     |
|                     | callab |
|                     | le     |
|                     | or an  |
|                     | array  |
|                     | of     |
|                     | number |
|                     | s      |
+---------------------+--------+
| test\_func          | Testin |
|                     | g      |
|                     | functi |
|                     | on     |
|                     | for    |
|                     | docstr |
|                     | ing    |
|                     | exampl |
|                     | es     |
+---------------------+--------+
| poch                | Comput |
|                     | es     |
|                     | the    |
|                     | Pochha |
|                     | mmer   |
|                     | symbol |
+---------------------+--------+
| GLcoeffs            | Determ |
|                     | ines   |
|                     | the    |
|                     | convol |
|                     | ution  |
|                     | filter |
|                     | compos |
|                     | ed     |
|                     | of     |
|                     | genera |
|                     | lized  |
|                     | binomi |
|                     | al     |
|                     | coeffi |
|                     | cients |
|                     | used   |
|                     | in the |
|                     | GL     |
|                     | algori |
|                     | thm    |
+---------------------+--------+
| RLcoeffs            | Calcul |
|                     | ates   |
|                     | the    |
|                     | coeffi |
|                     | cients |
|                     | used   |
|                     | in the |
|                     | RLpoin |
|                     | t      |
|                     | and RL |
|                     | algori |
|                     | thms   |
+---------------------+--------+
| RLmatrix            | Determ |
|                     | ines   |
|                     | the    |
|                     | matrix |
|                     | used   |
|                     | in the |
|                     | RL     |
|                     | algori |
|                     | thm    |
+---------------------+--------+

Contribute
----------

To contribute to this project, see the `contributing
guidelines <https://github.com/snimpids/differint/blob/master/CONTRIBUTING.md>`__.

Credits
-------

Baleanu, D., Diethelm, K., Scalas, E., & Trujillo, J.J. (2012).
Fractional Calculus: Models and Numerical Methods. World Scientific.

Oldham, K.B. & Spanier, J. (1974). The Fractional Calculus: Theory and
Applications of Differentiation and Integration to Arbitrary Order.
Academic Press Inc.

License
-------

MIT © `Matthew Adams <2018>`__
