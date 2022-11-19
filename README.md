## differint
This package is used for numerically calculating fractional derivatives and integrals (differintegrals). Options for varying definitions of the differintegral are available, including the Grunwald-Letnikov (GL), the 'improved' Grunwald-Letnikov (GLI), the Riemann-Liouville (RL), and the Caputo (L1, L2, and L2C). Through the API, you can compute differintegrals at a point or over an array of function values.

See below for an example of how to use this package, or check out the [wiki](https://github.com/differint/differint/wiki) for references, signatures, and examples for each function.

## Motivation
There is little in the way of readily available, easy-to-use code for numerical fractional calculus. What is currently available are functions that are generally either smart parts of a much larger package, or only offer one numerical algorithm. The *differint* package offers a variety of algorithms for computing differintegrals and several auxiliary functions relating to generalized binomial coefficients.

## Installation
This project requires Python 3+ and NumPy to run.

Installation from the Python Packaging index (https://pypi.python.org/pypi) is simple using pip.

```python
pip install differint
```

## Included Files
Core File | Description
--------- | -----------
differint/differint.py | Contains algorithms for fractional differentiation and integration.
tests/test.py | Testing suite containing all unit tests.

Both of the above files have corresponding `__init__.py` files.

Setup File | Description
---------- | -----------
.gitignore | List of files to ignore during `git` push/pull requests.
CONTRIBUTING.md | Instructions for potential contributors to the *differint* project.
LICENSE | MIT license agreement.
MANIFEST.in | Selects the README file for uploading to PyPi.
README.md | This README file.
README.rst | This README file in ReStructuredText format.
__init__.py | `__init__` file for overall package.
changelog.txt | List of updates to package.
setup.py | Script for downloading package from `pip`.

## Example Usage
Taking a fractional derivative is easy with the *differint* package. Let's take the 1/2 derivative of the square root function on the interval [0,1], using the Riemann-Liouville definition of the fractional derivative.

```python
import numpy as np
import differint.differint as df

def f(x):
   return x**0.5

DF = df.RL(0.5, f)
print(DF)
```

You can also specify the endpoints of the domain and the number of points used as follows.

```python
DF = df.RL(0.5, f, 0, 1, 128)
```

For a description of all functions, their signatures, and more usage examples, see the project's [wiki](https://github.com/differint/differint/wiki).

## Tests
All tests can be run with nose from the command line. Setup will automatically install nose if it is not present on your machine.

```python
python setup.py tests
```

Alternatively, you can run the test script directly.

```python
cd <file_path>/differint/tests/
python test.py
```

## API Reference
In this section we cover the usage of the various functions within the *differint* package.

Main Function | Usage
------------- | -----
[GLpoint](https://github.com/differint/differint/wiki/GLpoint) | Computes the GL differintegral at a point
[GL](https://github.com/differint/differint/wiki/GL) | Computes the GL differintegral over an entire array of function values using the Fast Fourier Transform
[GLI](https://github.com/differint/differint/wiki/GLI) | Computes the improved GL differintegral over an entire array of function values
[CRONE](https://github.com/differint/differint/wiki/CRONE) | Calculates the GL derivative approximation using the CRONE operator.
[RLpoint](https://github.com/differint/differint/wiki/RLpoint) | Computes the RL differintegral at a point
[RL](https://github.com/differint/differint/wiki/RL) | Computes the RL differintegral over an entire array of function values using matrix methods
[CaputoL1point](https://github.com/differint/differint/wiki/CaputoL1point) | Computes the Caputo differintegral at a point using the L1 algorithm
[CaputoL2point](https://github.com/differint/differint/wiki/CaputoL2point) | Computes the Caputo differintegral at a point using the L2 algorithm
[CaputoL2Cpoint](https://github.com/differint/differint/wiki/CaputoL2Cpoint) | Computes the Caputo differintegral at a point using the L2C algorithm
[PCsolver](https://github.com/differint/differint/wiki/PCsolver) | Solves IVPs for fractional ODEs of the form ${}^CD^\alpha[y(x)]=f(x,y(x))$ using the predictor-corrector method

Auxiliary Function | Usage
------------------ | -----
[isInteger](https://github.com/differint/differint/wiki/isInteger) | Determine if a number is an integer
[isPositiveInteger](https://github.com/differint/differint/wiki/isPositiveInteger) | Determine if a number is an integer, and if it is greater than 0
[checkValues](https://github.com/differint/differint/wiki/checkValues) | Used to check for valid algorithm input types
[GLIinterpolat](https://github.com/differint/differint/wiki/GLIinterpolat) | Define interpolating coefficients for the improved GL algorithm
[functionCheck](https://github.com/differint/differint/wiki/functionCheck) | Determines if algorithm function input is callable or an array of numbers
[poch](https://github.com/differint/differint/wiki/poch) | Computes the Pochhammer symbol
[Gamma](https://github.com/differint/differint/wiki/Gamma) | Computes the gamma function, an extension of the factorial to complex numbers
[Beta](https://github.com/differint/differint/wiki/Beta) | Computes the beta function, a function related to the binomial coefficient
[MittagLeffler](https://github.com/differint/differint/wiki/MittagLeffler) | Computes the two parameter Mittag-Leffler function, which is important in the solution of fractional ODEs
[GLcoeffs](https://github.com/differint/differint/wiki/GLcoeffs) | Determines the convolution filter composed of generalized binomial coefficients used in the GL algorithm
[RLcoeffs](https://github.com/differint/differint/wiki/RLcoeffs) | Calculates the coefficients used in the RLpoint and RL algorithms
[RLmatrix](https://github.com/differint/differint/wiki/RLmatrix) | Determines the matrix used in the RL algorithm
[PCcoeffs](https://github.com/differint/differint/wiki/PCcoeffs) | Determines the coefficients used in the PC algorithm

## Contribute
To contribute to this project, see the [contributing guidelines](https://github.com/snimpids/differint/blob/master/CONTRIBUTING.md).

## Credits
Baleanu, D., Diethelm, K., Scalas, E., & Trujillo, J.J. (2012). Fractional Calculus: Models and Numerical Methods. World Scientific.

Oldham, K.B. & Spanier, J. (1974). The Fractional Calculus: Theory and Applications of Differentiation and Integration to Arbitrary Order. Academic Press Inc. 

Karniadakis, G.E.. (2019). Handbook of Fractional Calculus with Applications Volume 3: Numerical Methods. De Gruyter.

## License

MIT © [Matthew Adams](2018)
