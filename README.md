## differint
This package is used for numerically calculating fractional derivatives and integrals (differintegrals). Options for varying definitions of the differintegral are available, including the Grunwald-Letnikov, the 'improved' Grunwald-Letnikov, the Riemann-Liouville, and the Caputo (coming soon!). Through the API, you can compute the fractional derivative at a point or over an array of function values.

## Motivation
There is little in the way of readily available, easy-to-use code for numerical fractional calculus. What is currently available are functions that are generally either smart parts of a much larger package, or only offer one numerical algorithm. The *differint* package offers a variety of algorithms for computing differintegrals and several auxiliary functions relating to generalized binomial coefficients.

## Features
What makes your project stand out?

## Installation
This project requires Python 3+ and NumPy to run.

Installation from the Python Packaging index (https://pypi.python.org/pypi) is simple using pip.

```python
pip install differint
```

## Example Usage
Taking a fractional derivative is easy with the *differint* package. Let's take the 1/2 derivative of the square root function on the interval [0,1], using the Riemann-Liouville definition of the fractional derivative.

```python
import numpy as np
import differint as df

def f(x):
   return x**0.5

DF = df.RL(0.5, f)
print(DF)
```

You can also specify the endpoints of the domain and the number of points used as follows.

```python
DF = df.RL(0.5, f, 0, 1, 128)
```

## API Reference

In this section we cover the usage of the various functions within the *differint* package.

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

## How to use?
If people like your project they’ll want to learn how they can use it. To do so include step by step guide to use your project.

## Contribute

To contribute to this project, see the [contributing guidelines](https://github.com/snimpids/differint/CONTRIBUTING.md).

## Credits
Baleanu, D., Diethelm, K., Scalas, E., & Trujillo, J.J. (2012). Fractional Calculus: Models and Numerical Methods. World Scientific.

Oldham, K.B. & Spanier, J. (1974). The Fractional Calculus: Theory and Applications of Differentiation and Integration to Arbitrary Order. 

## License

MIT © [Matthew Adams](2018)
