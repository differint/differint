Differint
=========

This package is used for numerically calculating fractional derivatives and integrals (differintegrals). Options for varying definitions of the differintegral are available, including the Grunwald-Letnikov, the 'improved' Grunwald-Letnikov, the Riemann-Liouville, and the Caputo (coming soon!).

Getting started
---------------

This package can be downloaded through the Python Package Index (PyPi) from the command line via pip.

.. code-block:: python

   pip install differint
   
Using the differint Package
---------------------------

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

   fracderiv_f_RL = differint.RLpoint(0.5, func)
   # Or, if the data is a list or array, do:
   fracderiv_f_RL = differint.RLpoint(0.5, f_values)

## differint
This package is used for numerically calculating fractional derivatives and integrals (differintegrals). Options for varying definitions of the differintegral are available, including the Grunwald-Letnikov, the 'improved' Grunwald-Letnikov, the Riemann-Liouville, and the Caputo (coming soon!). Through the API, you can compute the fractional derivative at a point or over an array of function values.

## Motivation
There is little in the way of readily available, easy-to-use code for numerical fractional calculus. The *differint* package offers a variety of algorithms for computing differintegrals and several auxiliary functions relating to generalized binomial coefficients.

## Build status
Build status of continus integration i.e. travis, appveyor etc. Ex. - 

[![Build Status](https://travis-ci.org/akashnimare/foco.svg?branch=master)](https://travis-ci.org/akashnimare/foco)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/github/akashnimare/foco?branch=master&svg=true)](https://ci.appveyor.com/project/akashnimare/foco/branch/master)

## Features
What makes your project stand out?

## Code Example
Show what the library does as concisely as possible, developers should be able to figure out **how** your project solves their problem by looking at the code example. Make sure the API you are showing off is obvious, and that your code is short and concise.

## Installation
Installation from the Python Packaging index (https://pypi.python.org/pypi) is simple using pip.

```python
pip install differint
```

## API Reference

Depending on the size of the project, if it is small and simple enough the reference docs can be added to the README. For medium size to larger projects it is important to at least provide a link to where the API reference docs live.

## Tests
Describe and show how to run the tests with code examples.

## How to use?
If people like your project they’ll want to learn how they can use it. To do so include step by step guide to use your project.

## Contribute

Let people know how they can contribute into your project. A [contributing guideline](https://github.com/zulip/zulip-electron/blob/master/CONTRIBUTING.md) will be a big plus.

## Credits
Give proper credits. This could be a link to any repo which inspired you to build this project, any blogposts or links to people who contrbuted in this project. 

#### Anything else that seems useful

## License
A short snippet describing the license (MIT, Apache etc)

MIT © [Yourname]()
