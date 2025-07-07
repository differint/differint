# differint changelog (parent package)

## 0.1

2017-11-22

### Added

  - 'GL' - calculate GL differintegral

## 0.2

2017-11-24

### Added

  - 'RLtrap' - calculate RL differintegral using trapezoid rule
  - 'RLcoeff' - calculate coefficients for use in 'RLtrap'
  - 'checkvalues' - check for valid inputs for algorithms

### Modified

  - 'GL' - now works with 'checkvalues' to validate inputs

2018-01-29

### Added

  - 'GLI' - calculate improved GL differintegral
  - Unit tests for all core and auxiliary functions
  - 'GLIinterpolat' calculate interpolation coefficients for GLI algorithm
  - 'isInteger' for integer typechecking
  - 'functionCheck' for checking function input
  - 'test_func' for unittesting
  - 'RL' for calculating RL algorithm over array of function values
  - 'GL' for calculating GL algorithm over array of function values

### Modified

  - 'GL' now uses DFT for faster computation
  - 'GL' changed to 'GLpoint'
  - 'RLtrap' now called 'RLpoint'

2018-05-03

### Added

  - 'CRONE' - calcualte GL differintegral operator using CRONE

2019-12-13

### Added

  - Added readme to github repository
  - Include future print function support

2022-11-13

### Added

  - 'isPositiveInteger' - check if a number is and integer and positive
  - 'PCsolver' - solve initial value fractional ODEs using the predictor-corrector method
  - Caputo differintegral using the L1, L2, and L2C methods

### Modified

 - Changed the way `isInteger` checks the given value
 - resolved bug in 'CRONE' method
 - updated Pochhammer symbol implementation to support a larger domain of inputs
 - 'checkValues' now supports numpy data types
 - 'Gamma' now uses `math.gamma` for real values to speed up computation


# differintP Changelog

## 0.0.x

### 0.0.2 (7/3/2025)

- Optimized core functions: `GL`, `RL`, and `RLpoint`
- Added GPU-accelerated `GL` (`GL_gpu`) via [CuPy](https://cupy.dev/) (optional dependency)
- Modernized package setup and build system

[Benchmark](https://github.com/iparsw/differintC/blob/main/BENCHMARK.md) Result for GL:

| count |  differintP | differint   |
| ----- | ----------- | ----------- |
| 1+e2  | 0.0782 ms   | 0.2159 ms   |
| 1+e3  | 0.0942 ms   | 0.6401 ms   |
| 1+e4  | 0.487 ms    | 3.8757 ms   |
| 1+e5  | 9.1574 ms   | 38.2753 ms  |
| 1+e6  | 97.8466 ms  | 403.2405 ms |
| 1+e7  | *           | *           |


[Benchmark](https://github.com/iparsw/differintC/blob/main/BENCHMARK.md) Result for RL

| count |  differintP | differint   |
| ----- | ----------- | ----------- |
| 1+e2  | 0.9003 ms   | 4.054 ms    |
| 1+e3  | 15.8491 ms  | 418.7378 ms |
| 1+e4  | 1372.8949 ms| 3.8757 ms   |


[Benchmark](https://github.com/iparsw/differintC/blob/main/BENCHMARK.md) Result for RLpoint

| count |  differintP | differint   |
| ----- | ----------- | ----------- |
| 1+e2  | 0.0874 ms   | 0.1039 ms   |
| 1+e3  | 0.1519 ms   | 0.9835 ms   |
| 1+e4  | 1.0374 ms   | 9.7152 ms   |
| 1+e5  | 12.5892 ms  | 108.2792 ms |
| 1+e6  | 126.6893 ms | 995.7722 ms |
| 1+e7  | 1257.0316 ms| *           |


### 0.0.3 (7/4/2025)

* **Redesigned `GLpoint`:**
  Replaced the previous GLpoint implementation with a highly efficient, single-pass C++-style recurrence. The new method uses a Numba JIT-compiled kernel for substantial speed gains, directly computing the Grünwald-Letnikov fractional derivative at the endpoint in one loop over the function values.
  Legacy and alternative endpoint methods (`GLpoint_direct`, `GLpoint_via_GL`) have been moved to `special.py` for reference and comparison.


[Benchmark](https://github.com/iparsw/differintC/blob/main/BENCHMARK.md) Result for GLpoint:

| count |  differintP | differint    | Speedup Factor |
| ----- | ----------- | ------------ | -------------- |
| 1+e2  | 0.0148 ms   | 0.0416 ms    | x2.81          |
| 1+e3  | 0.0276 ms   | 0.4213 ms    | x15.26         |
| 1+e4  | 0.0397 ms   | 3.8376 ms    | x96.66         |
| 1+e5  | 0.3285 ms   | 36.303 ms    | x110.5         |
| 1+e6  | 5.267 ms    | 385.9205 ms  | x73.27         |
| 1+e7  | 51.2645 ms  | 3674.3541 ms | x71.67         |


* **Optimized `GLI` implementation:**

  * Rewrote core algorithm to use a Numba-accelerated helper for the main rolling-window convolution, matching the original literature.
  * Ensured correctness by flipping GL coefficients within the moving-window convolution.
  * Replaced class-based coefficient calculation with direct variable computation for improved clarity and efficiency.
  * All critical loops are now compiled with Numba for major speed improvements, especially at high grid resolutions.
  * Modernized function interface and added robust array conversion and shape checks.

[Benchmark](https://github.com/iparsw/differintC/blob/main/BENCHMARK.md) Result for GLI:

| count |  differintP | differint    | Speedup Factor |
| ----- | ----------- | ------------ | -------------- |
| 1+e2  | 0.0760 ms   | 0.7582 ms    | x6.06          |
| 1+e3  | 2.0336 ms   | 7.391 ms     | x3.63          |
| 1+e4  | 155.7649 ms | 128.3709 ms  | x0.8           |
| 5+e4  | 3895.176 ms | 15181.82 ms  | x3.89          |

- Wiki Pages for
  - `CaputoL1point` `CaputoL2point` `CaputoL2Cpoint` `CaputoFromRLpoint`
  - `CRONE` `GL` `GLpoint` `GLI` `GL_gpu`
  - `RL` `RLpoint` `PCsolver` functions Namespaces

- Aditional test and test restructuring
- Restructuring the Namespaces


### 0.0.4

* **Added Weyl and Riesz Fractional Derivative Functions:**

  * **`Weyl`**: Implements the periodic, right-sided fractional derivative using an FFT-based method. For periodic functions, this operator produces phase-shifted analytical derivatives, and is especially well-suited for spectral methods.
  * **`Riesz`**: Implements the symmetric (two-sided) fractional derivative, also using FFT. For pure sines/cosines, the Riesz operator multiplies each Fourier mode by $-|k|^\alpha$, resulting in a real-valued, sign-flipped output (not a phase shift). Especially relevant in physics and PDEs.
* **Testing Update:**

  * Both `Weyl` and `Riesz` are covered by **unittest** and **pytest**-style tests for compatibility and reliability.
  * The test suite is transitioned from `unittest` to `pytest` for modern, flexible, and more expressive testing.
