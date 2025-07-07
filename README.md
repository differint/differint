- This is a fork of original [differint](https://github.com/differint/differint) project
- There is also a faster version (but more limited) implemented in c++ [diffeintC](https://github.com/iparsw/differintC)


# differintP

**differintP** is a modern, pure Python library for **fractional calculus**—that is, for numerical differentiation and integration of arbitrary (non-integer) order.
It is a high-performance fork and major modernization of the original [differint](https://github.com/differint/differint) library, with many optimizations, expanded features, and a strong emphasis on speed, completeness, and code clarity.

> **Highlights**
> - Orders of differentiation and integration: any real (fractional) value
> - Fast, vectorized, and JIT-accelerated core routines (NumPy + Numba)
> - GPU support (optional, via CuPy)
> - Array and pointwise (endpoint) fractional derivative functions
> - Covers Grünwald-Letnikov, Riemann-Liouville, Caputo, Weyl, Riesz, CRONE, and more
> - Accurate for both smooth and non-smooth functions
> - Extensive tests and clear, modern codebase
> - MIT License

---

## 🚀 Quick Install

```bash
pip install differintP
````

* Optional: for GPU acceleration, also install [CuPy](https://docs.cupy.dev/en/stable/install.html):

  ```bash
  pip install cupy-cuda12x  # or cupy-cuda11x, depending on your CUDA
  ```

---

## ✨ Features

* **Fast Fractional Derivatives**: Grünwald-Letnikov (`GL`, `GLpoint`), Riemann-Liouville (`RL`, `RLpoint`), Caputo (`CaputoL1point`, etc.)
* **Modern Implementations**: All methods rewritten for clarity, speed, and extensibility
* **Advanced/Research Methods**:

  * **Weyl** (Fourier/spectral, periodic)
  * **Riesz** (symmetric, physical applications)
  * **CRONE** (for edge detection and signal processing)
* **GPU-accelerated GL**: `GL_gpu` (if CuPy is installed)
* **High test coverage**: Pytest, with analytical ground truths for common functions

---

## 📚 Usage Example

```python
import numpy as np
from differintP import GL, GLpoint, Weyl, Riesz

# Fractional derivative of order 0.5 of sqrt(x) on [0, 1]
x = np.linspace(0, 1, 100)
df_gl = GL(0.5, np.sqrt, 0, 1, 100)         # Array version
df_point = GLpoint(0.5, np.sqrt, 0, 1, 100) # Endpoint value

# Weyl and Riesz derivatives of sin(x) on [0, 2*pi]
x2 = np.linspace(0, 2*np.pi, 256, endpoint=False)
df_weyl = Weyl(0.5, np.sin, 0, 2*np.pi, 256)
df_riesz = Riesz(0.5, np.sin, 0, 2*np.pi, 256)
```

See the [wiki](https://github.com/iparsw/differintP/wiki) for detailed documentation and advanced use.

---

## 📑 Documentation

* **Function Reference**: See the [Wiki](https://github.com/iparsw/differintP/wiki) for math, options, and examples for all supported methods.
* **Benchmarks**: [BENCHMARK.md](https://github.com/iparsw/differintC/blob/main/BENCHMARK.md) for speed/performance notes.
* **Importing**: Import main algorithms directly, See the [Wiki](https://github.com/iparsw/differintP/wiki/Namespaces).

Of course! Here’s a suggested section for accuracy tests, following your README’s style:

---

## ✅ Accuracy Tests

To verify the correctness and numerical accuracy of all algorithms, see the interactive **accuracy test notebook**:
[**accuracy\_test.ipynb**](https://github.com/iparsw/differintP/blob/master/accuracy_test.ipynb)

This notebook compares numerical results to analytical ground truth for a variety of functions and methods, and serves as a practical demonstration and benchmark reference.

---

## 🧑‍💻 Development and Testing

* **Python ≥ 3.10** required
* All code is in pure Python; C++ not required
* **All tests use [pytest](https://pytest.org/)** for fast, expressive, and modern testing

To run the full test suite:

```bash
pytest tests/
```

---



## 📝 License

**MIT License**
See [LICENSE](./LICENSE).

---

## 🔗 Related Projects

* [differint](https://github.com/differint/differint) (original library, now less maintained)
* [differintC](https://github.com/iparsw/differintC) (C/C++ accelerated version by @iparsw)
---

*For questions or help, open an [issue](https://github.com/iparsw/differintP/issues) or check the Discussions!*

