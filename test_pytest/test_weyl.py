import pytest

import numpy as np

import sys
import os

# Add the parent directory (containing differintP) to the import path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from differintP.weyl import Weyl, Riesz

class TestWeylAccuray:
    #######################
    # Weyl
    #######################

    def test_Weyl_sin_first_derivative(self):
        """Weyl(1, sin(x)) should approximate cos(x) for periodic domain"""
        alpha = 1
        domain_start = 0.0
        domain_end = 2 * np.pi
        num_points = 128
        x = np.linspace(domain_start, domain_end, num_points, endpoint=False)
        gt = np.cos(x)
        result = Weyl(alpha, np.sin, domain_start, domain_end, num_points)
        max_error = np.max(np.abs(result - gt))
        avg_error = np.mean(np.abs(result - gt))
        assert max_error < 2e-14   # Great Performance
        assert avg_error < 5e-15   # Great Performance

    def test_Weyl_sin_half_derivative(self):
        """Weyl(0.5, sin(x)) should match analytical fractional derivative"""
        alpha = 0.5
        domain_start = 0.0
        domain_end = 2 * np.pi
        num_points = 128
        x = np.linspace(domain_start, domain_end, num_points, endpoint=False)
        gt = np.sin(x + alpha * np.pi / 2)
        result = Weyl(alpha, np.sin, domain_start, domain_end, num_points)
        max_error = np.max(np.abs(result - gt))
        avg_error = np.mean(np.abs(result - gt))
        assert max_error < 1e-14   # Great Performance
        assert avg_error < 1e-15   # Great Performance

    def test_Weyl_sin_half_derivative_middle_third(self):
        """Weyl(0.5, sin(x)) should match analytical fractional derivative (middle third only)"""
        alpha = 0.5
        domain_start = 0.0
        domain_end = 2 * np.pi
        num_points = 1024
        x = np.linspace(domain_start, domain_end, num_points, endpoint=False)
        gt = np.sin(x + alpha * np.pi / 2)
        result = Weyl(alpha, np.sin, domain_start, domain_end, num_points)

        # Indices for the middle third
        start_idx = num_points // 3
        end_idx = 2 * num_points // 3

        mid_result = result[start_idx:end_idx]
        mid_gt = gt[start_idx:end_idx]

        max_error = np.max(np.abs(mid_result - mid_gt))
        avg_error = np.mean(np.abs(mid_result - mid_gt))

        assert max_error < 1e-14   # 6.99e-15
        assert avg_error < 1e-14   # 1.97e-15



#######################
# Riesz
#######################
@pytest.mark.parametrize(
    "alpha,max_error_thresh,avg_error_thresh",
    [
        (0.1, 1e-14, 1e-15),   # Riesz(0.5, sin(x)) -> -sin(x)
        (0.3, 1e-14, 1e-15),   # Riesz(0.5, sin(x)) -> -sin(x)
        (0.5, 1e-14, 1e-15),   # Riesz(0.5, sin(x)) -> -sin(x)
        (0.7, 2e-14, 5e-15),   # Riesz(0.5, sin(x)) -> -sin(x)
        (1.0, 2e-14, 5e-15),   # Riesz(1.0, sin(x)) -> -sin(x)
        (1.5, 1.3e-13, 5e-14),   # Riesz(1.5, sin(x)) -> -sin(x)
    ]
)
def test_Riesz_sin_derivative(alpha, max_error_thresh, avg_error_thresh):
    """Riesz(alpha, sin(x)) should approximate -sin(x) for periodic domain for various alphas."""
    domain_start = 0.0
    domain_end = 2 * np.pi
    num_points = 128
    x = np.linspace(domain_start, domain_end, num_points, endpoint=False)
    gt = -np.sin(x)
    # import Riesz however you do in your test file, e.g.:
    result = Riesz(alpha, np.sin, domain_start, domain_end, num_points) # type: ignore
    max_error = np.max(np.abs(result - gt))
    avg_error = np.mean(np.abs(result - gt))
    assert max_error < max_error_thresh
    assert avg_error < avg_error_thresh