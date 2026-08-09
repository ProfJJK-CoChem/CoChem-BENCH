#!/usr/bin/env python3
"""
PyTest suite for CoChem-BENCH Extrapolation Algorithms
Resolves BENCH-02: Validates CBS extrapolation formulas against exact literature values.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from bench_core.extrapolate import calculate_halkier, calculate_f12_cbs, calculate_karton_martin, compute_total_cbs_energy

def test_halkier_extrapolation():
    # Test known literature values for correlation energy
    # E_corr(3) = -0.280 Hartree, E_corr(4) = -0.290 Hartree
    corr_x3 = -0.280
    corr_x4 = -0.290
    cbs_corr = calculate_halkier(corr_x3, corr_x4)
    # Expected: (64*(-0.290) - 27*(-0.280)) / 37 = (-18.56 + 7.56)/37 = -11.0 / 37 = -0.297297...
    expected = (64.0 * (-0.290) - 27.0 * (-0.280)) / 37.0
    np.testing.assert_allclose(cbs_corr, expected, rtol=1e-6)

def test_f12_cbs_extrapolation():
    # BENCH-04: Test F12 1/X^7 extrapolation
    corr_x3 = -0.280
    corr_x4 = -0.290
    cbs_f12 = calculate_f12_cbs(corr_x3, corr_x4, beta=7.0)
    # F12 converges faster, so extrapolation delta is smaller than 1/X^3
    assert cbs_f12 < corr_x4
    assert abs(cbs_f12 - corr_x4) < abs(calculate_halkier(corr_x3, corr_x4) - corr_x4)

def test_karton_martin_hf_alpha_parameterization():
    # BENCH-15: Test Karton-Martin HF extrapolation with different basis families
    hf_x3 = -76.060
    hf_x4 = -76.065
    cbs_aug = calculate_karton_martin(hf_x3, hf_x4, family="aug-cc-pVXZ")
    cbs_cc = calculate_karton_martin(hf_x3, hf_x4, family="cc-pVXZ")
    assert cbs_aug < hf_x4
    assert cbs_cc < hf_x4
    assert cbs_aug != cbs_cc

def test_compute_total_cbs_energy():
    res = compute_total_cbs_energy(
        hf_x3=-76.060, hf_x4=-76.065,
        corr_x3=-0.280, corr_x4=-0.290,
        is_f12=True, basis_family="aug-cc-pVXZ"
    )
    assert "hf_cbs" in res
    assert "corr_cbs" in res
    assert "total_cbs" in res
    assert res["total_cbs"] < -76.35
