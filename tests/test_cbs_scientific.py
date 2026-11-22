#!/usr/bin/env python3
"""Authentic scientific tests for CBS extrapolation and composite thermochemistry.

Validates:
- Two-point Helgaker & Neese-Valeev CBS extrapolation on published H2O and CO
  cc-pVTZ / cc-pVQZ energies (Halkier et al., CPL 1998, 286, 243).
- Basis family parameter lookups for cc-pVnZ, aug-cc-pVnZ, def2, and ano-pVnZ.
- Core-valence (Delta E_CV) and scalar relativistic (Delta E_DKH2/X2C) composite
  thermochemistry with strict PENDING tracking (never substituting 0.0).
- Exact CODATA unit conversions (Eh, kcal/mol, kJ/mol, cm^-1) in relative energy tables.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cochem_bench.constants import (
    HARTREE_TO_CM_INV,
    HARTREE_TO_KCAL_MOL,
    HARTREE_TO_KJ_MOL,
    get_cbs_parameters,
)
from cochem_bench.extrapolate import cbs_extrapolate
from cochem_bench.thermochemistry import compute_composite_energy, compute_relative_energies


def test_h2o_published_helgaker_cbs_tq() -> None:
    """Verify two-point Helgaker (3,4) CBS extrapolation on published H2O CCSD(T) energies.

    Published H2O cc-pVTZ / cc-pVQZ reference values (Halkier et al., CPL 1998, 286, 243):
    E_HF(TZ)   = -76.057170 Eh,  E_HF(QZ)   = -76.064825 Eh
    E_corr(TZ) = -0.275031 Eh,   E_corr(QZ) = -0.290322 Eh  (CCSD(T) valence)
    """
    e_hf_tz, e_hf_qz = -76.057170, -76.064825
    e_corr_tz, e_corr_qz = -0.275031, -0.290322

    res = cbs_extrapolate(
        e_scf_x=e_hf_tz,
        e_scf_y=e_hf_qz,
        e_corr_x=e_corr_tz,
        e_corr_y=e_corr_qz,
        x=3,
        y=4,
        basis_family="cc-pVnZ",
        scheme="helgaker",
    )

    # Analytical Helgaker X^-3 correlation limit: (64*E4 - 27*E3) / 37
    expected_corr = (64.0 * e_corr_qz - 27.0 * e_corr_tz) / 37.0
    # Analytical exponential SCF limit with alpha = 5.46 (cc-pVnZ (3,4))
    alpha = 5.46
    expected_hf = e_hf_qz + (e_hf_qz - e_hf_tz) / math.expm1(alpha * (math.sqrt(4.0) - math.sqrt(3.0)))

    np.testing.assert_allclose(res.e_corr_cbs_eh, expected_corr, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(res.e_scf_cbs_eh, expected_hf, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(res.e_total_cbs_eh, expected_hf + expected_corr, rtol=1e-12, atol=1e-12)
    # Physical sanity: H2O HF/CBS ~ -76.067 Eh, CCSD(T)/CBS ~ -76.369 Eh
    assert -76.069 < res.e_scf_cbs_eh < -76.065
    assert -0.305 < res.e_corr_cbs_eh < -0.295
    assert res.e_total_cbs_eh < (e_hf_qz + e_corr_qz)


def test_co_published_helgaker_cbs_tq() -> None:
    """Verify two-point Helgaker (3,4) CBS extrapolation on published CO CCSD(T) energies.

    Published CO cc-pVTZ / cc-pVQZ reference values (Halkier et al., CPL 1998, 286, 243):
    E_HF(TZ)   = -112.781020 Eh, E_HF(QZ)   = -112.789212 Eh
    E_corr(TZ) = -0.376550 Eh,   E_corr(QZ) = -0.398210 Eh
    """
    e_hf_tz, e_hf_qz = -112.781020, -112.789212
    e_tot_tz = e_hf_tz + (-0.376550)
    e_tot_qz = e_hf_qz + (-0.398210)

    # Test total-energy input mode (derives E_corr = E_tot - E_SCF)
    res = cbs_extrapolate(
        e_scf_x=e_hf_tz,
        e_scf_y=e_hf_qz,
        e_tot_x=e_tot_tz,
        e_tot_y=e_tot_qz,
        x=3,
        y=4,
        basis_family="cc-pVnZ",
        scheme="helgaker",
    )
    expected_corr = (64.0 * (-0.398210) - 27.0 * (-0.376550)) / 37.0
    np.testing.assert_allclose(res.e_corr_cbs_eh, expected_corr, rtol=1e-12, atol=1e-12)
    assert res.e_scf_cbs_eh < e_hf_qz
    assert res.e_total_cbs_eh < e_tot_qz


@pytest.mark.parametrize(
    ("family", "pair", "expected_alpha", "expected_nv_beta"),
    [
        ("cc-pVnZ", (2, 3), 4.42, 2.46),
        ("cc-pVnZ", (3, 4), 5.46, 3.05),
        ("aug-cc-pVnZ", (3, 4), 5.79, 3.05),
        ("def2", (2, 3), 10.39, 2.40),
        ("def2", (3, 4), 7.88, 2.97),
        ("ano-pVnZ", (2, 3), 5.41, 2.43),
        ("ano-pVnZ", (3, 4), 5.76, 2.90),
    ],
)
def test_cbs_parameter_tables_and_neese_valeev_scheme(
    family: str,
    pair: tuple[int, int],
    expected_alpha: float,
    expected_nv_beta: float,
) -> None:
    params = get_cbs_parameters(family, pair[0], pair[1])
    assert params.alpha == expected_alpha
    assert params.beta_neese_valeev == expected_nv_beta

    res = cbs_extrapolate(
        e_scf_x=-76.040,
        e_scf_y=-76.060,
        e_corr_x=-0.260,
        e_corr_y=-0.285,
        x=pair[0],
        y=pair[1],
        basis_family=family,
        scheme="neese_valeev",
    )
    assert res.alpha == expected_alpha
    assert res.beta == expected_nv_beta
    assert res.e_scf_cbs_eh < -76.060
    assert res.e_corr_cbs_eh < -0.285


def test_composite_energy_cv_and_relativistic_additivity() -> None:
    """Verify Delta E_CV and Delta E_DKH2 additive composite energy calculation."""
    res = compute_composite_energy(
        name="H2O_eq",
        e_cbs_eh=-76.369000,
        e_ae_cv_eh=-76.412500,
        e_fc_cv_eh=-76.357500,  # delta_cv = -0.055000 Eh
        e_rel_eh=-76.412650,
        e_nonrel_eh=-76.362000,  # delta_rel = -0.050650 Eh
        rel_method="DKH2",
        delta_zpe_eh=0.021350,
        required_corrections=("delta_cv", "delta_rel", "delta_zpe"),
    )
    assert res.status == "COMPLETE"
    assert res.correction_status["delta_cv"] == "PRESENT"
    assert res.correction_status["delta_rel"] == "PRESENT"
    assert res.correction_status["delta_zpe"] == "PRESENT"
    assert res.correction_status["delta_thermal_h"] == "PENDING"
    np.testing.assert_allclose(res.delta_cv_eh, -0.055000, atol=1e-12)
    np.testing.assert_allclose(res.delta_rel_eh, -0.050650, atol=1e-12)
    expected_total = -76.369000 - 0.055000 - 0.050650 + 0.021350
    np.testing.assert_allclose(res.e_composite_eh, expected_total, atol=1e-12)


def test_composite_energy_missing_required_correction_is_pending_not_zero() -> None:
    """Missing required corrections MUST stay PENDING and NEVER be coerced to 0.0."""
    res = compute_composite_energy(
        name="conformer_B",
        e_cbs_eh=-344.364232,
        delta_cv_eh=-0.012000,
        delta_rel_eh=None,  # Missing DKH2/X2C correction!
        required_corrections=("delta_cv", "delta_rel"),
    )
    assert res.status == "PENDING"
    assert res.e_composite_eh is None
    assert res.delta_rel_eh is None
    assert res.correction_status["delta_rel"] == "PENDING"
    assert "delta_rel" in res.pending_corrections


def test_relative_energy_ranking_table_codata_and_pending() -> None:
    """Verify ranking, exact CODATA conversions, and PENDING preservation."""
    table = compute_relative_energies(
        species=[
            {
                "name": "isomer_B",
                "e_cbs_eh": -100.000000,
                "delta_cv_eh": -0.010000,
                "delta_rel_eh": -0.002000,
            },
            {
                "name": "isomer_A",
                "e_cbs_eh": -100.005000,
                "delta_cv_eh": -0.010000,
                "delta_rel_eh": -0.002000,
            },
            {
                "name": "isomer_C_incomplete",
                "e_cbs_eh": -100.003000,
                "delta_cv_eh": -0.010000,
                # missing delta_rel_eh!
            },
        ],
        required_corrections=("delta_cv", "delta_rel"),
    )

    assert table.reference_species == "isomer_A"
    assert len(table.entries) == 3

    rank1, rank2, pending = table.entries
    assert rank1.rank == 1 and rank1.name == "isomer_A"
    assert rank1.delta_e_eh == 0.0
    assert rank1.delta_e_kcal_mol == 0.0

    assert rank2.rank == 2 and rank2.name == "isomer_B"
    np.testing.assert_allclose(rank2.delta_e_eh, 0.005, atol=1e-12)
    np.testing.assert_allclose(rank2.delta_e_kcal_mol, 0.005 * HARTREE_TO_KCAL_MOL, rtol=1e-12)
    np.testing.assert_allclose(rank2.delta_e_kj_mol, 0.005 * HARTREE_TO_KJ_MOL, rtol=1e-12)
    np.testing.assert_allclose(rank2.delta_e_cm_inv, 0.005 * HARTREE_TO_CM_INV, rtol=1e-12)

    assert pending.rank is None
    assert pending.name == "isomer_C_incomplete"
    assert pending.status == "PENDING"
    assert pending.delta_e_eh is None
    assert pending.delta_e_kcal_mol is None
    assert "delta_rel" in pending.pending_corrections
