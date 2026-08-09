#!/usr/bin/env python3
"""
PyTest suite for CoChem-BENCH Core Modules
Validates basis mapping, CABS selection, DKH2 injection, charge/multiplicity,
counterpoise, multireference verification, geometry validation, thermochemical aggregation,
and SWMR HDF5 serialization.
"""

import tempfile
import pytest
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from bench_core.basis_manager import map_core_valence_basis_per_element, get_cabs_auxiliary_basis
from bench_core.orca_writer import generate_dlpno_ccsd_f12, generate_counterpoise_input
from bench_core.verifier import verify_multireference_gate, export_provenance_json
from bench_core.dispatcher import validate_geometry, compute_thermochemical_cbs, serialize_cbs_to_hdf5

def test_element_core_valence_basis_mapping():
    # BENCH-07: Light elements H, He mapped to aug-cc-pVTZ, heavy elements mapped to aug-cc-pCVTZ
    symbols = ['H', 'C', 'N', 'O', 'He']
    mapped = map_core_valence_basis_per_element(symbols, base_cv_basis="aug-cc-pCVTZ")
    assert mapped['H'] == "aug-cc-pVTZ"
    assert mapped['He'] == "aug-cc-pVTZ"
    assert mapped['C'] == "aug-cc-pCVTZ"

def test_cabs_auxiliary_basis_selection():
    # BENCH-10: CABS auxiliary basis selection
    assert get_cabs_auxiliary_basis("aug-cc-pVTZ") == "aug-cc-pVTZ-CABS"
    assert get_cabs_auxiliary_basis("aug-cc-pVQZ") == "aug-cc-pVQZ-CABS"

def test_heavy_element_dkh2_injection():
    # BENCH-08 & BENCH-09 & BENCH-14: Test ORCA input generation with heavy element Iodine (Z=53)
    coords = [("I", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.6)]
    inp = generate_dlpno_ccsd_f12(coords, basis="aug-cc-pVTZ", charge=0, mult=1, nprocs=12)
    assert "Relativistic DKH2" in inp
    assert "* xyz 0 1" in inp
    assert "%pal nprocs 12 end" in inp
    assert "CABS" in inp

def test_counterpoise_input_generation():
    # BENCH-05: Boys-Bernardi counterpoise input generation
    frag_a = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.9)]
    frag_b = [("H", 2.0, 0.0, 0.0), ("Cl", 3.2, 0.0, 0.0)]
    cp_res = generate_counterpoise_input(frag_a, frag_b, basis="aug-cc-pVTZ")
    assert "complex_input" in cp_res
    assert "monomer_a_input" in cp_res
    assert "H :" in cp_res["monomer_a_input"] or "Cl:" in cp_res["monomer_a_input"]

def test_multireference_verification():
    # BENCH-06: Multi-tiered diagnostic check
    orca_singlet = "T1 diagnostic : 0.0150\nD1 diagnostic : 0.0300\n"
    res_singlet = verify_multireference_gate(orca_singlet, spin_multiplicity=1)
    assert res_singlet["passed"] is True

    orca_radical = "T1 diagnostic : 0.0350\nD1 diagnostic : 0.0400\n"
    res_radical = verify_multireference_gate(orca_radical, spin_multiplicity=2)
    assert res_radical["passed"] is True

def test_geometry_validation():
    # BENCH-18: Clash detection
    valid_coords = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.95)]
    assert validate_geometry(valid_coords) is True

    clash_coords = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.2)]
    with pytest.raises(ValueError):
        validate_geometry(clash_coords)

def test_thermochemical_cbs_and_hdf5_serialization():
    # BENCH-13 & BENCH-16: Thermochemistry and HDF5 serialization
    thermo = compute_thermochemical_cbs(electronic_cbs=-76.350, zpe=0.02, thermal_h_corr=0.025, thermal_g_corr=0.005)
    assert "h_cbs" in thermo
    assert "g_cbs" in thermo

    with tempfile.TemporaryDirectory() as tmpdir:
        h5_p = Path(tmpdir) / "test_state.h5"
        serialize_cbs_to_hdf5(h5_p, thermo)
        assert h5_p.exists()
