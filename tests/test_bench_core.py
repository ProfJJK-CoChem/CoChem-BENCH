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

from cochem_bench.basis_manager import map_core_valence_basis_per_element, get_cabs_auxiliary_basis
from cochem_bench.orca_writer import generate_dlpno_ccsd_f12, generate_counterpoise_input
from cochem_bench.verifier import verify_multireference_gate, export_provenance_json
from cochem_bench.dispatcher import validate_geometry, compute_thermochemical_cbs, serialize_cbs_to_hdf5

def test_element_core_valence_basis_mapping() -> None:
    # BENCH-07: Light elements H, He mapped to aug-cc-pVTZ, heavy elements mapped to aug-cc-pCVTZ
    symbols = ['H', 'C', 'N', 'O', 'He']
    mapped = map_core_valence_basis_per_element(symbols, base_cv_basis="aug-cc-pCVTZ")
    assert mapped['H'] == "aug-cc-pVTZ"
    assert mapped['He'] == "aug-cc-pVTZ"
    assert mapped['C'] == "aug-cc-pCVTZ"

def test_cabs_auxiliary_basis_selection() -> None:
    # BENCH-10: CABS auxiliary basis selection
    assert get_cabs_auxiliary_basis("aug-cc-pVTZ") == "aug-cc-pVTZ/C"
    assert get_cabs_auxiliary_basis("aug-cc-pVQZ") == "aug-cc-pVQZ/C"

def test_heavy_element_dkh2_injection() -> None:
    # BENCH-08 & BENCH-09 & BENCH-14: Test ORCA input generation with heavy element Iodine (Z=53)
    coords = [("I", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.6)]
    inp = generate_dlpno_ccsd_f12(coords, basis="aug-cc-pVTZ", charge=0, mult=1, nprocs=12)
    assert "DKH2" in inp
    assert "* xyz 0 1" in inp
    assert "%pal nprocs 12 end" in inp
    assert "CABS" in inp

def test_counterpoise_input_generation() -> None:
    # BENCH-05: Boys-Bernardi counterpoise input generation
    frag_a = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.9)]
    frag_b = [("H", 2.0, 0.0, 0.0), ("Cl", 3.2, 0.0, 0.0)]
    cp_res = generate_counterpoise_input(frag_a, frag_b, basis="aug-cc-pVTZ")
    assert "complex_input" in cp_res
    assert "monomer_a_input" in cp_res
    assert "monomer_b_input" in cp_res
    # Complex input must have all real atoms and no ghost atoms (no ':')
    assert ":" not in cp_res["complex_input"]
    # Monomer A input must contain frag_b ghost atoms
    assert "Cl:" in cp_res["monomer_a_input"]
    # Monomer B input must contain frag_a ghost atoms ("O:" without space before colon)
    assert "O:" in cp_res["monomer_b_input"]

def test_multireference_verification() -> None:
    # BENCH-06: Multi-tiered diagnostic check
    orca_singlet = "T1 diagnostic : 0.0150\nD1 diagnostic : 0.0300\n"
    res_singlet = verify_multireference_gate(orca_singlet, spin_multiplicity=1)
    assert res_singlet["passed"] is True

    orca_radical = "T1 diagnostic : 0.0350\nD1 diagnostic : 0.0400\n"
    res_radical = verify_multireference_gate(orca_radical, spin_multiplicity=2)
    assert res_radical["passed"] is True

def test_geometry_validation() -> None:
    # BENCH-18: Clash detection
    valid_coords = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.95)]
    assert validate_geometry(valid_coords) is True

    clash_coords = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.2)]
    with pytest.raises(ValueError):
        validate_geometry(clash_coords)

def test_thermochemical_cbs_and_hdf5_serialization() -> None:
    # BENCH-13 & BENCH-16: Thermochemistry and HDF5 serialization
    thermo = compute_thermochemical_cbs(electronic_cbs=-76.350, zpe=0.02, thermal_h_corr=0.025, thermal_g_corr=0.005)
    assert "h_cbs" in thermo
    assert "g_cbs" in thermo

    with tempfile.TemporaryDirectory() as tmpdir:
        h5_p = Path(tmpdir) / "test_state.h5"
        serialize_cbs_to_hdf5(h5_p, thermo)
        assert h5_p.exists()

def test_dynamic_maxcore_and_tight_keywords() -> None:
    from cochem_bench.ram_estimator import calculate_dynamic_maxcore
    mc = calculate_dynamic_maxcore(nprocs=4, node_max_gb=16.0)
    assert mc >= 1000

    coords = [("C", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.1)]
    inp = generate_dlpno_ccsd_f12(coords, tight_scf=True, is_opt=True, nprocs=4)
    assert "TightSCF" in inp
    assert "TolMaxG" in inp
    assert "TightOPT" not in inp

def test_heavy_element_z36_core_blocks() -> None:
    # Suggestion 14: Iodine Z=53 heavy element core blocks and -DK basis
    coords = [("I", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.6)]
    inp = generate_dlpno_ccsd_f12(coords, basis="aug-cc-pVTZ", rel_mode="auto")
    assert "%core CVCut 1.0e-5 end" not in inp
    assert "-DK" in inp

def test_scratch_disk_space_exception_safety() -> None:
    from cochem_bench.ram_estimator import check_scratch_disk_space
    # Invalid path causing exception returns False
    res = check_scratch_disk_space(required_gb=10.0, scratch_path=str(Path("/nonexistent/path/for/testing")))
    assert res is False

def test_shm_checksum_verification() -> None:
    from cochem_bench.verifier import verify_openmpi_shm_checksum
    data = b"SampleOpenMPISharedMemoryTensorBufferData"
    res = verify_openmpi_shm_checksum(data)
    assert res["passed"] is True
    assert len(res["checksum"]) == 64

def test_provenance_json_dynamic_version() -> None:
    from cochem_bench.verifier import get_active_orca_version, export_provenance_json
    ver = get_active_orca_version()
    assert "ORCA" in ver
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        assert p.exists()
