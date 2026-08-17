#!/usr/bin/env python3
"""
PyTest suite for CoChem-BENCH Core Modules
Validates basis mapping, CABS selection, DKH2 injection, charge/multiplicity,
counterpoise, geometry validation, thermochemical aggregation,
and SWMR HDF5 serialization.
"""

import tempfile
import pytest
from pathlib import Path
import shutil

from cochem_bench.basis_manager import map_core_valence_basis_per_element, get_cabs_auxiliary_basis
from cochem_bench.orca_writer import generate_dlpno_ccsd_f12, generate_counterpoise_input
from cochem_bench.verifier import export_provenance_json, get_active_orca_version
from cochem_bench.dispatcher import validate_geometry, compute_thermochemical_cbs, serialize_cbs_to_hdf5
from cochem_bench.ram_estimator import calculate_dynamic_maxcore, check_scratch_disk_space

def test_element_core_valence_basis_mapping() -> None:
    # BENCH-07: Light elements H, He mapped to aug-cc-pVTZ, heavy elements mapped to aug-cc-pCVTZ
    symbols: list[str] = ['H', 'C', 'N', 'O', 'He']
    mapped: dict[str, str] = map_core_valence_basis_per_element(symbols, base_cv_basis="aug-cc-pCVTZ")
    assert mapped['H'] == "aug-cc-pVTZ"
    assert mapped['He'] == "aug-cc-pVTZ"
    assert mapped['C'] == "aug-cc-pCVTZ"

def test_cabs_auxiliary_basis_selection() -> None:
    # BENCH-10: CABS auxiliary basis selection
    assert get_cabs_auxiliary_basis("aug-cc-pVTZ") == "aug-cc-pVTZ/OptRI"
    assert get_cabs_auxiliary_basis("aug-cc-pVQZ") == "aug-cc-pVQZ/OptRI"

def test_heavy_element_dkh2_injection() -> None:
    # BENCH-08 & BENCH-09 & BENCH-14: Test ORCA input generation with heavy element Iodine (Z=53)
    coords: list[tuple[str, float, float, float]] = [("I", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.6)]
    inp: str = generate_dlpno_ccsd_f12(coords, basis="aug-cc-pVTZ", charge=0, mult=1, nprocs=12)
    assert "DKH2" in inp
    assert "* xyz 0 1" in inp
    assert "%pal nprocs 12 end" in inp
    assert "CABS" in inp

def test_counterpoise_input_generation() -> None:
    # BENCH-05: Boys-Bernardi counterpoise input generation
    frag_a: list[tuple[str, float, float, float]] = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.9)]
    frag_b: list[tuple[str, float, float, float]] = [("H", 2.0, 0.0, 0.0), ("Cl", 3.2, 0.0, 0.0)]
    cp_res: dict[str, str | float | bool] = generate_counterpoise_input(frag_a, frag_b, basis="aug-cc-pVTZ")
    assert "complex_input" in cp_res
    assert "monomer_a_input" in cp_res
    assert "monomer_b_input" in cp_res
    # Complex input must have all real atoms and no ghost atoms (no ':')
    assert isinstance(cp_res["complex_input"], str)
    assert ":" not in cp_res["complex_input"]
    # Monomer A input must contain frag_b ghost atoms
    assert isinstance(cp_res["monomer_a_input"], str)
    assert "Cl:" in cp_res["monomer_a_input"]
    # Monomer B input must contain frag_a ghost atoms ("O:" without space before colon)
    assert isinstance(cp_res["monomer_b_input"], str)
    assert "O:" in cp_res["monomer_b_input"]



def test_geometry_validation() -> None:
    # BENCH-18: Clash detection
    valid_coords: list[tuple[str, float, float, float]] = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.95)]
    assert validate_geometry(valid_coords) is True

    clash_coords: list[tuple[str, float, float, float]] = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.2)]
    with pytest.raises(ValueError):
        validate_geometry(clash_coords)

def test_thermochemical_cbs_and_hdf5_serialization() -> None:
    # BENCH-13 & BENCH-16: Thermochemistry and HDF5 serialization
    thermo: dict[str, str | float | bool] = compute_thermochemical_cbs(electronic_cbs=-76.350, zpe=0.02, thermal_h_corr=0.025, thermal_g_corr=0.005)
    assert "h_cbs" in thermo
    assert "g_cbs" in thermo

    with tempfile.TemporaryDirectory() as tmpdir:
        h5_p: Path = Path(tmpdir) / "test_state.h5"
        serialize_cbs_to_hdf5(h5_p, thermo)
        assert h5_p.exists()

def test_dynamic_maxcore_and_tight_keywords() -> None:
    mc: int = calculate_dynamic_maxcore(nprocs=4, node_max_gb=16.0)
    assert mc >= 1000

    coords: list[tuple[str, float, float, float]] = [("C", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.1)]
    inp: str = generate_dlpno_ccsd_f12(coords, tight_scf=True, is_opt=True, nprocs=4)
    assert "TightSCF" in inp
    assert "TolMaxG" in inp
    assert "TightOPT" not in inp

def test_heavy_element_z36_core_blocks() -> None:
    # Suggestion 14: Iodine Z=53 heavy element core blocks and -DK basis
    coords: list[tuple[str, float, float, float]] = [("I", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.6)]
    inp: str = generate_dlpno_ccsd_f12(coords, basis="aug-cc-pVTZ", rel_mode="auto")
    assert "%core CVCut 1.0e-5 end" not in inp
    assert "-DK" in inp

def test_scratch_disk_space_exception_safety() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        nonexistent_path: Path = Path(tmpdir) / "nonexistent_subdir"
        with pytest.raises(FileNotFoundError):
            check_scratch_disk_space(required_gb=10.0, scratch_path=str(nonexistent_path))



@pytest.mark.skipif(shutil.which('orca') is None, reason="ORCA not installed")
def test_provenance_json_dynamic_version() -> None:
    ver: str = get_active_orca_version()
    assert "ORCA" in ver
    with tempfile.TemporaryDirectory() as tmpdir:
        p: Path = Path(tmpdir) / "prov.json"
        
        input_dict: dict[str, str | float | bool] = {"method": "DLPNO-CCSD(T)", "basis": "aug-cc-pVTZ"}
        results_dict: dict[str, str | float | bool] = {"energy": -76.350, "passed": True}
        provenance_tags: dict[str, str] = {"energy": "Hartree"}
        
        export_provenance_json(input_dict, results_dict, provenance_tags, output_path=p)
        assert p.exists()
