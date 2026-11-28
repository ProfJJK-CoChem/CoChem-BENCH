"""Zero-mock regression tests for FRIC-10: ORCA 6.1.1 double-hybrid SCF parsing & rel-energy keys."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from cochem_bench.cli import main as cli_main
from cochem_bench.orca_parser import parse_orca_energy, parse_orca_output_text
from cochem_bench.thermochemistry import compute_relative_energies


def test_parse_orca_double_hybrid_does_not_overwrite_scf_with_mp2_total() -> None:
    """FRIC-10: 'MP2 TOTAL ENERGY:' must NOT overwrite SCF/DFT 'Total Energy :'."""
    sample = """
                                 Program Version 6.1.1  -  RELEASE   -
    SCF CONVERGED AFTER   7 CYCLES
    Total Energy       :        -76.21152855600619 Eh           -2073.82112 eV
     RI-MP2 CORRELATION ENERGY:     -0.134235954 Eh
    MP2 TOTAL ENERGY:      -76.345764510 Eh
    Dispersion correction           -0.000164397
    FINAL SINGLE POINT ENERGY       -76.345928907188
    ****ORCA TERMINATED NORMALLY****
    """
    res = parse_orca_output_text(sample)
    np.testing.assert_allclose(res.scf_energy_eh, -76.21152855600619, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(res.correlation_energy_eh, -0.134235954, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(res.final_single_point_energy_eh, -76.345928907188, rtol=0.0, atol=1e-12)
    # SCF + MP2 correlation must equal MP2 TOTAL ENERGY (zero double counting!),
    # and SCF + MP2 correlation + D4 dispersion must equal FINAL SINGLE POINT ENERGY.
    np.testing.assert_allclose(
        res.scf_energy_eh + res.correlation_energy_eh,
        -76.345764510,
        rtol=0.0,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        res.scf_energy_eh + res.correlation_energy_eh + (-0.000164397),
        res.final_single_point_energy_eh,
        rtol=0.0,
        atol=1e-8,
    )


def test_parse_real_s1_revdsd_scf_does_not_double_count() -> None:
    win_p = Path(r"D:\_orca_jobs\hcooh-dme\S1_revdsd_opt.out")
    wsl_p = Path("/mnt/d/_orca_jobs/hcooh-dme/S1_revdsd_opt.out")
    target = win_p if win_p.is_file() else wsl_p
    if not target.is_file():
        pytest.skip("Real S1_revdsd_opt.out fixture unavailable")
    res = parse_orca_energy(target)
    assert res.scf_energy_eh is not None and res.correlation_energy_eh is not None
    assert res.final_single_point_energy_eh is not None
    # scf_energy_eh must be strictly higher (less negative) than final_single_point_energy_eh
    assert res.scf_energy_eh > res.final_single_point_energy_eh + 0.5
    # scf + corr must sit within dispersion (~0.02 Eh) of final_single_point_energy_eh
    assert abs((res.scf_energy_eh + res.correlation_energy_eh) - res.final_single_point_energy_eh) < 0.05


def test_rel_energy_accepts_final_single_point_energy_and_gives_clear_error(tmp_path: Path) -> None:
    """FRIC-10: rel-energy accepts uniform final_single_point_energy_eh and rejects mixed keys or levels of theory."""
    table = compute_relative_energies(
        [
            {
                "name": "conf_A",
                "final_single_point_energy_eh": -176.71734061,
                "method": "revDSD-PBEP86-D4/2021",
                "basis": "def2-TZVPP",
            },
            {
                "name": "conf_B",
                "final_single_point_energy_eh": -176.71500000,
                "method": "revDSD-PBEP86-D4/2021",
                "basis": "def2-TZVPP",
            },
        ]
    )
    assert table.reference_species == "conf_A"
    assert table.base_energy_key == "final_single_point_energy_eh"
    assert table.protocol_signature == {"method": "revDSD-PBEP86-D4/2021", "basis": "def2-TZVPP"}
    assert table.entries[1].delta_e_kcal_mol > 1.4

    with pytest.raises(ValueError, match="missing a finite base energy"):
        compute_relative_energies([{"name": "empty_entry"}])

    with pytest.raises(ValueError, match="non-finite base energy"):
        compute_relative_energies([{"name": "nan_entry", "final_single_point_energy_eh": float("nan")}])

    # Mixing incompatible base energy keys (e.g. e_cbs_eh vs final_single_point_energy_eh) MUST fail closed
    with pytest.raises(ValueError, match="Inconsistent base energy keys"):
        compute_relative_energies(
            [
                {"name": "A", "final_single_point_energy_eh": -10.0},
                {"name": "B", "e_cbs_eh": -9.99},
            ]
        )

    # Mixing incompatible method/basis protocols with the same key MUST fail closed
    with pytest.raises(ValueError, match="Inconsistent level of theory"):
        compute_relative_energies(
            [
                {"name": "A", "final_single_point_energy_eh": -10.0, "method": "B3LYP", "basis": "def2-SVP"},
                {"name": "B", "final_single_point_energy_eh": -9.99, "method": "CCSD", "basis": "cc-pVTZ"},
            ]
        )

    json_file = tmp_path / "species.json"
    json_file.write_text(
        json.dumps(
            [
                {"name": "A", "final_single_point_energy_eh": -10.0},
                {"name": "B", "final_single_point_energy_eh": -9.99},
            ]
        ),
        encoding="utf-8",
    )
    assert cli_main(["rel-energy", str(json_file)]) == 0
