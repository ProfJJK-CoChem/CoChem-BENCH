#!/usr/bin/env python3
"""Zero-mock tests for real ORCA 6.1.1 output parsing, diagnostic audits, CLI, and bench_provider."""

from __future__ import annotations

import importlib.metadata
import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from cochem_bench.cli import main as cli_main
from cochem_bench.orca_parser import parse_orca_energy, parse_orca_output_text
from cochem_bench.provider import OPERATIONS, bench_provider, dispatch_operation, execute_handoff
from cochem_bench.verifier import audit_diagnostics

def _d_drive_path(*parts: str) -> Path:
    win_path = Path("D:\\", *parts)
    if win_path.is_file():
        return win_path
    return Path("/mnt/d", *parts)


H2O_ORCA_OUT = _d_drive_path(
    "__CoChem", "GitHub-Repo", "CoChem-GEOM", "benchmarks", "h2o_t_geom_01", "orca_ccsd_cc_pvtz", "h2o_ccsd_pvtz.out"
)
HCOOH_DME_ORCA_OUT = _d_drive_path("_orca_jobs", "hcooh-dme", "S1_revdsd_opt.out")


def test_parse_real_orca_h2o_ccsd_pvtz() -> None:
    """Verify real ORCA 6.1.1 CCSD/cc-pVTZ output parsing against exact converged values."""
    assert H2O_ORCA_OUT.is_file(), f"Real ORCA benchmark file missing: {H2O_ORCA_OUT}"
    res = parse_orca_energy(H2O_ORCA_OUT)

    assert res.orca_version == "6.1.1"
    np.testing.assert_allclose(res.final_single_point_energy_eh, -76.324556570992, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(res.t1_diagnostic, 0.006523723, rtol=0.0, atol=1e-9)
    np.testing.assert_allclose(res.scf_energy_eh, -76.05713391183889, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(res.correlation_energy_eh, -0.267422659, rtol=0.0, atol=1e-9)
    assert res.scf_converged is True
    assert res.opt_converged is True
    assert res.terminated_normally is True
    assert res.slow_conv_detected is False
    assert res.n_single_points >= 1


def test_parse_real_orca_s1_revdsd_opt() -> None:
    """Verify real ORCA 6.1.1 revDSD-PBEP86-D4 output parsing against exact converged values."""
    assert HCOOH_DME_ORCA_OUT.is_file(), f"Real ORCA benchmark file missing: {HCOOH_DME_ORCA_OUT}"
    res = parse_orca_energy(HCOOH_DME_ORCA_OUT)

    assert res.orca_version == "6.1.1"
    np.testing.assert_allclose(res.final_single_point_energy_eh, -344.364232550375, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(res.correlation_energy_eh, -0.648018417, rtol=0.0, atol=1e-9)
    assert res.scf_converged is True
    assert res.opt_converged is True
    assert res.terminated_normally is True
    assert res.slow_conv_detected is False


def test_audit_diagnostics_real_h2o_and_slowconv_detection() -> None:
    """Verify audit_diagnostics passes on converged H2O and flags SlowConv + high T1."""
    audit_h2o = audit_diagnostics(H2O_ORCA_OUT, spin_multiplicity=1)
    assert audit_h2o["passed"] is True
    assert audit_h2o["t1_status"] == "PASSED"
    assert audit_h2o["d1_status"] == "NOT_REPORTED"
    assert audit_h2o["slow_conv_detected"] is False

    troubled_orca_text = """
    Program Version 6.1.1  -  RELEASE   -
    |  1> ! CCSD(T) cc-pVTZ SlowConv
    SCF CONVERGED AFTER 45 CYCLES
    T1 diagnostic ... 0.031500
    D1 diagnostic ... 0.082000
    FINAL SINGLE POINT ENERGY -76.123456789000
    ****ORCA TERMINATED NORMALLY****
    """
    parsed_bad = parse_orca_output_text(troubled_orca_text)
    assert parsed_bad.slow_conv_detected is True
    audit_bad = audit_diagnostics(troubled_orca_text, spin_multiplicity=1)
    assert audit_bad["passed"] is False
    assert audit_bad["slow_conv_detected"] is True
    assert audit_bad["t1_status"] == "MULTIREFERENCE_WARNING"
    assert audit_bad["d1_status"] == "MULTIREFERENCE_WARNING"


def test_bench_provider_contract_and_all_operations() -> None:
    """Verify bench_provider metadata, entry point discovery, and all 4 canonical operations."""
    meta = bench_provider()
    assert meta["scientific_api_version"] == "cochem.module-provider/1"
    assert meta["module_id"] == "bench"
    assert meta["distribution"] == "cochem-bench"
    assert meta["operations"] == list(OPERATIONS)
    assert meta["operations"] == [
        "cbs_extrapolate",
        "relative_energies",
        "parse_orca_energy",
        "audit_diagnostics",
    ]

    # 1. cbs_extrapolate
    cbs_res = dispatch_operation(
        "cbs_extrapolate",
        {
            "e_scf_x": -76.057170,
            "e_scf_y": -76.064825,
            "e_corr_x": -0.275031,
            "e_corr_y": -0.290322,
            "x": 3,
            "y": 4,
            "basis_family": "cc-pVnZ",
        },
    )
    assert cbs_res["e_total_cbs_eh"] < -76.36

    # 2. relative_energies
    rel_res = dispatch_operation(
        "relative_energies",
        {
            "species": [
                {"name": "S1", "e_cbs_eh": -344.364232, "delta_cv_eh": -0.02},
                {"name": "S2", "e_cbs_eh": -344.360000, "delta_cv_eh": -0.02},
            ],
            "required_corrections": ["delta_cv"],
        },
    )
    assert rel_res["reference_species"] == "S1"
    assert rel_res["entries"][1]["delta_e_kcal_mol"] > 0.0

    # 3. parse_orca_energy
    orca_res = dispatch_operation("parse_orca_energy", {"path": str(H2O_ORCA_OUT)})
    np.testing.assert_allclose(orca_res["final_single_point_energy_eh"], -76.324556570992, atol=1e-12)

    # 4. audit_diagnostics
    audit_res = dispatch_operation("audit_diagnostics", {"path": str(H2O_ORCA_OUT)})
    assert audit_res["passed"] is True

    # Unknown operation raises ValueError
    with pytest.raises(ValueError):
        dispatch_operation("nonexistent_op", {})


def test_entry_point_and_handoff_and_cli() -> None:
    """Verify importlib.metadata entry point, execute_handoff, and CLI subcommands."""
    eps = importlib.metadata.entry_points(group="cochem.modules")
    bench_eps = [ep for ep in eps if ep.name == "bench"]
    assert len(bench_eps) == 1
    loaded = bench_eps[0].load()
    assert loaded is bench_provider
    assert loaded()["module_id"] == "bench"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        handoff_file = tmp / "handoff.json"
        handoff_file.write_text(
            json.dumps(
                {
                    "operation": "parse_orca_energy",
                    "parameters": {"path": str(HCOOH_DME_ORCA_OUT)},
                }
            ),
            encoding="utf-8",
        )
        out = execute_handoff(handoff_file, tmp / "out")
        assert Path(out["output_path"]).is_file()
        np.testing.assert_allclose(
            out["result"]["final_single_point_energy_eh"],
            -344.364232550375,
            atol=1e-12,
        )

    assert cli_main(["info"]) == 0
    assert cli_main(["parse-orca", str(H2O_ORCA_OUT)]) == 0
    assert cli_main(["audit", str(H2O_ORCA_OUT)]) == 0
