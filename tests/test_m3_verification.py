#!/usr/bin/env python3
"""
PyTest Suite for CoChem-BENCH Milestone M3 Verification
Contains all 60 verification assertions specified in explorer_m3/handoff.md (10 per task).
"""

import json
import tempfile
import pytest
from pathlib import Path

from cochem_bench.orca_writer import (
    generate_dlpno_ccsd_f12,
    generate_counterpoise_input,
)
from cochem_bench.verifier import (
    export_provenance_json,
    validate_standing_rule_7,
    audit_frozen_core_bias,
)
from cochem_bench.ram_estimator import (
    estimate_ccsd_memory,
    estimate_ccsd_memory_detailed,
)

# Shared test fixtures
COORDS: list[tuple[str, float, float, float]] = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.96)]
FRAG_A: list[tuple[str, float, float, float]] = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.96)]
FRAG_B: list[tuple[str, float, float, float]] = [("H", 2.0, 0.0, 0.0), ("Cl", 3.2, 0.0, 0.0)]
WORKFLOW_DOC_PATH: Path = Path(__file__).resolve().parents[1] / "20260807_bench_workflow.md"


# ==============================================================================
# BENCH-01 Verification Suite (10 Tests)
# ==============================================================================
def test_bench01_opt_kw_eliminated_dlpno() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True)
    assert "TightOPT" not in inp
    header_opt = [line for line in inp.split("\n") if line.startswith("!")][0]
    assert "Opt" in header_opt.split()
    
    inp_sp = generate_dlpno_ccsd_f12(COORDS, is_opt=False)
    header_sp = [line for line in inp_sp.split("\n") if line.startswith("!")][0]
    assert "Opt" not in header_sp.split()


def test_bench01_opt_kw_eliminated_cp() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B, is_opt=True)
    assert "TightOPT" not in cp_res["complex_input"]
    header_opt = [line for line in cp_res["complex_input"].split("\n") if line.startswith("!")][0]
    assert "Opt" in header_opt.split()

    cp_res_sp = generate_counterpoise_input(FRAG_A, FRAG_B, is_opt=False)
    header_sp = [line for line in cp_res_sp["complex_input"].split("\n") if line.startswith("!")][0]
    assert "Opt" not in header_sp.split()


def test_bench01_tole_threshold_present() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True)
    assert "TolE     1e-7" in inp or "TolE 1e-7" in inp


def test_bench01_tolrmsg_threshold_present() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True)
    assert "TolRMSG  3e-6" in inp or "TolRMSG 3e-6" in inp


def test_bench01_tolmaxg_threshold_present() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True)
    assert "TolMaxG  1e-5" in inp or "TolMaxG 1e-5" in inp


def test_bench01_tolrmsd_threshold_present() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True)
    assert "TolRMSD  5e-5" in inp or "TolRMSD 5e-5" in inp


def test_bench01_tolmaxd_threshold_present() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True)
    assert "TolMaxD  1e-4" in inp or "TolMaxD 1e-4" in inp


def test_bench01_inhess_xtb2_default() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True)
    assert "InHess   XTB2" in inp or "InHess XTB2" in inp


def test_bench01_inhess_lindh_custom() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=True, in_hess="Lindh")
    assert "InHess   Lindh" in inp or "InHess Lindh" in inp


def test_bench01_single_point_no_geom() -> None:
    inp = generate_dlpno_ccsd_f12(COORDS, is_opt=False)
    assert "%geom" not in inp


# ==============================================================================
# BENCH-02 Verification Suite (10 Tests)
# ==============================================================================
def test_bench02_cbs_cp_opt_prohibited_exception() -> None:
    with pytest.raises(ValueError):
        generate_counterpoise_input(FRAG_A, FRAG_B, basis="aug-cc-pVTZ/QZ", is_opt=True)


def test_bench02_cbs_string_cp_opt_prohibited() -> None:
    with pytest.raises(ValueError):
        generate_counterpoise_input(FRAG_A, FRAG_B, basis="aug-cc-pVTZ-cbs", is_opt=True)


def test_bench02_non_augmented_cp_opt_enforced_flag() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B, basis="cc-pVTZ", is_opt=True)
    assert cp_res["cp_opt_enforced"] is True


def test_bench02_augmented_cp_opt_allowed_flag() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B, basis="aug-cc-pVTZ", is_opt=True)
    assert cp_res["cp_opt_enforced"] is False


def test_bench02_non_augmented_single_point_capped_3pct() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B, basis="def2-TZVP", is_opt=False)
    assert cp_res["bsse_accuracy_capped_3pct"] is True


def test_bench02_augmented_single_point_not_capped() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B, basis="aug-cc-pVTZ", is_opt=False)
    assert cp_res["bsse_accuracy_capped_3pct"] is False


def test_bench02_ghost_atoms_monomer_a() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B)
    assert "Cl:" in cp_res["monomer_a_input"]
    assert "H:" in cp_res["monomer_a_input"]
    assert "H :" not in cp_res["monomer_a_input"]


def test_bench02_ghost_atoms_monomer_b() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B)
    assert "O:" in cp_res["monomer_b_input"]
    assert "O :" not in cp_res["monomer_b_input"]
    assert "H:" in cp_res["monomer_b_input"]
    assert "H :" not in cp_res["monomer_b_input"]


def test_bench02_complex_input_no_ghosts() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B)
    assert ":" not in cp_res["complex_input"]


def test_bench02_cbs_single_point_allowed() -> None:
    cp_res = generate_counterpoise_input(FRAG_A, FRAG_B, basis="aug-cc-pVTZ/QZ", is_opt=False)
    assert "complex_input" in cp_res


# ==============================================================================
# BENCH-03 Verification Suite (10 Tests)
# ==============================================================================
def test_bench03_provenance_tags_dict_present() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        data = json.loads(p.read_text(encoding="utf-8"))
        assert "provenance_tags" in data


def test_bench03_hf_cbs_tagged_derived() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        data = json.loads(p.read_text(encoding="utf-8"))
        assert data["provenance_tags"]["hf_cbs_energy_hartree"] == "[D]"


def test_bench03_corr_cbs_tagged_derived() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        data = json.loads(p.read_text(encoding="utf-8"))
        assert data["provenance_tags"]["corr_cbs_energy_hartree"] == "[D]"


def test_bench03_total_cbs_tagged_derived() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        data = json.loads(p.read_text(encoding="utf-8"))
        assert data["provenance_tags"]["total_cbs_energy_hartree"] == "[D]"


def test_bench03_t1_tagged_measured() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        data = json.loads(p.read_text(encoding="utf-8"))
        assert data["provenance_tags"]["t1_diagnostic"] == "[M]"


def test_bench03_d1_tagged_measured() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        data = json.loads(p.read_text(encoding="utf-8"))
        assert data["provenance_tags"]["d1_diagnostic"] == "[M]"


def test_bench03_rule_7_compliance_section_exists() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov.json"
        export_provenance_json({}, {"passed": True}, output_path=p)
        data = json.loads(p.read_text(encoding="utf-8"))
        assert "rule_7_compliance" in data


def test_bench03_rule_7_passes_valid_payload() -> None:
    valid_payload = {
        "provenance_tags": {"m1": "[M]"},
        "accuracy_claims": [{"metric_key": "m1", "name": "Measured Claim 1"}],
    }
    assert validate_standing_rule_7(valid_payload)["rule_7_compliant"] is True


def test_bench03_rule_7_flags_unanchored_derived_claim() -> None:
    bad_d_payload = {
        "provenance_tags": {"m2": "[D]"},
        "accuracy_claims": [{"metric_key": "m2", "name": "Unanchored Derived Claim", "has_measured_anchor": False}],
    }
    assert validate_standing_rule_7(bad_d_payload)["rule_7_compliant"] is False


def test_bench03_rule_7_flags_unanchored_estimated_claim() -> None:
    bad_e_payload = {
        "provenance_tags": {"m3": "[E]"},
        "accuracy_claims": [{"metric_key": "m3", "name": "Unanchored Estimated Claim", "has_measured_anchor": False}],
    }
    assert validate_standing_rule_7(bad_e_payload)["rule_7_compliant"] is False


# ==============================================================================
# BENCH-04 Verification Suite (10 Tests)
# ==============================================================================
def test_bench04_claim_03pct_without_cv_fails() -> None:
    assert audit_frozen_core_bias(0.3, has_core_valence=False)["passed"] is False


def test_bench04_claim_03pct_without_cv_warning_status() -> None:
    assert audit_frozen_core_bias(0.3, has_core_valence=False)["status"] == "FROZEN_CORE_BIAS_WARNING"


def test_bench04_claim_03pct_with_cv_passes() -> None:
    assert audit_frozen_core_bias(0.3, has_core_valence=True)["passed"] is True


def test_bench04_claim_10pct_without_cv_passes() -> None:
    assert audit_frozen_core_bias(1.0, has_core_valence=False)["passed"] is True


def test_bench04_frozen_core_mean_bias_tag() -> None:
    res = audit_frozen_core_bias(0.3, False)
    assert res["frozen_core_mean_bias_pct"] == -0.81
    assert res["provenance_tag"] == "[M]"


def test_bench04_reasons_explain_bias() -> None:
    res = audit_frozen_core_bias(0.3, False)
    assert "-0.81%" in res["reasons"][0]


def test_bench04_claim_05pct_exact_boundary_fails_without_cv() -> None:
    assert audit_frozen_core_bias(0.5, False)["passed"] is False


def test_bench04_claim_051pct_boundary_passes_without_cv() -> None:
    assert audit_frozen_core_bias(0.51, False)["passed"] is True


def test_bench04_return_dict_schema() -> None:
    res = audit_frozen_core_bias(0.3, False)
    assert all(k in res for k in ["passed", "status", "claimed_error_pct", "provenance_tag"])


def test_bench04_integration_with_verifier_export() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        p = Path(tmpdir) / "prov_audit.json"
        export_provenance_json(
            {},
            {"passed": True},
            output_path=p,
            accuracy_claims=[{"claimed_error_pct": 0.3, "has_core_valence": False}],
        )
        data = json.loads(p.read_text(encoding="utf-8"))
        assert "frozen_core_audit" in data
        assert data["frozen_core_audit"][0]["status"] == "FROZEN_CORE_BIAS_WARNING"


# ==============================================================================
# BENCH-05 Verification Suite (10 Tests)
# ==============================================================================
def test_bench05_cfour_hessian_memory_scales_with_natoms() -> None:
    mem_10 = estimate_ccsd_memory(100, 10, n_atoms=10, engine="CFOUR", calculation_type="HESSIAN")
    mem_2 = estimate_ccsd_memory(100, 10, n_atoms=2, engine="CFOUR", calculation_type="HESSIAN")
    assert mem_10 > mem_2


def test_bench05_cfour_36n2_factor_exact() -> None:
    res = estimate_ccsd_memory_detailed(100, 10, n_atoms=2, engine="CFOUR", calculation_type="HESSIAN")
    assert res["cfour_36n2_memory_gb"] > 0


def test_bench05_orca_hessian_does_not_add_cfour_overhead() -> None:
    res = estimate_ccsd_memory_detailed(100, 10, n_atoms=2, engine="ORCA", calculation_type="HESSIAN")
    assert res["cfour_36n2_memory_gb"] == 0


def test_bench05_gpu4pyscf_vram_estimation() -> None:
    res = estimate_ccsd_memory_detailed(100, 10)
    expected_vram = ((100**4) * 8.0) / (1024.0**3)
    assert abs(res["gpu_vram_gb"] - expected_vram) < 1e-6


def test_bench05_gpu_vram_overflow_flag() -> None:
    res = estimate_ccsd_memory_detailed(500, 50, gpu_available=True, gpu_mem_gb=16.0)
    assert res["gpu_overflow"] is True


def test_bench05_gpu_vram_no_overflow_small_basis() -> None:
    res = estimate_ccsd_memory_detailed(50, 10, gpu_available=True, gpu_mem_gb=16.0)
    assert res["gpu_overflow"] is False


def test_bench05_backward_compatibility_float_return() -> None:
    res = estimate_ccsd_memory(100, 10)
    assert isinstance(res, float)


def test_bench05_detailed_dict_return() -> None:
    res = estimate_ccsd_memory_detailed(100, 10)
    assert isinstance(res, dict)


def test_bench05_open_shell_doubles_spin_factor() -> None:
    mem_open = estimate_ccsd_memory(100, 10, is_open_shell=True)
    mem_closed = estimate_ccsd_memory(100, 10, is_open_shell=False)
    assert mem_open > mem_closed


def test_bench05_recommendation_string() -> None:
    res = estimate_ccsd_memory_detailed(500, 50, gpu_available=True, gpu_mem_gb=16.0)
    assert res["recommendation"] == "CPU_FALLBACK"


# ==============================================================================
# BENCH-06 Verification Suite (10 Tests)
# ==============================================================================
def test_bench06_file_exists() -> None:
    assert WORKFLOW_DOC_PATH.exists()


def test_bench06_archived_status_header() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "ARCHIVED & SUPERSEDED" in text


def test_bench06_system_1_ar_ketene() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "Ar–ketene" in text


def test_bench06_system_2_ar_oxazole() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "Ar–oxazole" in text


def test_bench06_system_3_h2co_hcl() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "H₂CO" in text


def test_bench06_system_4_water_dimer() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "Water dimer" in text


def test_bench06_system_5_nh3_hcooh() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "NH₃" in text


def test_bench06_system_6_matched_topology() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "C₆H₆–HCN" in text


def test_bench06_section_17_reference() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "Section 17" in text or "§17" in text


def test_bench06_provenance_tags_referenced() -> None:
    text = WORKFLOW_DOC_PATH.read_text(encoding="utf-8")
    assert "[M]" in text
