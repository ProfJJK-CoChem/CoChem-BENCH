#!/usr/bin/env python3
"""
CoChem-BENCH: Diagnostic Verifier & SCRIBE Provenance Logger
Evaluates multireference T1 and D1 diagnostics differentiated by spin multiplicity
and exports structured JSON provenance logs for SCRIBE.
"""

import re
import json
import shutil
import hashlib
import logging
import subprocess
from pathlib import Path

logger = logging.getLogger("CoChem-BENCH.verifier")

def get_active_orca_version() -> str:
    """
    MOCK-04 / Suggestion 12: Probes host system or project configuration for actual ORCA binary version.
    """
    orca_bin = shutil.which("orca")
    if orca_bin:
        try:
            res = subprocess.run([orca_bin, "--version"], capture_output=True, text=True, timeout=5)
            output = (res.stdout or res.stderr or "").strip()
            match = re.search(r'(\d+\.\d+(\.\d+)?)', output)
            if match:
                return f"ORCA {match.group(1)}"
        except Exception as e:
            logger.debug(f"ORCA binary version check failed: {e}")
            
    root_cfg = Path(__file__).resolve().parents[2] / "cochem_system_config.json"
    if root_cfg.exists():
        try:
            with open(root_cfg, "r", encoding="utf-8") as f:
                cfg = json.load(f)
                v = cfg.get("orca_version") or cfg.get("engines", {}).get("orca", {}).get("version")
                if v and v not in ("auto", "None"):
                    return f"ORCA {v}"
        except Exception as e:
            logger.debug(f"Config orca_version check failed: {e}")
            
    return "ORCA 6.1.1 (Auto Fallback)"

def parse_multireference_diagnostics(orca_output_text: str) -> dict:
    """
    Parses T1 and D1(CCSD) diagnostic values from ORCA log file output.
    Parses the final converged occurrence in multi-step files.
    """
    diagnostics = {"t1_diagnostic": None, "d1_diagnostic": None}
    
    t1_matches = re.findall(r'T1\s+diagnostic\s+:\s+([\d\.-]+)', orca_output_text, re.IGNORECASE)
    if t1_matches:
        diagnostics["t1_diagnostic"] = float(t1_matches[-1])

    d1_matches = re.findall(r'D1\s+diagnostic\s+:\s+([\d\.-]+)', orca_output_text, re.IGNORECASE)
    if d1_matches:
        diagnostics["d1_diagnostic"] = float(d1_matches[-1])

    return diagnostics

def verify_multireference_gate(orca_output_text: str, spin_multiplicity: int = 1) -> dict:
    """
    Resolves BENCH-06 & Suggestion 16: Multi-tiered diagnostic verification evaluating T1 and D1(CCSD)
    differentiated by spin multiplicity 2S+1.
    - Singlets (2S+1=1): T1 > 0.02 or D1 > 0.05 flags multireference character.
    - Open-shell radicals (2S+1 >= 2): T1 up to 0.045 is acceptable if D1 < 0.05.
    Sets 'status': 'MULTIREFERENCE_WARNING' and 'extrapolation_valid': False when thresholds are exceeded.
    """
    diag = parse_multireference_diagnostics(orca_output_text)
    t1 = diag["t1_diagnostic"]
    d1 = diag["d1_diagnostic"]

    if t1 is None:
        return {
            "passed": True,
            "is_multireference": False,
            "extrapolation_valid": True,
            "status": "VALIDATED",
            "t1": 0.0,
            "d1": 0.0,
            "reason": "No T1 reported"
        }

    if spin_multiplicity == 1:
        t1_threshold = 0.02
        d1_threshold = 0.05
    else:
        t1_threshold = 0.045
        d1_threshold = 0.05

    is_multireference = False
    reasons = []

    if t1 > t1_threshold:
        is_multireference = True
        reasons.append(f"T1 diagnostic ({t1:.4f}) exceeds threshold ({t1_threshold:.4f}) for spin multiplicity {spin_multiplicity}")

    if d1 is not None and d1 > d1_threshold:
        is_multireference = True
        reasons.append(f"D1 diagnostic ({d1:.4f}) exceeds threshold ({d1_threshold:.4f})")

    status_str = "MULTIREFERENCE_WARNING" if is_multireference else "VALIDATED"
    result = {
        "passed": not is_multireference,
        "is_multireference": is_multireference,
        "extrapolation_valid": not is_multireference,
        "status": status_str,
        "t1": t1,
        "d1": d1 if d1 is not None else 0.0,
        "spin_multiplicity": spin_multiplicity,
        "reasons": reasons
    }

    if is_multireference:
        logger.warning(f"Multireference gate triggered: {reasons}")
        
    return result

def parse_orca_output_file(output_path: Path) -> dict:
    """
    Parses an ORCA output file and returns diagnostics and multireference verification.
    """
    if not output_path.exists():
        raise FileNotFoundError(f"ORCA output file not found: {output_path}")
    text = output_path.read_text(encoding="utf-8", errors="replace")
    mult = 1
    mult_match = re.search(r'Multiplicity\s+::\s+(\d+)', text)
    if mult_match:
        mult = int(mult_match.group(1))
    return verify_multireference_gate(text, spin_multiplicity=mult)

def verify_openmpi_shm_checksum(data_buffer, expected_checksum: str = None) -> dict:
    """
    Suggestion 17: Performs byte-level SHA-256 checksum verification of OpenMPI shared-memory tensor buffers
    to detect silent bus corruption or memory tearing (§1.2.8).
    """
    if not isinstance(data_buffer, (bytes, bytearray)):
        if hasattr(data_buffer, "tobytes"):
            data_buffer = data_buffer.tobytes()
        elif hasattr(data_buffer, "buffer"):
            data_buffer = bytes(data_buffer.buffer)
        else:
            data_buffer = bytes(data_buffer)

    actual_checksum = hashlib.sha256(data_buffer).hexdigest()
    passed = True
    err_msg = ""
    
    if expected_checksum and actual_checksum.lower() != expected_checksum.lower():
        passed = False
        err_msg = f"SHM Checksum Mismatch: expected {expected_checksum[:16]}..., got {actual_checksum[:16]}..."
        logger.error(err_msg)

    return {
        "passed": passed,
        "checksum": actual_checksum,
        "error": err_msg,
        "extrapolation_valid": passed
    }

def audit_frozen_core_bias(claimed_error_pct: float, has_core_valence: bool, basis_set: str = "") -> dict:
    """
    Audits benchmark accuracy claims against Frozen-Core Systematic Bias statistics (§4.8).
    fc-CCSD(T) carries -0.81% mean bias [M] in B_e. Any claim <= 0.5% error requires core-valence evaluation.
    """
    FROZEN_CORE_MEAN_BIAS_PCT = -0.81  # [M]
    
    requires_cv = claimed_error_pct <= 0.5
    passed = not (requires_cv and not has_core_valence)
    
    status = "VALIDATED" if passed else "FROZEN_CORE_BIAS_WARNING"
    reasons = []
    if not passed:
        reasons.append(
            f"Claimed error {claimed_error_pct:.2f}% <= 0.5% in B_e requires core-valence correlation evaluation (cc-pCVTZ). "
            f"fc-CCSD(T) carries a -0.81% systematic mean bias [M] (§4.8)."
        )
        
    return {
        "passed": passed,
        "status": status,
        "claimed_error_pct": claimed_error_pct,
        "has_core_valence": has_core_valence,
        "frozen_core_mean_bias_pct": FROZEN_CORE_MEAN_BIAS_PCT,
        "provenance_tag": "[M]",
        "reasons": reasons
    }

def validate_standing_rule_7(payload: dict) -> dict:
    """
    Audits JSON payload for Rule 7 compliance (§12.5):
    No derived [D] or estimated [E] value may solely support a hardware exclusion or accuracy claim.
    """
    violations = []
    tags = payload.get("provenance_tags", {})
    claims = payload.get("accuracy_claims", []) + payload.get("hardware_exclusions", [])
    
    for claim in claims:
        key = claim.get("metric_key")
        tag = tags.get(key, claim.get("tag"))
        if tag in ["[D]", "[E]"] and not claim.get("has_measured_anchor", False):
            violations.append(
                f"Rule 7 Violation: Claim '{claim.get('name')}' relies solely on tag {tag} for metric '{key}' without [M] measured anchor."
            )
            
    is_compliant = len(violations) == 0
    return {
        "rule_7_compliant": is_compliant,
        "violations": violations
    }

def export_provenance_json(
    cbs_results: dict, 
    verification_results: dict, 
    output_path: Path = Path("bench_provenance.json"), 
    software_version: str = None,
    accuracy_claims: list = None,
    hardware_exclusions: list = None,
    has_core_valence: bool = False
) -> Path:
    """
    Resolves BENCH-17 & MOCK-04 & BENCH-03/04: Exports structured JSON provenance logging for SCRIBE.
    Contains dynamic software version (ORCA binary check), basis set citations, extrapolation parameters,
    diagnostic metrics, provenance tags, Rule 7 compliance, and frozen-core bias audits.
    """
    soft_ver = software_version or get_active_orca_version()
    status_val = verification_results.get("status")
    if not status_val:
        status_val = "VALIDATED" if verification_results.get("passed", True) else "MULTIREFERENCE_WARNING"

    provenance_tags = {
        "hf_cbs_energy_hartree": "[D]",
        "corr_cbs_energy_hartree": "[D]",
        "total_cbs_energy_hartree": "[D]",
        "t1_diagnostic": "[M]",
        "d1_diagnostic": "[M]",
        "status": "[M]"
    }

    provenance = {
        "software": soft_ver,
        "citation_dispersion": "Grimme D4",
        "cbs_extrapolation": {
            "hf_method": "Karton-Martin Exponential [D]",
            "corr_method": "F12 Power-Law (1/X^7) [D]" if cbs_results.get("is_f12", True) else "1/X^3 Halkier [D]",
            "basis_family": cbs_results.get("basis_family", "aug-cc-pVXZ"),
            "hf_cbs_energy_hartree": cbs_results.get("hf_cbs", 0.0),
            "corr_cbs_energy_hartree": cbs_results.get("corr_cbs", 0.0),
            "total_cbs_energy_hartree": cbs_results.get("total_cbs", 0.0)
        },
        "multireference_diagnostics": verification_results,
        "provenance_tags": provenance_tags,
        "status": status_val
    }

    payload_for_rule7 = {
        "provenance_tags": provenance_tags,
        "accuracy_claims": accuracy_claims or [],
        "hardware_exclusions": hardware_exclusions or []
    }
    provenance["rule_7_compliance"] = validate_standing_rule_7(payload_for_rule7)

    # Evaluate frozen core bias audit if accuracy claims are provided
    frozen_core_audits = []
    for claim in (accuracy_claims or []):
        err_pct = claim.get("claimed_error_pct")
        if err_pct is not None:
            cv_flag = claim.get("has_core_valence", has_core_valence)
            aud_res = audit_frozen_core_bias(err_pct, cv_flag, basis_set=cbs_results.get("basis_family", ""))
            frozen_core_audits.append(aud_res)
    if frozen_core_audits:
        provenance["frozen_core_audit"] = frozen_core_audits

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=4)

    logger.info(f"Exported BENCH provenance JSON to {output_path.absolute()}")
    return output_path
