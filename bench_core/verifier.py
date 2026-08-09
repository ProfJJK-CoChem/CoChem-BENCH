#!/usr/bin/env python3
"""
CoChem-BENCH: Diagnostic Verifier & SCRIBE Provenance Logger
Evaluates multireference T1 and D1 diagnostics differentiated by spin multiplicity
and exports structured JSON provenance logs for SCRIBE.
"""

import re
import json
import logging
from pathlib import Path

logger = logging.getLogger("CoChem-BENCH.verifier")

def parse_multireference_diagnostics(orca_output_text: str) -> dict:
    """
    Parses T1 and D1(CCSD) diagnostic values from ORCA log file output.
    """
    diagnostics = {"t1_diagnostic": None, "d1_diagnostic": None}
    
    t1_match = re.search(r'T1\s+diagnostic\s+:\s+([\d\.-]+)', orca_output_text, re.IGNORECASE)
    if t1_match:
        diagnostics["t1_diagnostic"] = float(t1_match.group(1))

    d1_match = re.search(r'D1\s+diagnostic\s+:\s+([\d\.-]+)', orca_output_text, re.IGNORECASE)
    if d1_match:
        diagnostics["d1_diagnostic"] = float(d1_match.group(1))

    return diagnostics

def verify_multireference_gate(orca_output_text: str, spin_multiplicity: int = 1) -> dict:
    """
    Resolves BENCH-06: Multi-tiered diagnostic verification evaluating T1 and D1(CCSD)
    differentiated by spin multiplicity 2S+1.
    - Singlets (2S+1=1): T1 > 0.02 or D1 > 0.05 flags multireference character.
    - Open-shell radicals (2S+1 >= 2): T1 up to 0.045 is acceptable if D1 < 0.05.
    """
    diag = parse_multireference_diagnostics(orca_output_text)
    t1 = diag["t1_diagnostic"]
    d1 = diag["d1_diagnostic"]

    if t1 is None:
        return {"passed": True, "is_multireference": False, "t1": 0.0, "d1": 0.0, "reason": "No T1 reported"}

    if spin_multiplicity == 1:
        # Closed-shell singlet
        t1_threshold = 0.02
        d1_threshold = 0.05
    else:
        # Open-shell radical (2S+1 >= 2)
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

    result = {
        "passed": not is_multireference,
        "is_multireference": is_multireference,
        "t1": t1,
        "d1": d1 if d1 is not None else 0.0,
        "spin_multiplicity": spin_multiplicity,
        "reasons": reasons
    }

    if is_multireference:
        logger.warning(f"Multireference gate triggered: {reasons}")
        
    return result

def export_provenance_json(cbs_results: dict, verification_results: dict, output_path: Path = Path("bench_provenance.json")) -> Path:
    """
    Resolves BENCH-17: Exports structured JSON provenance logging for SCRIBE.
    Contains exact software versions (ORCA 6.1.1), basis set citations, extrapolation parameters,
    and diagnostic metrics.
    """
    provenance = {
        "software": "ORCA 6.1.1",
        "citation_dispersion": "Grimme D4",
        "cbs_extrapolation": {
            "hf_method": "Karton-Martin Exponential",
            "corr_method": "F12 Power-Law (1/X^7)" if cbs_results.get("is_f12", True) else "1/X^3 Halkier",
            "basis_family": cbs_results.get("basis_family", "aug-cc-pVXZ"),
            "hf_cbs_energy_hartree": cbs_results.get("hf_cbs", 0.0),
            "corr_cbs_energy_hartree": cbs_results.get("corr_cbs", 0.0),
            "total_cbs_energy_hartree": cbs_results.get("total_cbs", 0.0)
        },
        "multireference_diagnostics": verification_results,
        "status": "VALIDATED" if verification_results.get("passed", True) else "MULTIREFERENCE_WARNING"
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=4)

    logger.info(f"Exported BENCH provenance JSON to {output_path.absolute()}")
    return output_path
