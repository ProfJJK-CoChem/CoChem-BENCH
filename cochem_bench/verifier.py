"""CoChem-BENCH Verifier & Compliance Engine.

Validates multi-reference diagnostics (T1/D1), SlowConv/convergence audits,
Rule 7 provenance discipline, frozen-core bias audits, and OpenMPI SHM checksums.
"""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from .orca_parser import OrcaEnergyResult, parse_orca_energy, parse_orca_output_text

DEFAULT_PROVENANCE_TAGS: dict[str, str] = {
    "hf_cbs_energy_hartree": "[D]",
    "corr_cbs_energy_hartree": "[D]",
    "total_cbs_energy_hartree": "[D]",
    "t1_diagnostic": "[M]",
    "d1_diagnostic": "[M]",
}


class ProvenancePayload(BaseModel):
    input: dict[str, Any]
    results: dict[str, Any]
    provenance_tags: dict[str, str]
    accuracy_claims: list[dict[str, Any]]
    orca_version: str
    rule_7_compliance: dict[str, Any] = Field(default_factory=dict)
    frozen_core_audit: list[dict[str, Any]] = Field(default_factory=list)


def verify_multireference_gate(orca_output: str, spin_multiplicity: int = 1) -> dict[str, Any]:
    """Evaluate T1 and D1 diagnostic values against threshold gates."""
    t1_match = re.search(r"T1\s+diagnostic\s*(?:\.+|:)\s*([0-9\.]+)", orca_output, re.IGNORECASE)
    d1_match = re.search(r"D1\s+diagnostic\s*(?:\.+|:)\s*([0-9\.]+)", orca_output, re.IGNORECASE)

    if not t1_match or not d1_match:
        raise ValueError("T1 or D1 diagnostic not found in ORCA output")

    t1 = float(t1_match.group(1))
    d1 = float(d1_match.group(1))

    t1_thresh = 0.045 if spin_multiplicity > 1 else 0.020
    d1_thresh = 0.070 if spin_multiplicity > 1 else 0.050

    passed = (t1 <= t1_thresh) and (d1 <= d1_thresh)
    return {
        "passed": passed,
        "t1": t1,
        "d1": d1,
        "t1_threshold": t1_thresh,
        "d1_threshold": d1_thresh,
        "status": "PASSED" if passed else "MULTIREFERENCE_WARNING",
    }


def audit_diagnostics(
    orca_source: str | Path | OrcaEnergyResult | dict[str, Any],
    *,
    spin_multiplicity: int = 1,
    claimed_error_pct: float | None = None,
    has_core_valence: bool = False,
    s2_expected: float | None = None,
    s2_actual: float | None = None,
) -> dict[str, Any]:
    """Comprehensive diagnostic audit for an ORCA calculation.

    Accepts a file path, raw ORCA output text, ``OrcaEnergyResult``, or dict.
    Handles outputs that report T1 without D1 (e.g., canonical CCSD) without
    fabricating missing diagnostics.
    """
    if isinstance(orca_source, OrcaEnergyResult):
        parsed = orca_source
    elif isinstance(orca_source, dict):
        parsed = OrcaEnergyResult.model_validate(orca_source)
    elif isinstance(orca_source, Path):
        parsed = parse_orca_energy(orca_source)
    elif isinstance(orca_source, str):
        maybe_path = Path(orca_source)
        if "\n" not in orca_source and maybe_path.is_file():
            parsed = parse_orca_energy(maybe_path)
        else:
            parsed = parse_orca_output_text(orca_source)
    else:
        raise TypeError(f"Unsupported orca_source type: {type(orca_source).__name__}")

    t1_thresh = 0.045 if spin_multiplicity > 1 else 0.020
    d1_thresh = 0.070 if spin_multiplicity > 1 else 0.050

    findings: list[str] = list(parsed.warnings)
    passed = parsed.terminated_normally and parsed.scf_converged and (parsed.opt_converged is not False)

    t1_status = "NOT_REPORTED"
    if parsed.t1_diagnostic is not None:
        if parsed.t1_diagnostic <= t1_thresh:
            t1_status = "PASSED"
        else:
            t1_status = "MULTIREFERENCE_WARNING"
            passed = False
            findings.append(
                f"T1 diagnostic {parsed.t1_diagnostic:.6f} exceeds threshold {t1_thresh:.3f}."
            )

    d1_status = "NOT_REPORTED"
    if parsed.d1_diagnostic is not None:
        if parsed.d1_diagnostic <= d1_thresh:
            d1_status = "PASSED"
        else:
            d1_status = "MULTIREFERENCE_WARNING"
            passed = False
            findings.append(
                f"D1 diagnostic {parsed.d1_diagnostic:.6f} exceeds threshold {d1_thresh:.3f}."
            )

    if parsed.slow_conv_detected:
        passed = False

    fc_audit: dict[str, Any] | None = None
    if claimed_error_pct is not None:
        fc_audit = audit_frozen_core_bias(claimed_error_pct, has_core_valence=has_core_valence)
        if not fc_audit["passed"]:
            passed = False
            findings.extend(fc_audit["reasons"])

    spin_audit: dict[str, Any] | None = None
    if s2_expected is not None and s2_actual is not None:
        try:
            check_spin_contamination(s2_expected, s2_actual)
            spin_audit = {"passed": True, "s2_expected": s2_expected, "s2_actual": s2_actual}
        except ValueError as exc:
            passed = False
            findings.append(str(exc))
            spin_audit = {
                "passed": False,
                "s2_expected": s2_expected,
                "s2_actual": s2_actual,
                "error": str(exc),
            }

    return {
        "passed": passed,
        "status": "PASSED" if passed else "AUDIT_WARNING",
        "orca_version": parsed.orca_version,
        "scf_converged": parsed.scf_converged,
        "opt_converged": parsed.opt_converged,
        "terminated_normally": parsed.terminated_normally,
        "slow_conv_detected": parsed.slow_conv_detected,
        "t1_diagnostic": parsed.t1_diagnostic,
        "t1_threshold": t1_thresh,
        "t1_status": t1_status,
        "d1_diagnostic": parsed.d1_diagnostic,
        "d1_threshold": d1_thresh,
        "d1_status": d1_status,
        "frozen_core_audit": fc_audit,
        "spin_audit": spin_audit,
        "findings": findings,
        "provenance_tag": "[M]",
    }


def audit_frozen_core_bias(claimed_error_pct: float, has_core_valence: bool = False) -> dict[str, Any]:
    """Audit accuracy claims <= 0.5% to ensure core-valence basis was evaluated.

    Unaccounted frozen-core correlation introduces a systematic -0.81% bias in rotational constants.
    """
    if claimed_error_pct <= 0.5 and not has_core_valence:
        return {
            "passed": False,
            "status": "FROZEN_CORE_BIAS_WARNING",
            "claimed_error_pct": claimed_error_pct,
            "frozen_core_mean_bias_pct": -0.81,
            "provenance_tag": "[M]",
            "reasons": [
                "Error claim <= 0.5% without core-valence basis incurs mean frozen-core bias "
                "of -0.81% in rotational constants"
            ],
        }
    return {
        "passed": True,
        "status": "PASSED",
        "claimed_error_pct": claimed_error_pct,
        "frozen_core_mean_bias_pct": 0.0 if has_core_valence else -0.81,
        "provenance_tag": "[M]",
        "reasons": [],
    }


def check_spin_contamination(s2_expected: float, s2_actual: float) -> None:
    """Mandate S-squared check for open-shell systems; halt if deviation > 10%."""
    if s2_expected > 0:
        deviation = abs(s2_actual - s2_expected) / s2_expected
        if deviation > 0.10:
            raise ValueError(
                f"Spin contamination exceeds 10%: Expected S^2={s2_expected}, Actual={s2_actual}"
            )


def check_dispersion_correction(method: str, is_weak_complex: bool) -> None:
    """Reject DFT optimizations of weak complexes lacking D3/D4."""
    is_dft = any(d in method.upper() for d in ["B3LYP", "PBE", "M06", "WB97X"])
    has_dispersion = "D3" in method.upper() or "D4" in method.upper()
    if is_weak_complex and is_dft and not has_dispersion:
        raise ValueError(
            "Dispersion correction (D3/D4) is mandatory for DFT optimizations of weak complexes. "
            f"Method: {method}"
        )


def validate_standing_rule_7(payload: dict[str, Any]) -> dict[str, Any]:
    """Enforce Standing Rule 7: Derived [D] and Estimated [E] claims must have explicit [M] anchor."""
    tags = payload.get("provenance_tags", {})
    claims = payload.get("accuracy_claims", [])
    violations: list[str] = []

    for c in claims:
        key = c.get("metric_key", "")
        tag = tags.get(key, "")
        if tag in ["[D]", "[E]"]:
            if not c.get("has_measured_anchor", False):
                violations.append(f"Unanchored {tag} claim for metric '{key}' ({c.get('name', '')})")

    return {
        "rule_7_compliant": len(violations) == 0,
        "violations": violations,
    }


def get_active_orca_version(orca_output_path: str | Path | None = None) -> str:
    """Return active ORCA runtime version stamp from output file or installed binary."""
    if orca_output_path is not None:
        parsed = parse_orca_energy(orca_output_path)
        if parsed.orca_version:
            return f"ORCA {parsed.orca_version}"

    orca_exe = shutil.which("orca")
    if not orca_exe:
        return "ORCA (not installed)"

    # Check installation directory name (e.g. C:\ORCA_6.1.1\orca.EXE)
    path_ver = re.search(r"orca[_-]?([0-9]+\.[0-9]+(?:\.[0-9]+)?)", orca_exe, re.IGNORECASE)
    if path_ver:
        return f"ORCA {path_ver.group(1)}"

    # Fallback: query binary output
    try:
        result = subprocess.run([orca_exe], capture_output=True, text=True, timeout=10)
        output = (result.stdout or "") + "\n" + (result.stderr or "")
        version_match = re.search(
            r"(?:ORCA\s+version|Program\s+Version)\s+([0-9]+\.[0-9]+(?:\.[0-9]+)?)",
            output,
            re.IGNORECASE,
        )
        if version_match:
            return f"ORCA {version_match.group(1)}"
    except OSError:
        pass

    return f"ORCA ({Path(orca_exe).name})"


def verify_openmpi_shm_checksum(data: bytes, expected_checksum: str) -> dict[str, Any]:
    """Compute and verify SHA256 integrity checksum for shared-memory tensor buffers."""
    sha = hashlib.sha256(data).hexdigest()
    passed = sha.lower() == expected_checksum.lower()
    return {
        "passed": passed,
        "checksum": sha,
        "expected_checksum": expected_checksum,
        "bytes_checked": len(data),
    }


def export_provenance_json(
    input_dict: dict[str, Any],
    results_dict: dict[str, Any],
    provenance_tags: dict[str, str] | None = None,
    output_path: Path | None = None,
    accuracy_claims: list[dict[str, Any]] | None = None,
    orca_output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Export full provenance JSON package with Method Matrix v4 tagging discipline."""
    merged_tags = dict(DEFAULT_PROVENANCE_TAGS)
    if provenance_tags:
        merged_tags.update(provenance_tags)

    payload = ProvenancePayload(
        input=input_dict,
        results=results_dict,
        provenance_tags=merged_tags,
        accuracy_claims=accuracy_claims or [],
        orca_version=get_active_orca_version(orca_output_path=orca_output_path),
    )

    payload.rule_7_compliance = validate_standing_rule_7(payload.model_dump())

    if accuracy_claims:
        payload.frozen_core_audit = [
            audit_frozen_core_bias(
                c.get("claimed_error_pct", 1.0),
                c.get("has_core_valence", False),
            )
            for c in accuracy_claims
        ]

    if output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        out_p.write_text(payload.model_dump_json(indent=2), encoding="utf-8")

    return payload.model_dump()
