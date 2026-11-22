"""Real ORCA 6.x output parser for energies, convergence, T1/D1 diagnostics, and SlowConv.

Verified directly against converged ORCA 6.1.1 outputs (CCSD/cc-pVTZ and revDSD-PBEP86-D4).
"""

from __future__ import annotations

import re
from pathlib import Path

from pydantic import BaseModel, Field

_FLOAT_RE = r"[-+]?\d+(?:\.\d+)?(?:[eEdD][-+]?\d+)?"

_PROGRAM_VERSION_RE = re.compile(r"Program\s+Version\s+([0-9]+\.[0-9]+\.[0-9]+(?:\S*)?)", re.IGNORECASE)
_FINAL_SPE_RE = re.compile(rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+({_FLOAT_RE})", re.IGNORECASE)
_SCF_TOTAL_ENERGY_RE = re.compile(rf"Total\s+Energy\s*:\s*({_FLOAT_RE})\s+Eh", re.IGNORECASE)
_E0_RE = re.compile(rf"E\(0\)\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_ECORR_MDCI_RE = re.compile(rf"E\(CORR\)\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_ECORR_MP2_RE = re.compile(
    rf"(?:RI-MP2\s+)?CORRELATION\s+ENERGY\s*:\s*({_FLOAT_RE})\s+Eh",
    re.IGNORECASE,
)
_T1_RE = re.compile(rf"T1\s+diagnostic\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_D1_RE = re.compile(rf"D1\s+diagnostic\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_SCF_CONVERGED_RE = re.compile(r"SCF\s+CONVERGED\s+AFTER\s+(\d+)\s+CYCLES", re.IGNORECASE)
_SCF_NOT_CONVERGED_RE = re.compile(r"SCF\s+NOT\s+CONVERGED|Error.*SCF\s+not\s+converged", re.IGNORECASE)
_OPT_STARTED_RE = re.compile(r"ORCA\s+GEOMETRY\s+RELAXATION\s+STEP|GEOMETRY\s+OPTIMIZATION\s+CYCLE", re.IGNORECASE)
_OPT_CONVERGED_RE = re.compile(r"THE\s+OPTIMIZATION\s+HAS\s+CONVERGED", re.IGNORECASE)
_TERMINATED_NORMALLY_RE = re.compile(r"ORCA\s+TERMINATED\s+NORMALLY", re.IGNORECASE)
_SLOW_CONV_RE = re.compile(r"\bSlowConv\b|VerySlowConv|SCF\s+convergence\s+was\s+slow", re.IGNORECASE)


class OrcaEnergyResult(BaseModel):
    """Structured telemetry parsed from an ORCA output file."""

    source_path: str | None = None
    orca_version: str | None = None
    final_single_point_energy_eh: float | None = None
    scf_energy_eh: float | None = None
    correlation_energy_eh: float | None = None
    t1_diagnostic: float | None = None
    d1_diagnostic: float | None = None
    scf_converged: bool = False
    scf_cycles: int | None = None
    opt_converged: bool | None = None
    terminated_normally: bool = False
    slow_conv_detected: bool = False
    n_single_points: int = 0
    warnings: list[str] = Field(default_factory=list)
    provenance_tag: str = "[M]"


def parse_orca_output_text(text: str, *, source_path: str | None = None) -> OrcaEnergyResult:
    """Parse an ORCA output string into an ``OrcaEnergyResult``."""
    ver_match = _PROGRAM_VERSION_RE.search(text)
    orca_version = ver_match.group(1) if ver_match else None

    spe_matches = _FINAL_SPE_RE.findall(text)
    final_spe = float(spe_matches[-1].replace("D", "E").replace("d", "e")) if spe_matches else None

    # Prefer MDCI E(0) if present (exact Hartree-Fock reference for CC), else SCF Total Energy
    e0_matches = _E0_RE.findall(text)
    scf_matches = _SCF_TOTAL_ENERGY_RE.findall(text)
    if scf_matches:
        scf_energy = float(scf_matches[-1].replace("D", "E").replace("d", "e"))
    elif e0_matches:
        scf_energy = float(e0_matches[-1].replace("D", "E").replace("d", "e"))
    else:
        scf_energy = None

    ecorr_mdci = _ECORR_MDCI_RE.findall(text)
    ecorr_mp2 = _ECORR_MP2_RE.findall(text)
    if ecorr_mdci:
        corr_energy = float(ecorr_mdci[-1].replace("D", "E").replace("d", "e"))
    elif ecorr_mp2:
        corr_energy = float(ecorr_mp2[-1].replace("D", "E").replace("d", "e"))
    else:
        corr_energy = None

    t1_matches = _T1_RE.findall(text)
    t1_val = float(t1_matches[-1]) if t1_matches else None

    d1_matches = _D1_RE.findall(text)
    d1_val = float(d1_matches[-1]) if d1_matches else None

    scf_conv_matches = _SCF_CONVERGED_RE.findall(text)
    scf_not_conv = bool(_SCF_NOT_CONVERGED_RE.search(text))
    scf_converged = bool(scf_conv_matches) and not scf_not_conv
    scf_cycles = int(scf_conv_matches[-1]) if scf_conv_matches else None

    opt_ran = bool(_OPT_STARTED_RE.search(text)) or bool(_OPT_CONVERGED_RE.search(text))
    opt_converged: bool | None = None
    if opt_ran:
        opt_converged = bool(_OPT_CONVERGED_RE.search(text))

    terminated_normally = bool(_TERMINATED_NORMALLY_RE.search(text))
    slow_conv = bool(_SLOW_CONV_RE.search(text))

    warnings: list[str] = []
    if slow_conv:
        warnings.append("SlowConv keyword or slow SCF convergence detected.")
    if scf_not_conv:
        warnings.append("SCF convergence failure detected.")
    if not terminated_normally:
        warnings.append("ORCA did not report normal termination.")
    if opt_converged is False:
        warnings.append("Geometry optimization did not converge.")

    return OrcaEnergyResult(
        source_path=source_path,
        orca_version=orca_version,
        final_single_point_energy_eh=final_spe,
        scf_energy_eh=scf_energy,
        correlation_energy_eh=corr_energy,
        t1_diagnostic=t1_val,
        d1_diagnostic=d1_val,
        scf_converged=scf_converged,
        scf_cycles=scf_cycles,
        opt_converged=opt_converged,
        terminated_normally=terminated_normally,
        slow_conv_detected=slow_conv,
        n_single_points=len(spe_matches),
        warnings=warnings,
    )


def parse_orca_energy(path: str | Path) -> OrcaEnergyResult:
    """Read and parse a real ORCA output file from disk."""
    out_path = Path(path)
    if not out_path.is_file():
        raise FileNotFoundError(f"ORCA output file not found: {out_path}")
    text = out_path.read_text(encoding="utf-8", errors="replace")
    return parse_orca_output_text(text, source_path=str(out_path))
