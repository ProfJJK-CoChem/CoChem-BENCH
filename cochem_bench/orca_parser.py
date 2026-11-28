"""Real ORCA 6.x output parser for energies, convergence, T1/D1 diagnostics, and SlowConv.

Verified directly against converged ORCA 6.1.1 outputs (CCSD/cc-pVTZ and revDSD-PBEP86-D4).
"""

from __future__ import annotations

import re
from pathlib import Path

from pydantic import BaseModel, Field

_FLOAT_RE = r"[-+]?\d+(?:\.\d+)?(?:[eEdD][-+]?\d+)?"

_PROGRAM_VERSION_RE = re.compile(
    r"Program\s+Version\s+([0-9]+\.[0-9]+\.[0-9]+(?:\S*)?)",
    re.IGNORECASE,
)
_FINAL_SPE_RE = re.compile(rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+({_FLOAT_RE})", re.IGNORECASE)
_SCF_TOTAL_ENERGY_RE = re.compile(
    rf"^[ \t]*Total\s+Energy\s*:\s*({_FLOAT_RE})\s+Eh\b",
    re.IGNORECASE | re.MULTILINE,
)
_MP2_TOTAL_ENERGY_RE = re.compile(
    rf"^[ \t]*(?:(?:RI|SCS|OO)-)?MP2\s+TOTAL\s+ENERGY\s*:\s*({_FLOAT_RE})\s+Eh\b",
    re.IGNORECASE | re.MULTILINE,
)
_DFT_DOUBLE_HYBRID_EVIDENCE_RE = re.compile(
    r"^[ \t]*Total\s+energy\s+after\s+final\s+integration\s*:"
    r"|^[ \t]*DFT\s+DISPERSION\s+CORRECTION\b"
    r"|^[ \t]*\|\s*\d+>\s*!.*\b(?:rev)?dsd-[a-z0-9/_-]+"
    r"|^[ \t]*\|\s*\d+>\s*!.*\b(?:b2g?p?-?plyp|mpw2plyp|pwpb95|pbe0-dh|pbe-qidh)"
    r"|^[ \t]*\|\s*\d+>\s*!.*\b(?:wb2g?p?-?plyp|w?pr2scan)",
    re.IGNORECASE | re.MULTILINE,
)
_E0_RE = re.compile(rf"E\(0\)\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_FINAL_CORR_F12_RE = re.compile(
    rf"Final\s+correlation\s+energy\s*\(\s*with\s+F12\s*\)\s*(?:\.+|:)\s*({_FLOAT_RE})",
    re.IGNORECASE,
)
_FINAL_CORR_RE = re.compile(
    rf"Final\s+correlation\s+energy\s*(?:\.+|:)\s*({_FLOAT_RE})",
    re.IGNORECASE,
)
_FINAL_F12_CORR_RE = re.compile(
    rf"Final\s+F12\s+correlation\s+energy\s*(?:\.+|:)\s*({_FLOAT_RE})",
    re.IGNORECASE,
)
_ECORR_CORRECTED_RE = re.compile(
    rf"E\(CORR\)\s*\(\s*corrected\s*\)\s*(?:\.+|:)\s*({_FLOAT_RE})",
    re.IGNORECASE,
)
_ECORR_MDCI_RE = re.compile(rf"E\(CORR\)\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_TRIPLES_RE = re.compile(
    rf"Triples\s+Correction\s*\(T\)\s*(?:\.+|:)\s*({_FLOAT_RE})",
    re.IGNORECASE,
)
_SCALED_TRIPLES_RE = re.compile(
    rf"Scaled\s+triples\s+correction\s*\(T\)\s*(?:\.+|:)\s*({_FLOAT_RE})",
    re.IGNORECASE,
)
_ECORR_MP2_RE = re.compile(
    rf"(?:RI-MP2\s+)?CORRELATION\s+ENERGY\s*:\s*({_FLOAT_RE})\s+Eh",
    re.IGNORECASE,
)
_T1_RE = re.compile(rf"T1\s+diagnostic\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_D1_RE = re.compile(rf"D1\s+diagnostic\s*(?:\.+|:)\s*({_FLOAT_RE})", re.IGNORECASE)
_SCF_CONVERGED_RE = re.compile(r"SCF\s+CONVERGED\s+AFTER\s+(\d+)\s+CYCLES", re.IGNORECASE)
_SCF_NOT_CONVERGED_RE = re.compile(
    r"SCF\s+NOT\s+CONVERGED|Error.*SCF\s+not\s+converged",
    re.IGNORECASE,
)
_OPT_STARTED_RE = re.compile(
    r"ORCA\s+GEOMETRY\s+RELAXATION\s+STEP|GEOMETRY\s+OPTIMIZATION\s+CYCLE",
    re.IGNORECASE,
)
_OPT_CONVERGED_RE = re.compile(r"THE\s+OPTIMIZATION\s+HAS\s+CONVERGED", re.IGNORECASE)
_TERMINATED_NORMALLY_RE = re.compile(r"ORCA\s+TERMINATED\s+NORMALLY", re.IGNORECASE)
_SLOW_CONV_RE = re.compile(
    r"\bSlowConv\b|VerySlowConv|SCF\s+convergence\s+was\s+slow",
    re.IGNORECASE,
)


class OrcaEnergyResult(BaseModel):
    """Structured telemetry parsed from an ORCA output file."""

    source_path: str | None = None
    orca_version: str | None = None
    final_single_point_energy_eh: float | None = None
    scf_energy_eh: float | None = None
    mp2_total_energy_eh: float | None = None
    double_hybrid_total_energy_eh: float | None = None
    correlation_energy_eh: float | None = None
    triples_correction_eh: float | None = None
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


def _parse_float(raw: str) -> float:
    return float(raw.replace("D", "E").replace("d", "e"))


def parse_orca_output_text(text: str, *, source_path: str | None = None) -> OrcaEnergyResult:
    """Parse an ORCA output string into an ``OrcaEnergyResult``."""
    ver_match = _PROGRAM_VERSION_RE.search(text)
    orca_version = ver_match.group(1) if ver_match else None

    spe_matches = _FINAL_SPE_RE.findall(text)
    final_spe = _parse_float(spe_matches[-1]) if spe_matches else None

    # Prefer anchored SCF Total Energy if present, else MDCI E(0) reference.
    # Anchoring prevents post-SCF `MP2 TOTAL ENERGY:` lines (e.g. revDSD-PBEP86-D4)
    # from overwriting `scf_energy_eh`.
    e0_matches = _E0_RE.findall(text)
    scf_matches = _SCF_TOTAL_ENERGY_RE.findall(text)
    if scf_matches:
        scf_energy = _parse_float(scf_matches[-1])
    elif e0_matches:
        scf_energy = _parse_float(e0_matches[-1])
    else:
        scf_energy = None

    mp2_total_matches = _MP2_TOTAL_ENERGY_RE.findall(text)
    mp2_total_energy = _parse_float(mp2_total_matches[-1]) if mp2_total_matches else None

    # Correlation energy priority:
    # 1. Final correlation energy (with F12)  -> CCSD(T)-F12 / DLPNO-CCSD(T)-F12
    # 2. Final correlation energy             -> CCSD(T) / DLPNO-CCSD(T) (includes (T) triples)
    # 3. Final F12 correlation energy         -> CCSD-F12 / DLPNO-CCSD-F12 without triples
    # 4. E(CORR)(corrected)                   -> DLPNO-CCSD without triples
    # 5. E(CORR)                              -> canonical CCSD without triples
    # 6. [RI-MP2] CORRELATION ENERGY          -> MP2 / double-hybrid DFT
    corr_energy: float | None = None
    matched_mp2_corr = False
    for pattern in (
        _FINAL_CORR_F12_RE,
        _FINAL_CORR_RE,
        _FINAL_F12_CORR_RE,
        _ECORR_CORRECTED_RE,
        _ECORR_MDCI_RE,
        _ECORR_MP2_RE,
    ):
        matches = pattern.findall(text)
        if matches:
            corr_energy = _parse_float(matches[-1])
            matched_mp2_corr = pattern is _ECORR_MP2_RE
            break

    is_double_hybrid = (mp2_total_energy is not None or matched_mp2_corr) and bool(
        _DFT_DOUBLE_HYBRID_EVIDENCE_RE.search(text)
    )
    double_hybrid_total_energy = (
        (final_spe if final_spe is not None else mp2_total_energy)
        if is_double_hybrid
        else None
    )

    triples_matches = _SCALED_TRIPLES_RE.findall(text) or _TRIPLES_RE.findall(text)
    triples_corr = _parse_float(triples_matches[-1]) if triples_matches else None

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
        mp2_total_energy_eh=mp2_total_energy,
        double_hybrid_total_energy_eh=double_hybrid_total_energy,
        correlation_energy_eh=corr_energy,
        triples_correction_eh=triples_corr,
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
