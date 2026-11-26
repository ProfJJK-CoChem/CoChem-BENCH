"""Authentic ORCA 6.x execution runner for CoChem-BENCH.

Resolves a real local ORCA 6 executable (ORCA_BIN / COCHEM_ORCA_BIN / PATH /
~/CoChem_Workstation_Direct), runs it without mocks or fabricated output, enforces
rc == 0 + normal termination + real Program Version 6.x banner, parses energies and
T1/D1 diagnostics via ``orca_parser``, and feeds directly into ``verifier.py`` with
SHA-256 provenance. Fails closed with ``FileNotFoundError`` when ORCA is absent.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from .orca_parser import OrcaEnergyResult, parse_orca_energy
from .orca_writer import generate_dlpno_ccsd_f12_with_metadata
from .verifier import audit_diagnostics, export_provenance_json

_ORCA_NORMAL_TERMINATION_MARKER = "****ORCA TERMINATED NORMALLY****"
_ORCA_VERSION_6_BANNER_RE = re.compile(r"Program\s+Version\s+(6\.[0-9]+(?:\.[0-9]+)?(?:\S*)?)", re.IGNORECASE)
_ENV_ORCA_KEYS: tuple[str, ...] = ("ORCA_BIN", "COCHEM_ORCA_BIN", "ORCA_CMD", "TOPOS_ORCA_EXECUTABLE")


class OrcaExecutionError(RuntimeError):
    """Raised when an ORCA 6 calculation fails execution or strict output verification."""


class OrcaRunResult(BaseModel):
    """Structured authentic execution receipt and verified telemetry from ORCA 6."""

    job_name: str
    workdir: str
    input_path: str
    output_path: str
    stderr_path: str
    orca_binary: str
    orca_binary_sha256: str
    input_sha256: str
    output_sha256: str
    returncode: int
    elapsed_seconds: float
    orca_version: str
    final_single_point_energy_eh: float
    scf_energy_eh: float | None = None
    correlation_energy_eh: float | None = None
    t1_diagnostic: float | None = None
    d1_diagnostic: float | None = None
    parsed: OrcaEnergyResult
    audit: dict[str, Any]
    provenance: dict[str, Any]
    provenance_tags: dict[str, str] = Field(default_factory=dict)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_orca_path_from_activate_script(script_path: Path) -> Path | None:
    """Extract COCHEM_ORCA_BIN / ORCA_DIR from activate_cochem_wsl.sh without executing it."""
    if not script_path.is_file():
        return None
    try:
        text = script_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    orca_dir_match = re.search(r'(?m)^\s*export\s+ORCA_DIR="([^"]+)"', text)
    if orca_dir_match:
        candidate = Path(orca_dir_match.group(1)) / "orca"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    orca_bin_match = re.search(r'(?m)^\s*export\s+COCHEM_ORCA_BIN="([^"$]+)"', text)
    if orca_bin_match:
        candidate = Path(orca_bin_match.group(1))
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    return None


def _discover_workstation_orca_candidates() -> list[Path]:
    """Return known CoChem workstation ORCA 6 installation candidates."""
    candidates: list[Path] = []
    home = Path.home()

    orca_dir_env = os.environ.get("ORCA_DIR")
    if orca_dir_env:
        candidates.append(Path(orca_dir_env) / "orca")
        candidates.append(Path(orca_dir_env) / "orca.exe")

    ws_root = home / "CoChem_Workstation_Direct"
    activate_candidate = _parse_orca_path_from_activate_script(ws_root / "activate_cochem_wsl.sh")
    if activate_candidate is not None:
        candidates.append(activate_candidate)

    for pattern in (
        ws_root / "engines",
        home / ".local" / "share" / "cochem-topos",
    ):
        if pattern.is_dir():
            for found in sorted(pattern.glob("**/orca_6_*/orca"), reverse=True):
                candidates.append(found)
            for found in sorted(pattern.glob("**/orca"), reverse=True):
                if found.is_file():
                    candidates.append(found)

    return candidates


def resolve_orca_binary(orca_bin: str | Path | None = None) -> Path:
    """Resolve a real ORCA executable path or fail closed with ``FileNotFoundError``."""
    if orca_bin is not None:
        requested = Path(orca_bin).expanduser()
        if requested.is_file():
            return requested.resolve()
        found = shutil.which(str(orca_bin))
        if found:
            return Path(found).resolve()
        raise FileNotFoundError(f"Explicitly requested ORCA binary not found: {orca_bin}")

    for env_key in _ENV_ORCA_KEYS:
        env_val = os.environ.get(env_key, "").strip()
        if not env_val:
            continue
        candidate = Path(env_val).expanduser()
        if candidate.is_file():
            return candidate.resolve()
        found = shutil.which(env_val)
        if found:
            return Path(found).resolve()
        if env_key == "ORCA_BIN":
            raise FileNotFoundError(f"ORCA_BIN points to a missing executable: {env_val}")

    on_path = shutil.which("orca")
    if on_path:
        return Path(on_path).resolve()

    for candidate in _discover_workstation_orca_candidates():
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()

    raise FileNotFoundError(
        "ORCA 6 executable not found via ORCA_BIN, PATH, or ~/CoChem_Workstation_Direct. "
        "Install ORCA 6 or set ORCA_BIN; simulated or mocked execution is forbidden."
    )


def _build_orca_env(orca_path: Path, deck_text: str) -> dict[str, str]:
    """Prepare subprocess environment for ORCA 6 including OpenMPI shared libraries."""
    env = os.environ.copy()
    orca_dir = str(orca_path.parent)

    path_parts = [orca_dir]
    if env.get("PATH"):
        path_parts.append(env["PATH"])
    env["PATH"] = os.pathsep.join(path_parts)

    if os.name == "posix":
        ld_parts: list[str] = []
        cochem_ld = env.get("COCHEM_ORCA_LD_LIBRARY_PATH", "").strip()
        if cochem_ld:
            ld_parts.extend(p for p in cochem_ld.split(":") if p)
        # Auto-discover sibling OpenMPI runtime lib next to orca-native-* directories.
        engines_root = orca_path.parent.parent.parent
        if engines_root.is_dir():
            for mpi_lib in sorted(engines_root.glob("openmpi*/lib")):
                if mpi_lib.is_dir():
                    ld_parts.append(str(mpi_lib))
        ld_parts.append(orca_dir)
        if env.get("LD_LIBRARY_PATH"):
            ld_parts.extend(p for p in env["LD_LIBRARY_PATH"].split(":") if p)
        deduped_ld = list(dict.fromkeys(ld_parts))
        env["LD_LIBRARY_PATH"] = ":".join(deduped_ld)

    # Guard against OpenMP oversubscription when %pal nprocs > 1.
    pal_match = re.search(r"(?i)%pal\s+nprocs\s+(\d+)\s+end", deck_text)
    if pal_match and int(pal_match.group(1)) > 1:
        env["OMP_NUM_THREADS"] = "1"
    else:
        env.setdefault("OMP_NUM_THREADS", "2")
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    return env


def run_orca_input(
    deck_text: str,
    workdir: str | Path,
    *,
    job_name: str = "bench_job",
    orca_bin: str | Path | None = None,
    timeout_seconds: float = 600.0,
    spin_multiplicity: int = 1,
    claimed_error_pct: float | None = None,
    has_core_valence: bool = False,
    require_opt_converged: bool | None = None,
    provenance_output_path: str | Path | None = None,
    input_metadata: dict[str, Any] | None = None,
) -> OrcaRunResult:
    """Execute a real ORCA 6 input deck, verify output strictly, and return SHA-256 provenance."""
    if not deck_text or not deck_text.strip():
        raise ValueError("ORCA input deck must be a non-empty string.")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", job_name):
        raise ValueError(f"Invalid job_name {job_name!r}; must match [A-Za-z0-9._-]+.")

    binary = resolve_orca_binary(orca_bin)
    binary_sha256 = _sha256_file(binary)

    work_path = Path(workdir).resolve()
    work_path.mkdir(parents=True, exist_ok=True)
    inp_path = work_path / f"{job_name}.inp"
    out_path = work_path / f"{job_name}.out"
    err_path = work_path / f"{job_name}.err"

    inp_path.write_text(deck_text, encoding="utf-8")
    input_sha256 = _sha256_file(inp_path)

    env = _build_orca_env(binary, deck_text)
    run_kwargs: dict[str, Any] = {
        "cwd": str(work_path),
        "env": env,
        "capture_output": True,
        "text": True,
        "timeout": timeout_seconds,
    }
    if os.name == "nt":
        # Workstation rule: avoid popping console windows for background subprocesses on Windows.
        run_kwargs["creationflags"] = 0x08000000

    t0 = time.monotonic()
    completed = subprocess.run([str(binary), inp_path.name], **run_kwargs)
    elapsed = time.monotonic() - t0

    stdout_text = completed.stdout or ""
    stderr_text = completed.stderr or ""
    out_path.write_text(stdout_text, encoding="utf-8")
    err_path.write_text(stderr_text, encoding="utf-8")
    output_sha256 = _sha256_file(out_path)

    if completed.returncode != 0:
        tail = "\n".join(stdout_text.splitlines()[-20:])
        raise OrcaExecutionError(
            f"ORCA exited with rc={completed.returncode} for job {job_name!r}.\n"
            f"STDERR: {stderr_text.strip()[:500]}\nSTDOUT tail:\n{tail}"
        )

    if _ORCA_NORMAL_TERMINATION_MARKER not in stdout_text:
        raise OrcaExecutionError(
            f"ORCA exited rc=0 but lacked {_ORCA_NORMAL_TERMINATION_MARKER!r} in {out_path}."
        )

    banner_match = _ORCA_VERSION_6_BANNER_RE.search(stdout_text)
    if banner_match is None:
        raise OrcaExecutionError(
            f"ORCA output {out_path} lacks an authentic 'Program Version 6.' banner."
        )

    parsed = parse_orca_energy(out_path)
    if parsed.orca_version is None or not parsed.orca_version.startswith("6."):
        raise OrcaExecutionError(
            f"Parsed ORCA version {parsed.orca_version!r} is not ORCA 6.x in {out_path}."
        )
    if parsed.final_single_point_energy_eh is None or not math.isfinite(parsed.final_single_point_energy_eh):
        raise OrcaExecutionError(
            f"No finite FINAL SINGLE POINT ENERGY found in ORCA output {out_path}."
        )

    is_opt_expected = require_opt_converged
    if is_opt_expected is None:
        is_opt_expected = parsed.opt_converged is not None
    if is_opt_expected and parsed.opt_converged is not True:
        raise OrcaExecutionError(
            f"ORCA geometry optimization did not converge in {out_path} (opt_converged={parsed.opt_converged!r})."
        )

    audit = audit_diagnostics(
        out_path,
        spin_multiplicity=spin_multiplicity,
        claimed_error_pct=claimed_error_pct,
        has_core_valence=has_core_valence,
    )

    provenance_tags: dict[str, str] = {
        "final_single_point_energy_eh": "[M]",
        "scf_energy_eh": "[M]",
        "correlation_energy_eh": "[M]",
        "t1_diagnostic": "[M]",
        "d1_diagnostic": "[M]",
    }
    prov_input: dict[str, Any] = {
        "job_name": job_name,
        "input_path": str(inp_path),
        "input_sha256": input_sha256,
        "orca_binary": str(binary),
        "orca_binary_sha256": binary_sha256,
        **(input_metadata or {}),
    }
    prov_results: dict[str, Any] = {
        "output_path": str(out_path),
        "output_sha256": output_sha256,
        "returncode": completed.returncode,
        "elapsed_seconds": elapsed,
        "final_single_point_energy_eh": parsed.final_single_point_energy_eh,
        "scf_energy_eh": parsed.scf_energy_eh,
        "correlation_energy_eh": parsed.correlation_energy_eh,
        "t1_diagnostic": parsed.t1_diagnostic,
        "d1_diagnostic": parsed.d1_diagnostic,
        "audit_passed": audit["passed"],
        "audit_status": audit["status"],
    }
    prov_out_path = Path(provenance_output_path) if provenance_output_path is not None else (work_path / f"{job_name}_provenance.json")
    provenance = export_provenance_json(
        input_dict=prov_input,
        results_dict=prov_results,
        provenance_tags=provenance_tags,
        output_path=prov_out_path,
        orca_output_path=out_path,
    )

    return OrcaRunResult(
        job_name=job_name,
        workdir=str(work_path),
        input_path=str(inp_path),
        output_path=str(out_path),
        stderr_path=str(err_path),
        orca_binary=str(binary),
        orca_binary_sha256=binary_sha256,
        input_sha256=input_sha256,
        output_sha256=output_sha256,
        returncode=completed.returncode,
        elapsed_seconds=elapsed,
        orca_version=parsed.orca_version,
        final_single_point_energy_eh=parsed.final_single_point_energy_eh,
        scf_energy_eh=parsed.scf_energy_eh,
        correlation_energy_eh=parsed.correlation_energy_eh,
        t1_diagnostic=parsed.t1_diagnostic,
        d1_diagnostic=parsed.d1_diagnostic,
        parsed=parsed,
        audit=audit,
        provenance=provenance,
        provenance_tags=provenance_tags,
    )


def run_orca_calculation(
    coords: list[tuple[str, float, float, float] | list[Any]],
    workdir: str | Path,
    *,
    basis: str = "cc-pVDZ-F12",
    method: str = "DLPNO-CCSD(T)-F12",
    is_opt: bool = False,
    in_hess: str = "XTB2",
    rel_mode: str = "auto",
    tight_scf: bool = False,
    charge: int = 0,
    mult: int = 1,
    nprocs: int = 1,
    pno_setting: str = "TightPNO",
    maxcore: int | None = None,
    job_name: str = "bench_job",
    orca_bin: str | Path | None = None,
    timeout_seconds: float = 600.0,
) -> OrcaRunResult:
    """Generate a valid ORCA 6.1.1 deck and execute it live through ``run_orca_input``."""
    deck, basis_meta = generate_dlpno_ccsd_f12_with_metadata(
        coords=coords,
        basis=basis,
        is_opt=is_opt,
        in_hess=in_hess,
        rel_mode=rel_mode,
        tight_scf=tight_scf,
        charge=charge,
        mult=mult,
        nprocs=nprocs,
        method=method,
        pno_setting=pno_setting,
        maxcore=maxcore,
    )
    return run_orca_input(
        deck_text=deck,
        workdir=workdir,
        job_name=job_name,
        orca_bin=orca_bin,
        timeout_seconds=timeout_seconds,
        spin_multiplicity=mult,
        require_opt_converged=is_opt,
        input_metadata={
            "method": method,
            "basis": basis,
            "effective_basis": basis_meta["effective_basis"],
            "basis_swapped": basis_meta["basis_swapped"],
            "basis_warnings": basis_meta["warnings"],
            "is_opt": is_opt,
            "rel_mode": rel_mode,
            "charge": charge,
            "mult": mult,
            "nprocs": nprocs,
        },
    )
