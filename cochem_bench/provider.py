"""CoChem-BASE `cochem.modules` entry point provider (`bench_provider`).

Implements the versioned `cochem.module-provider/1` and `cochem.module-handoff/1`
contracts for CoChem-BASE integration.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .constants import HARTREE_TO_CM_INV, HARTREE_TO_KCAL_MOL, HARTREE_TO_KJ_MOL
from .extrapolate import cbs_extrapolate
from .orca_parser import parse_orca_energy, parse_orca_output_text
from .thermochemistry import compute_relative_energies
from .verifier import audit_diagnostics

SCIENTIFIC_API_VERSION = "cochem.module-provider/1"
INTEGRATION_CONTRACT = "cochem.module-handoff/1"
OPERATIONS: tuple[str, ...] = (
    "cbs_extrapolate",
    "relative_energies",
    "parse_orca_energy",
    "audit_diagnostics",
)


def metadata() -> dict[str, Any]:
    """Return side-effect-free CoChem-BASE module capability metadata."""
    return {
        "scientific_api_version": SCIENTIFIC_API_VERSION,
        "module_id": "bench",
        "distribution": "cochem-bench",
        "version": "0.2.0",
        "integration_contract": INTEGRATION_CONTRACT,
        "operations": list(OPERATIONS),
        "execute_handoff": "cochem_bench.provider:execute_handoff",
        "supported_basis_families": ["cc-pVnZ", "aug-cc-pVnZ", "def2", "ano-pVnZ"],
        "codata_constants": {
            "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
            "hartree_to_kj_mol": HARTREE_TO_KJ_MOL,
            "hartree_to_cm_inv": HARTREE_TO_CM_INV,
        },
    }


def dispatch_operation(operation: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Execute a declared canonical CoChem-BENCH operation."""
    if operation not in OPERATIONS:
        raise ValueError(
            f"Unsupported CoChem-BENCH operation {operation!r}; valid: {OPERATIONS}"
        )

    if operation == "cbs_extrapolate":
        res = cbs_extrapolate(
            e_scf_x=float(payload["e_scf_x"]),
            e_scf_y=float(payload["e_scf_y"]),
            e_corr_x=payload.get("e_corr_x"),
            e_corr_y=payload.get("e_corr_y"),
            e_tot_x=payload.get("e_tot_x"),
            e_tot_y=payload.get("e_tot_y"),
            x=int(payload.get("x", 3)),
            y=int(payload.get("y", 4)),
            basis_family=str(payload.get("basis_family", "cc-pVnZ")),
            scheme=payload.get("scheme", "helgaker"),
            alpha=payload.get("alpha"),
            beta=payload.get("beta"),
        )
        return res.model_dump(mode="json")

    if operation == "relative_energies":
        table = compute_relative_energies(
            species=payload["species"],
            reference=payload.get("reference"),
            required_corrections=payload.get("required_corrections", ()),
            energy_field=payload.get("energy_field", "e_composite_eh"),
        )
        return table.model_dump(mode="json")

    if operation == "parse_orca_energy":
        if "path" in payload and payload["path"] is not None:
            return parse_orca_energy(Path(payload["path"])).model_dump(mode="json")
        if "text" in payload and payload["text"] is not None:
            return parse_orca_output_text(str(payload["text"])).model_dump(mode="json")
        raise ValueError("parse_orca_energy requires either 'path' or 'text' in payload.")

    if operation == "audit_diagnostics":
        source = payload.get("path") or payload.get("text") or payload.get("parsed")
        if source is None:
            raise ValueError(
                "audit_diagnostics requires 'path', 'text', or 'parsed' in payload."
            )
        return audit_diagnostics(
            source,
            spin_multiplicity=int(payload.get("spin_multiplicity", 1)),
            claimed_error_pct=payload.get("claimed_error_pct"),
            has_core_valence=bool(payload.get("has_core_valence", False)),
            s2_expected=payload.get("s2_expected"),
            s2_actual=payload.get("s2_actual"),
        )

    raise ValueError(f"Unhandled operation: {operation!r}")


def execute_handoff(
    handoff_path: str | Path,
    output_directory: str | Path,
    registry: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Execute a validated CoChem-BASE handoff JSON envelope and write structured output."""
    handoff_file = Path(handoff_path)
    out_dir = Path(output_directory)
    out_dir.mkdir(parents=True, exist_ok=True)

    envelope = json.loads(handoff_file.read_text(encoding="utf-8"))
    operation = str(envelope.get("operation", "")).strip()
    params = envelope.get("parameters", {}) or {}
    result = dispatch_operation(operation, params)

    out_file = out_dir / f"{operation}_result.json"
    out_file.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return {
        "module_id": "bench",
        "operation": operation,
        "output_path": str(out_file),
        "result": result,
    }


class BenchProvider:
    """CoChem-BASE `cochem.modules` provider object for CoChem-BENCH."""

    module_id = "bench"
    distribution = "cochem-bench"
    scientific_api_version = SCIENTIFIC_API_VERSION
    integration_contract = INTEGRATION_CONTRACT
    operations = OPERATIONS
    metadata = staticmethod(metadata)
    dispatch_operation = staticmethod(dispatch_operation)
    execute_handoff = staticmethod(execute_handoff)

    def __call__(self) -> dict[str, Any]:
        return metadata()


bench_provider = BenchProvider()
