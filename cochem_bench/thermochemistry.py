"""Additive composite thermochemistry and relative energy ranking tables.

Supports:
- Core-valence correlation correction (Delta E_CV = E_AE - E_FC)
- Scalar relativistic correction (Delta E_DKH2 / Delta E_X2C = E_rel - E_nonrel)
- Zero-point vibrational energy (Delta ZPE) and thermal H/G corrections
- Relative energy ranking in Eh, kcal/mol, kJ/mol, and cm^-1 using exact CODATA constants.

Strict Integrity Policy:
Missing corrections are NEVER silently substituted with 0.0; they are explicitly
recorded as ``PENDING`` in ``correction_status`` and ``pending_corrections``.
"""

from __future__ import annotations

import math
from typing import Any, Literal

from pydantic import BaseModel, Field

from .constants import HARTREE_TO_CM_INV, HARTREE_TO_KCAL_MOL, HARTREE_TO_KJ_MOL

CorrectionState = Literal["PRESENT", "PENDING"]


class CompositeEnergyResult(BaseModel):
    """Additive composite energy breakdown with explicit PENDING tracking."""

    name: str
    e_cbs_eh: float
    base_energy_key: str = "e_cbs_eh"
    delta_cv_eh: float | None = None
    delta_rel_eh: float | None = None
    rel_method: str | None = None
    delta_zpe_eh: float | None = None
    delta_thermal_h_eh: float | None = None
    delta_thermal_g_eh: float | None = None
    correction_status: dict[str, CorrectionState] = Field(default_factory=dict)
    pending_corrections: list[str] = Field(default_factory=list)
    required_corrections: list[str] = Field(default_factory=list)
    status: Literal["COMPLETE", "PENDING"] = "COMPLETE"
    e_composite_eh: float | None = None
    h_composite_eh: float | None = None
    g_composite_eh: float | None = None


class RelativeEnergyEntry(BaseModel):
    """Single row in a relative energy ranking table."""

    rank: int | None
    name: str
    status: Literal["COMPLETE", "PENDING"]
    pending_corrections: list[str] = Field(default_factory=list)
    e_total_eh: float | None
    delta_e_eh: float | None
    delta_e_kcal_mol: float | None
    delta_e_kj_mol: float | None
    delta_e_cm_inv: float | None
    correction_status: dict[str, CorrectionState] = Field(default_factory=dict)


class RelativeEnergyTable(BaseModel):
    """Ranked relative energy table across conformers or isomers."""

    reference_species: str | None
    reference_energy_eh: float | None
    energy_field: str
    base_energy_key: str | None = None
    protocol_signature: dict[str, str] = Field(default_factory=dict)
    entries: list[RelativeEnergyEntry]
    codata_constants: dict[str, float] = Field(
        default_factory=lambda: {
            "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
            "hartree_to_kj_mol": HARTREE_TO_KJ_MOL,
            "hartree_to_cm_inv": HARTREE_TO_CM_INV,
        }
    )


def compute_composite_energy(
    *,
    name: str,
    e_cbs_eh: float,
    base_energy_key: str = "e_cbs_eh",
    delta_cv_eh: float | None = None,
    e_ae_cv_eh: float | None = None,
    e_fc_cv_eh: float | None = None,
    delta_rel_eh: float | None = None,
    e_rel_eh: float | None = None,
    e_nonrel_eh: float | None = None,
    rel_method: str | None = None,
    delta_zpe_eh: float | None = None,
    delta_thermal_h_eh: float | None = None,
    delta_thermal_g_eh: float | None = None,
    required_corrections: tuple[str, ...] | list[str] = (),
) -> CompositeEnergyResult:
    """Compute additive composite energy with explicit PENDING status for missing terms.

    Core-valence correction can be supplied either directly as ``delta_cv_eh`` or
    via ``(e_ae_cv_eh, e_fc_cv_eh)`` where ``delta_cv_eh = e_ae_cv_eh - e_fc_cv_eh``.
    Scalar relativistic correction (DKH2/X2C) can be supplied directly as
    ``delta_rel_eh`` or via ``(e_rel_eh, e_nonrel_eh)``.
    """
    resolved_cv = delta_cv_eh
    if resolved_cv is None and e_ae_cv_eh is not None and e_fc_cv_eh is not None:
        resolved_cv = float(e_ae_cv_eh) - float(e_fc_cv_eh)

    resolved_rel = delta_rel_eh
    if resolved_rel is None and e_rel_eh is not None and e_nonrel_eh is not None:
        resolved_rel = float(e_rel_eh) - float(e_nonrel_eh)

    terms: dict[str, float | None] = {
        "delta_cv": float(resolved_cv) if resolved_cv is not None else None,
        "delta_rel": float(resolved_rel) if resolved_rel is not None else None,
        "delta_zpe": float(delta_zpe_eh) if delta_zpe_eh is not None else None,
        "delta_thermal_h": float(delta_thermal_h_eh) if delta_thermal_h_eh is not None else None,
        "delta_thermal_g": float(delta_thermal_g_eh) if delta_thermal_g_eh is not None else None,
    }

    correction_status: dict[str, CorrectionState] = {
        key: ("PRESENT" if val is not None else "PENDING")
        for key, val in terms.items()
    }

    req_list = list(required_corrections)
    pending_required = [k for k in req_list if correction_status.get(k, "PENDING") == "PENDING"]
    all_pending = [k for k, st in correction_status.items() if st == "PENDING"]

    if pending_required:
        return CompositeEnergyResult(
            name=name,
            e_cbs_eh=float(e_cbs_eh),
            base_energy_key=base_energy_key,
            delta_cv_eh=terms["delta_cv"],
            delta_rel_eh=terms["delta_rel"],
            rel_method=rel_method,
            delta_zpe_eh=terms["delta_zpe"],
            delta_thermal_h_eh=terms["delta_thermal_h"],
            delta_thermal_g_eh=terms["delta_thermal_g"],
            correction_status=correction_status,
            pending_corrections=pending_required,
            required_corrections=req_list,
            status="PENDING",
            e_composite_eh=None,
            h_composite_eh=None,
            g_composite_eh=None,
        )

    e_comp = float(e_cbs_eh)
    if terms["delta_cv"] is not None:
        e_comp += terms["delta_cv"]
    if terms["delta_rel"] is not None:
        e_comp += terms["delta_rel"]
    if terms["delta_zpe"] is not None:
        e_comp += terms["delta_zpe"]

    h_comp = (
        e_comp + terms["delta_thermal_h"]
        if terms["delta_thermal_h"] is not None
        else None
    )
    g_comp = (
        e_comp + terms["delta_thermal_g"]
        if terms["delta_thermal_g"] is not None
        else None
    )

    return CompositeEnergyResult(
        name=name,
        e_cbs_eh=float(e_cbs_eh),
        base_energy_key=base_energy_key,
        delta_cv_eh=terms["delta_cv"],
        delta_rel_eh=terms["delta_rel"],
        rel_method=rel_method,
        delta_zpe_eh=terms["delta_zpe"],
        delta_thermal_h_eh=terms["delta_thermal_h"],
        delta_thermal_g_eh=terms["delta_thermal_g"],
        correction_status=correction_status,
        pending_corrections=all_pending,
        required_corrections=req_list,
        status="COMPLETE",
        e_composite_eh=e_comp,
        h_composite_eh=h_comp,
        g_composite_eh=g_comp,
    )


def compute_relative_energies(
    species: list[dict[str, Any] | CompositeEnergyResult],
    *,
    reference: str | None = None,
    required_corrections: tuple[str, ...] | list[str] = (),
    energy_field: Literal["e_composite_eh", "h_composite_eh", "g_composite_eh"] = "e_composite_eh",
) -> RelativeEnergyTable:
    """Build a relative energy ranking table in Eh, kcal/mol, kJ/mol, and cm^-1.

    Species missing required corrections remain explicitly ``PENDING`` with ``None``
    relative energies (never coerced to 0.0).
    """
    if not species:
        raise ValueError("species list must not be empty.")

    evaluated: list[CompositeEnergyResult] = []
    uniform_energy_key: str | None = None
    protocol_signatures: list[tuple[str, dict[str, str]]] = []
    _PROTOCOL_KEYS = ("method", "basis", "level_of_theory", "functional")
    for item in species:
        if isinstance(item, CompositeEnergyResult):
            item_key = item.base_energy_key
            if uniform_energy_key is not None and uniform_energy_key != item_key:
                raise ValueError(
                    f"Inconsistent base energy keys across species ({uniform_energy_key!r} vs {item_key!r}); "
                    "relative energies require a uniform energy definition across all species."
                )
            uniform_energy_key = item_key
            evaluated.append(item)
        else:
            item_req = item.get("required_corrections", required_corrections)
            name_val = str(item.get("name") or item.get("source_path") or "")
            if not name_val:
                raise ValueError("Each species entry must include a nonempty 'name' field.")
            matched = next(
                (
                    (key, item[key])
                    for key in ("e_cbs_eh", "final_single_point_energy_eh", "e_total_eh", "energy_eh")
                    if key in item and item[key] is not None
                ),
                None,
            )
            if matched is None:
                raise ValueError(
                    f"Species {name_val!r} is missing a finite base energy "
                    "(expected 'e_cbs_eh' or 'final_single_point_energy_eh')."
                )
            matched_key, raw_energy = matched
            energy_float = float(raw_energy)
            if not math.isfinite(energy_float):
                raise ValueError(f"Species {name_val!r} has non-finite base energy ({raw_energy!r}).")
            if uniform_energy_key is None:
                uniform_energy_key = matched_key
            elif matched_key != uniform_energy_key:
                raise ValueError(
                    f"Inconsistent base energy keys across species ({uniform_energy_key!r} vs {matched_key!r}); "
                    "relative energies require a uniform energy definition across all species."
                )
            sig = {
                k: str(item[k]).strip()
                for k in _PROTOCOL_KEYS
                if item.get(k) is not None and str(item[k]).strip()
            }
            protocol_signatures.append((name_val, sig))
            evaluated.append(
                compute_composite_energy(
                    name=name_val,
                    e_cbs_eh=energy_float,
                    base_energy_key=matched_key,
                    delta_cv_eh=item.get("delta_cv_eh"),
                    e_ae_cv_eh=item.get("e_ae_cv_eh"),
                    e_fc_cv_eh=item.get("e_fc_cv_eh"),
                    delta_rel_eh=item.get("delta_rel_eh"),
                    e_rel_eh=item.get("e_rel_eh"),
                    e_nonrel_eh=item.get("e_nonrel_eh"),
                    rel_method=item.get("rel_method"),
                    delta_zpe_eh=item.get("delta_zpe_eh"),
                    delta_thermal_h_eh=item.get("delta_thermal_h_eh"),
                    delta_thermal_g_eh=item.get("delta_thermal_g_eh"),
                    required_corrections=item_req,
                )
            )

    uniform_protocol: dict[str, str] = {}
    if any(sig for _, sig in protocol_signatures):
        first_name, first_sig = protocol_signatures[0]
        for sp_name, sig in protocol_signatures[1:]:
            if sig != first_sig:
                raise ValueError(
                    f"Inconsistent level of theory across species ({first_name!r}: {first_sig} vs "
                    f"{sp_name!r}: {sig}); relative energies require matching method/basis protocol."
                )
        uniform_protocol = dict(first_sig)

    def _extract_energy(res: CompositeEnergyResult) -> float | None:
        if res.status != "COMPLETE":
            return None
        return getattr(res, energy_field)

    complete = [(res, _extract_energy(res)) for res in evaluated if _extract_energy(res) is not None]
    pending = [res for res in evaluated if _extract_energy(res) is None]

    if not complete:
        entries = [
            RelativeEnergyEntry(
                rank=None,
                name=res.name,
                status="PENDING",
                pending_corrections=res.pending_corrections,
                e_total_eh=None,
                delta_e_eh=None,
                delta_e_kcal_mol=None,
                delta_e_kj_mol=None,
                delta_e_cm_inv=None,
                correction_status=res.correction_status,
            )
            for res in pending
        ]
        return RelativeEnergyTable(
            reference_species=None,
            reference_energy_eh=None,
            energy_field=energy_field,
            base_energy_key=uniform_energy_key,
            protocol_signature=uniform_protocol,
            entries=entries,
        )

    complete.sort(key=lambda pair: float(pair[1]))  # type: ignore[arg-type]

    if reference is not None:
        ref_matches = [pair for pair in complete if pair[0].name == reference]
        if not ref_matches:
            raise ValueError(
                f"Reference species {reference!r} is missing or PENDING; cannot anchor relative energies."
            )
        ref_name, ref_energy = ref_matches[0][0].name, float(ref_matches[0][1])  # type: ignore[arg-type]
    else:
        ref_name, ref_energy = complete[0][0].name, float(complete[0][1])  # type: ignore[arg-type]

    entries: list[RelativeEnergyEntry] = []
    for idx, (res, energy) in enumerate(complete, start=1):
        assert energy is not None
        de_eh = energy - ref_energy
        entries.append(
            RelativeEnergyEntry(
                rank=idx,
                name=res.name,
                status="COMPLETE",
                pending_corrections=res.pending_corrections,
                e_total_eh=energy,
                delta_e_eh=de_eh,
                delta_e_kcal_mol=de_eh * HARTREE_TO_KCAL_MOL,
                delta_e_kj_mol=de_eh * HARTREE_TO_KJ_MOL,
                delta_e_cm_inv=de_eh * HARTREE_TO_CM_INV,
                correction_status=res.correction_status,
            )
        )

    for res in pending:
        entries.append(
            RelativeEnergyEntry(
                rank=None,
                name=res.name,
                status="PENDING",
                pending_corrections=res.pending_corrections,
                e_total_eh=None,
                delta_e_eh=None,
                delta_e_kcal_mol=None,
                delta_e_kj_mol=None,
                delta_e_cm_inv=None,
                correction_status=res.correction_status,
            )
        )

    return RelativeEnergyTable(
        reference_species=ref_name,
        reference_energy_eh=ref_energy,
        energy_field=energy_field,
        base_energy_key=uniform_energy_key,
        protocol_signature=uniform_protocol,
        entries=entries,
    )
