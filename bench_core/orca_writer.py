#!/usr/bin/env python3
"""
CoChem-BENCH: ORCA Input Writer
Generates publication-ready ORCA execution input files supporting DLPNO-CCSD(T)-F12,
mandatory %geom convergence blocks, relativistic DKH2 injections for heavy elements,
and Boys-Bernardi counterpoise corrections with basis set linear dependency inspection.
"""

import logging
from bench_core.basis_manager import ATOMIC_NUMBERS, get_cabs_auxiliary_basis

logger = logging.getLogger("CoChem-BENCH.orca_writer")

def inspect_basis_set_cp_capability(basis: str, is_opt: bool = False) -> dict:
    """
    Inspects basis set for counterpoise (CP) capability and linear dependency risks.
    - Prohibits CP optimization on CBS extrapolation basis sets (containing '/' or 'cbs').
    - Evaluates augmentation status ('aug-').
    - Returns metadata regarding CP optimization enforcement, BSSE accuracy capping, and linear dependency risks.
    """
    mb_lower = basis.lower()
    is_cbs = ("/" in mb_lower) or ("cbs" in mb_lower)
    is_augmented = "aug" in mb_lower

    if is_opt and is_cbs:
        raise ValueError(
            f"Counterpoise optimization (CP-OPT) is prohibited on CBS extrapolation basis sets ('{basis}')."
        )

    cp_opt_enforced = is_opt and (not is_augmented)
    bsse_accuracy_capped_3pct = (not is_opt) and (not is_augmented)
    has_linear_dependency_risk = is_augmented and ("q" in mb_lower or "5" in mb_lower)

    return {
        "basis": basis,
        "is_cbs": is_cbs,
        "is_augmented": is_augmented,
        "cp_opt_enforced": cp_opt_enforced,
        "bsse_accuracy_capped_3pct": bsse_accuracy_capped_3pct,
        "has_linear_dependency_risk": has_linear_dependency_risk,
    }

def generate_dlpno_ccsd_f12(
    coords: list,
    basis: str = "aug-cc-pVTZ",
    charge: int = 0,
    mult: int = 1,
    nprocs: int = 1,
    tight_scf: bool = True,
    is_opt: bool = False,
    in_hess: str = None,
    rel_mode: str = "auto",
    maxcore: int = None,
) -> str:
    """
    Generates high-precision ORCA input for DLPNO-CCSD(T)-F12 calculations.
    Enforces mandatory 5-threshold %geom convergence for optimizations,
    DKH2 relativistic injection for heavy elements (Z > 36), and explicit CABS auxiliary space.
    """
    # Detect heavy elements (Z > 36)
    has_heavy = False
    for item in coords:
        sym_clean = item[0].replace(":", "").strip().upper()
        z = ATOMIC_NUMBERS.get(sym_clean, 1)
        if z > 36:
            has_heavy = True
            break

    rel_lower = rel_mode.lower()
    is_dkh2 = (rel_lower == "dkh2") or (rel_lower == "auto" and has_heavy)
    is_zora = (rel_lower == "zora")

    eff_basis = basis
    if (is_dkh2 or is_zora) and has_heavy:
        if eff_basis.lower() == "aug-cc-pvtz":
            eff_basis = "aug-cc-pwCVTZ-DK"
        elif not eff_basis.lower().endswith("-dk"):
            eff_basis = f"{basis}-DK"

    # Construct header ! line
    header_keywords = ["DLPNO-CCSD(T)-F12", eff_basis, "Grid5 FinalGrid6"]
    if is_zora:
        header_keywords.append("ZORA")
    if tight_scf:
        header_keywords.append("TightSCF")
    if is_opt:
        header_keywords.append("Opt")

    lines = ["! " + " ".join(header_keywords)]

    if is_dkh2:
        lines.append("# Relativistic DKH2 treatment enabled")
        lines.append("%rel Relativistic DKH2 end")
    elif is_zora:
        lines.append("# Relativistic ZORA treatment enabled")
        lines.append("%rel Relativistic ZORA end")

    if has_heavy or is_dkh2 or is_zora:
        lines.append("%core CVCut 1.0e-5 end")

    if nprocs > 1:
        lines.append(f"%pal nprocs {nprocs} end")

    if maxcore is not None:
        lines.append(f"%maxcore {maxcore}")

    cabs_basis = get_cabs_auxiliary_basis(basis)
    lines.append("%basis")
    lines.append(f'  CABS "{cabs_basis}"')
    lines.append("end")

    if is_opt:
        hess_val = in_hess if in_hess else "XTB2"
        lines.append("%geom")
        lines.append("  TolE 1e-7")
        lines.append("  TolRMSG 3e-6")
        lines.append("  TolMaxG 1e-5")
        lines.append("  TolRMSD 5e-5")
        lines.append("  TolMaxD 1e-4")
        lines.append(f"  InHess {hess_val}")
        lines.append("end")

    lines.append(f"* xyz {charge} {mult}")
    for item in coords:
        sym = item[0]
        x, y, z = item[1], item[2], item[3]
        lines.append(f"{sym:<2} {x:12.6f} {y:12.6f} {z:12.6f}")
    lines.append("*")

    return "\n".join(lines)

def generate_counterpoise_input(
    frag_a_coords: list,
    frag_b_coords: list,
    basis: str = "aug-cc-pVTZ",
    charge: int = 0,
    mult: int = 1,
    nprocs: int = 1,
    tight_scf: bool = True,
    is_opt: bool = False,
    in_hess: str = None,
    rel_mode: str = "auto",
) -> dict:
    """
    Generates ORCA input files for Boys-Bernardi counterpoise correction calculations.
    Returns a dictionary containing input strings for complex, monomer A with ghost B, monomer B with ghost A,
    and metadata flags cp_opt_enforced and bsse_accuracy_capped_3pct.
    """
    cp_info = inspect_basis_set_cp_capability(basis, is_opt=is_opt)

    def _build_coords(real_coords, ghost_coords):
        res = []
        for sym, x, y, z in real_coords:
            res.append((sym, x, y, z))
        for sym, x, y, z in ghost_coords:
            res.append((f"{sym}:", x, y, z))
        return res

    complex_coords = _build_coords(frag_a_coords + frag_b_coords, [])
    monomer_a_coords = _build_coords(frag_a_coords, frag_b_coords)
    monomer_b_coords = _build_coords(frag_b_coords, frag_a_coords)

    def _gen_input(coords):
        return generate_dlpno_ccsd_f12(
            coords,
            basis=basis,
            charge=charge,
            mult=mult,
            nprocs=nprocs,
            tight_scf=tight_scf,
            is_opt=is_opt,
            in_hess=in_hess,
            rel_mode=rel_mode,
        )

    complex_input = _gen_input(complex_coords)
    monomer_a_input = _gen_input(monomer_a_coords)
    monomer_b_input = _gen_input(monomer_b_coords)

    return {
        "complex_input": complex_input,
        "monomer_a_input": monomer_a_input,
        "monomer_b_input": monomer_b_input,
        "cp_opt_enforced": cp_info["cp_opt_enforced"],
        "bsse_accuracy_capped_3pct": cp_info["bsse_accuracy_capped_3pct"],
    }

def compute_counterpoise_corrected_energy(
    e_complex: float,
    e_monomer_a_cp: float,
    e_monomer_b_cp: float,
    e_monomer_a_raw: float,
    e_monomer_b_raw: float,
) -> dict:
    """
    Computes Boys-Bernardi counterpoise corrected interaction energy and BSSE.
    """
    delta_e_cp = e_complex - (e_monomer_a_cp + e_monomer_b_cp)
    e_bsse = (e_monomer_a_raw - e_monomer_a_cp) + (e_monomer_b_raw - e_monomer_b_cp)
    delta_e_raw = e_complex - (e_monomer_a_raw + e_monomer_b_raw)
    return {
        "delta_e_cp_hartree": delta_e_cp,
        "e_bsse_hartree": e_bsse,
        "delta_e_raw_hartree": delta_e_raw,
    }
