#!/usr/bin/env python3
"""
CoChem-BENCH: Complete Basis Set (CBS) Extrapolation Algorithms
Implements Karton-Martin HF extrapolation and Halkier / F12 correlation extrapolation formulas.
"""

import numpy as np

def calculate_halkier(corr_x3: float, corr_x4: float, x3: int = 3, x4: int = 4) -> float:
    """
    Standard orbital $1/X^3$ Halkier CBS correlation extrapolation for canonical CCSD(T).
    Formula: E_corr(CBS) = (x4^3 * E_corr(x4) - x3^3 * E_corr(x3)) / (x4^3 - x3^3)
    """
    denominator = (x4**3 - x3**3)
    cbs_corr = (x4**3 * corr_x4 - x3**3 * corr_x3) / denominator
    return float(cbs_corr)

def calculate_f12_cbs(corr_x3: float, corr_x4: float, x3: int = 3, x4: int = 4, beta: float = 7.0) -> float:
    """
    Resolves BENCH-04: F12-specific correlation energy extrapolation.
    Explicitly correlated DLPNO-CCSD(T)-F12 correlation energy converges as ~1/X^7 or exponential.
    Applying 1/X^3 to F12 over-extrapolates.
    Formula: E_corr(CBS) = (x4^beta * E_corr(x4) - x3^beta * E_corr(x3)) / (x4^beta - x3^beta)
    """
    denominator = (x4**beta - x3**beta)
    cbs_corr = (x4**beta * corr_x4 - x3**beta * corr_x3) / denominator
    return float(cbs_corr)

def calculate_karton_martin(hf_x3: float, hf_x4: float, family: str = "aug-cc-pVXZ", x3: int = 3, x4: int = 4) -> float:
    """
    Resolves BENCH-15: Parameterized Karton-Martin Hartree-Fock SCF CBS extrapolation.
    E_HF(X) = E_HF(CBS) + A * exp(-alpha * sqrt(X)) or E_HF(CBS) + A * X^-alpha.
    Alpha exponents by basis set family:
    - cc-pVXZ: alpha = 1.43
    - aug-cc-pVXZ: alpha = 1.56
    - cc-pCVXZ / aug-cc-pCVTZ: alpha = 1.50
    - def2-TZVP / def2-QZVP: alpha = 1.62
    """
    family_lower = family.lower()
    if "aug-cc-pcv" in family_lower or "pcv" in family_lower:
        alpha = 1.50
    elif "aug-cc" in family_lower:
        alpha = 1.56
    elif "def2" in family_lower:
        alpha = 1.62
    elif "cc-p" in family_lower:
        alpha = 1.43
    else:
        alpha = 1.50

    # Karton-Martin exponential extrapolation:
    # (E_HF(x4) - E_HF(x3)) = A * (exp(-alpha*sqrt(x4)) - exp(-alpha*sqrt(x3)))
    denom = np.exp(-alpha * np.sqrt(x4)) - np.exp(-alpha * np.sqrt(x3))
    if abs(denom) < 1e-12:
        return hf_x4
        
    a_param = (hf_x4 - hf_x3) / denom
    hf_cbs = hf_x4 - a_param * np.exp(-alpha * np.sqrt(x4))
    return float(hf_cbs)

def compute_total_cbs_energy(hf_x3: float, hf_x4: float, corr_x3: float, corr_x4: float, 
                             is_f12: bool = True, basis_family: str = "aug-cc-pVXZ") -> dict:
    """
    Combines HF SCF CBS and correlation CBS energies.
    Returns dictionary with HF_CBS, Corr_CBS, and Total_CBS.
    """
    hf_cbs = calculate_karton_martin(hf_x3, hf_x4, family=basis_family)
    
    if is_f12:
        corr_cbs = calculate_f12_cbs(corr_x3, corr_x4)
    else:
        corr_cbs = calculate_halkier(corr_x3, corr_x4)
        
    total_cbs = hf_cbs + corr_cbs
    return {
        "hf_cbs": hf_cbs,
        "corr_cbs": corr_cbs,
        "total_cbs": total_cbs,
        "is_f12": is_f12,
        "basis_family": basis_family
    }
