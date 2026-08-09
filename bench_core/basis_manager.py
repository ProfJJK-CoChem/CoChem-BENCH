#!/usr/bin/env python3
"""
CoChem-BENCH: Basis Set Manager
Enforces element-wise core-valence basis set mappings and CABS space auxiliary assignments.
"""

import logging

logger = logging.getLogger("CoChem-BENCH.basis_manager")

# Atomic numbers table
ATOMIC_NUMBERS = {
    'H': 1, 'HE': 2, 'LI': 3, 'BE': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NE': 10,
    'NA': 11, 'MG': 12, 'AL': 13, 'SI': 14, 'P': 15, 'S': 16, 'CL': 17, 'AR': 18,
    'K': 19, 'CA': 20, 'SC': 21, 'TI': 22, 'V': 23, 'CR': 24, 'MN': 25, 'FE': 26,
    'CO': 27, 'NI': 28, 'CU': 29, 'ZN': 30, 'GA': 31, 'GE': 32, 'AS': 33, 'SE': 34,
    'BR': 35, 'KR': 36, 'RB': 37, 'SR': 38, 'Y': 39, 'ZR': 40, 'NB': 41, 'MO': 42,
    'TC': 43, 'RU': 44, 'RH': 45, 'PD': 46, 'AG': 47, 'CD': 48, 'IN': 49, 'SN': 50,
    'SB': 51, 'TE': 52, 'I': 53, 'XE': 54
}

def map_core_valence_basis_per_element(symbols: list, base_cv_basis: str = "aug-cc-pCVTZ") -> dict:
    """
    Resolves BENCH-07: Atomic-number aware core-valence basis set mapping.
    Heavy elements (Z >= 3) receive requested core-valence basis set (e.g. aug-cc-pCVTZ).
    Light elements H (Z=1) and He (Z=2) have no core electrons and no aug-cc-pCVTZ definition;
    they are mapped to standard valence basis sets (e.g. aug-cc-pVTZ).
    """
    element_basis_map = {}
    
    # Determine corresponding valence fallback basis name
    valence_fallback = base_cv_basis.replace("-pCV", "-pV").replace("CV", "V")
    
    for sym in symbols:
        sym_clean = sym.strip().upper()
        z = ATOMIC_NUMBERS.get(sym_clean, 6)  # Default Z=6 if unknown symbol
        
        if z < 3:
            # H (Z=1) and He (Z=2): map to aug-cc-pVTZ
            element_basis_map[sym] = valence_fallback
            logger.info(f"Element {sym} (Z={z}) mapped to valence basis {valence_fallback} (no core electrons).")
        else:
            # Heavy elements Z >= 3: map to aug-cc-pCVTZ
            element_basis_map[sym] = base_cv_basis
            
    return element_basis_map

def get_cabs_auxiliary_basis(main_basis: str) -> str:
    """
    Resolves BENCH-10: Assigns explicit CABS space auxiliary basis set for F12 calculations.
    Matching:
    - aug-cc-pVTZ / aug-cc-pVTZ-F12 -> aug-cc-pVTZ-CABS / cc-pVTZ-F12-CABS
    - aug-cc-pVQZ / aug-cc-pVQZ-F12 -> aug-cc-pVQZ-CABS / cc-pVQZ-F12-CABS
    """
    mb_lower = main_basis.lower()
    if "qv" in mb_lower or "qz" in mb_lower or "4" in mb_lower:
        return "aug-cc-pVQZ-CABS"
    elif "tv" in mb_lower or "tz" in mb_lower or "3" in mb_lower:
        return "aug-cc-pVTZ-CABS"
    else:
        return "cc-pVTZ-F12-CABS"
