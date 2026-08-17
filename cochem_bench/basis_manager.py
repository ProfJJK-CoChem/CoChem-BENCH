"""
CoChem-BENCH Basis Set Manager
Handles core-valence basis mapping, CABS selection, relativistic variants (-DK).
"""
from typing import List, Dict

HEAVY_ATOMIC_NUMBERS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54
}

def map_core_valence_basis_per_element(symbols: List[str], base_cv_basis: str = "aug-cc-pCVTZ") -> Dict[str, str]:
    """
    Maps light elements (H, He) to standard correlation-consistent basis without CV,
    and heavier elements to base_cv_basis.
    """
    mapped = {}
    light_basis = base_cv_basis.replace("pCV", "pV")
    for s in symbols:
        elem = s.capitalize()
        if elem in ["H", "He"]:
            mapped[elem] = light_basis
        else:
            mapped[elem] = base_cv_basis
    return mapped

def get_cabs_auxiliary_basis(orbital_basis: str) -> str:
    """
    Returns complementary auxiliary basis set (CABS) for F12 calculations.
    """
    return f"{orbital_basis}/C"
