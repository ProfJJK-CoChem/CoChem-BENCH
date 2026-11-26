"""
CoChem-BENCH Basis Set Manager
Handles full 118-element periodic-table lookup, core-valence basis mapping,
CABS auxiliary basis selection, and scalar-relativistic (-DK / DKH / ZORA / x2c) resolution.
"""

from __future__ import annotations

# Complete IUPAC periodic table: Z = 1 (H) through Z = 118 (Og)
PERIODIC_TABLE_ATOMIC_NUMBERS: dict[str, int] = {
    "H": 1, "He": 2,
    "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38,
    "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56,
    "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88,
    "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96,
    "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
    "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
    "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
}

# Backward-compatible export covering all 118 elements (previously truncated at Xe Z=54).
HEAVY_ATOMIC_NUMBERS: dict[str, int] = PERIODIC_TABLE_ATOMIC_NUMBERS

TRANSITION_METAL_3D_MIN_Z: int = 21
TRANSITION_METAL_3D_MAX_Z: int = 30
HEAVY_ELEMENT_MIN_Z: int = 36

# Official ORCA 6.1 manual Table 2.37 coverage for cc-pVXZ-F12-CABS keywords.
ORCA6_EXPLICIT_F12_CABS_ELEMENTS: frozenset[str] = frozenset({
    "H", "B", "C", "N", "O", "F", "Ne", "Al", "Si", "P", "S", "Cl", "Ar",
})


def normalize_element_symbol(symbol: str) -> str:
    """Normalize an element or ghost-atom symbol (e.g. 'br:', 'Fe') to canonical IUPAC form."""
    if not isinstance(symbol, str):
        raise TypeError(f"Symbol must be a string, got {type(symbol).__name__}")
    cleaned = symbol.strip().rstrip(":").strip()
    if not cleaned:
        raise ValueError(f"Empty chemical symbol: {symbol!r}")
    canonical = cleaned[0].upper() + cleaned[1:].lower()
    if canonical not in PERIODIC_TABLE_ATOMIC_NUMBERS:
        raise ValueError(f"Unrecognized chemical symbol: {symbol}")
    return canonical


def get_atomic_number(symbol: str) -> int:
    """Return the IUPAC atomic number Z (1..118) for a chemical symbol or ghost symbol."""
    return PERIODIC_TABLE_ATOMIC_NUMBERS[normalize_element_symbol(symbol)]


def is_3d_transition_metal(symbol: str) -> bool:
    """Return True for first-row 3d transition metals Sc (Z=21) through Zn (Z=30)."""
    z = get_atomic_number(symbol)
    return TRANSITION_METAL_3D_MIN_Z <= z <= TRANSITION_METAL_3D_MAX_Z


def is_heavy_element(symbol: str) -> bool:
    """Return True for heavy elements Kr (Z=36) through Og (Z=118)."""
    return get_atomic_number(symbol) >= HEAVY_ELEMENT_MIN_Z


def requires_scalar_relativistic(symbol: str) -> bool:
    """Return True when scalar-relativistic treatment (DKH2/ZORA) is scientifically required.

    Covers all heavy elements Z >= 36 through Z = 118, bromine (Z=35, heavy 4p halogen),
    and 3d transition metals (Z = 21..30).
    """
    z = get_atomic_number(symbol)
    return z >= 35 or (TRANSITION_METAL_3D_MIN_Z <= z <= TRANSITION_METAL_3D_MAX_Z)


def map_core_valence_basis_per_element(symbols: list[str], base_cv_basis: str = "aug-cc-pCVTZ") -> dict[str, str]:
    """
    Maps light elements (H, He) to standard correlation-consistent basis without CV,
    and heavier elements (Li Z=3 through Og Z=118) to base_cv_basis.
    """
    if not isinstance(symbols, list):
        raise TypeError("symbols must be a list of strings.")

    mapped: dict[str, str] = {}
    light_basis = base_cv_basis.replace("pCV", "pV")
    for s in symbols:
        elem = normalize_element_symbol(s)
        if elem in ("H", "He"):
            mapped[elem] = light_basis
        else:
            mapped[elem] = base_cv_basis
    return mapped


def get_cabs_auxiliary_basis(orbital_basis: str) -> str:
    """
    Returns complementary auxiliary basis set (CABS) for F12 calculations.
    """
    if not isinstance(orbital_basis, str):
        raise TypeError(f"orbital_basis must be a string, got {type(orbital_basis).__name__}")

    cabs_registry: dict[str, str] = {
        "cc-pVDZ-F12": "cc-pVDZ-F12-CABS",
        "cc-pVTZ-F12": "cc-pVTZ-F12-CABS",
        "cc-pVQZ-F12": "cc-pVQZ-F12-CABS",
        "aug-cc-pVDZ": "aug-cc-pVDZ/OptRI",
        "aug-cc-pVTZ": "aug-cc-pVTZ/OptRI",
        "aug-cc-pVQZ": "aug-cc-pVQZ/OptRI",
        "def2-SVP": "def2-SVP/C",
        "def2-TZVP": "def2-TZVP/C",
        "def2-QZVP": "def2-QZVP/C",
    }

    if orbital_basis not in cabs_registry:
        raise ValueError(f"No CABS auxiliary basis registered for {orbital_basis}")

    return cabs_registry[orbital_basis]
