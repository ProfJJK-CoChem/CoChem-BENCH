"""Exact CODATA 2018/2022 physical constants and verified CBS parameter tables.

References:
- CODATA 2018/2022 Recommended Values of the Fundamental Physical Constants (NIST).
- Helgaker, T.; Klopper, W.; Koch, H.; Noga, J. J. Chem. Phys. 1997, 106, 9639-9646.
- Halkier, A.; Helgaker, T.; Jorgensen, P.; Klopper, W.; Koch, H.; Olsen, J.; Wilson, A. K.
  Chem. Phys. Lett. 1998, 286, 243-252.
- Karton, A.; Martin, J. M. L. Theor. Chem. Acc. 2006, 115, 330-333.
- Neese, F.; Valeev, E. F. J. Chem. Theory Comput. 2011, 7, 33-43 (Table 1; ORCA defaults).
"""

from __future__ import annotations

from dataclasses import dataclass

# Exact CODATA 2018/2022 conversion factors from Hartree (E_h)
HARTREE_TO_KCAL_MOL: float = 627.5094740631
HARTREE_TO_KJ_MOL: float = 2625.4996394799
HARTREE_TO_CM_INV: float = 219474.63136320
HARTREE_TO_EV: float = 27.211386245988


@dataclass(frozen=True)
class CBSBasisParams:
    """Two-point CBS extrapolation parameters for E_SCF(X) and E_corr(X).

    E_SCF(X)  = E_SCF(inf) + A * exp(-alpha * sqrt(X))
    E_corr(X) = E_corr(inf) + B * X^(-beta)
    """

    family: str
    cardinal_pair: tuple[int, int]
    alpha: float
    beta_neese_valeev: float
    beta_helgaker: float = 3.0
    reference: str = "Neese & Valeev, JCTC 2011, 7, 33-43"


# Verified (alpha, beta) parameters keyed by (canonical_family, (X, Y))
CBS_PARAMETERS: dict[tuple[str, tuple[int, int]], CBSBasisParams] = {
    # Dunning correlation-consistent cc-pVnZ
    ("cc-pVnZ", (2, 3)): CBSBasisParams(
        family="cc-pVnZ",
        cardinal_pair=(2, 3),
        alpha=4.42,
        beta_neese_valeev=2.46,
    ),
    ("cc-pVnZ", (3, 4)): CBSBasisParams(
        family="cc-pVnZ",
        cardinal_pair=(3, 4),
        alpha=5.46,
        beta_neese_valeev=3.05,
    ),
    ("cc-pVnZ", (4, 5)): CBSBasisParams(
        family="cc-pVnZ",
        cardinal_pair=(4, 5),
        alpha=5.46,
        beta_neese_valeev=3.00,
        reference="Karton & Martin, TCA 2006, 115, 330; Helgaker et al., JCP 1997, 106, 9639",
    ),
    # Augmented Dunning aug-cc-pVnZ
    ("aug-cc-pVnZ", (2, 3)): CBSBasisParams(
        family="aug-cc-pVnZ",
        cardinal_pair=(2, 3),
        alpha=4.30,
        beta_neese_valeev=2.51,
    ),
    ("aug-cc-pVnZ", (3, 4)): CBSBasisParams(
        family="aug-cc-pVnZ",
        cardinal_pair=(3, 4),
        alpha=5.79,
        beta_neese_valeev=3.05,
    ),
    ("aug-cc-pVnZ", (4, 5)): CBSBasisParams(
        family="aug-cc-pVnZ",
        cardinal_pair=(4, 5),
        alpha=5.46,
        beta_neese_valeev=3.00,
        reference="Karton & Martin, TCA 2006, 115, 330; Helgaker et al., JCP 1997, 106, 9639",
    ),
    # Karlsruhe def2-nZVP (SVP=2, TZVP=3, QZVP=4)
    ("def2", (2, 3)): CBSBasisParams(
        family="def2",
        cardinal_pair=(2, 3),
        alpha=10.39,
        beta_neese_valeev=2.40,
    ),
    ("def2", (3, 4)): CBSBasisParams(
        family="def2",
        cardinal_pair=(3, 4),
        alpha=7.88,
        beta_neese_valeev=2.97,
    ),
    # Atomic natural orbital ano-pVnZ
    ("ano-pVnZ", (2, 3)): CBSBasisParams(
        family="ano-pVnZ",
        cardinal_pair=(2, 3),
        alpha=5.41,
        beta_neese_valeev=2.43,
    ),
    ("ano-pVnZ", (3, 4)): CBSBasisParams(
        family="ano-pVnZ",
        cardinal_pair=(3, 4),
        alpha=5.76,
        beta_neese_valeev=2.90,
    ),
}

_FAMILY_ALIASES: dict[str, str] = {
    "cc-pvnz": "cc-pVnZ",
    "cc-pvxz": "cc-pVnZ",
    "cc": "cc-pVnZ",
    "aug-cc-pvnz": "aug-cc-pVnZ",
    "aug-cc-pvxz": "aug-cc-pVnZ",
    "aug-cc": "aug-cc-pVnZ",
    "def2": "def2",
    "def2-nzvp": "def2",
    "def2-xzvp": "def2",
    "def2-tzvp/qzvp": "def2",
    "def2-svp/tzvp": "def2",
    "ano-pvnz": "ano-pVnZ",
    "ano-pvxz": "ano-pVnZ",
    "ano": "ano-pVnZ",
}


def canonicalize_basis_family(family: str) -> str:
    """Normalize a basis family string to its canonical table key."""
    key = family.strip().lower()
    if key in _FAMILY_ALIASES:
        return _FAMILY_ALIASES[key]
    raise ValueError(
        f"Unsupported basis family {family!r}. Supported families: "
        "cc-pVnZ, aug-cc-pVnZ, def2, ano-pVnZ."
    )


def get_cbs_parameters(family: str, x: int, y: int) -> CBSBasisParams:
    """Look up verified two-point CBS (alpha, beta) parameters for (family, (x, y))."""
    if not (isinstance(x, int) and isinstance(y, int)) or x < 2 or y <= x:
        raise ValueError(f"Cardinal numbers must be integers with 2 <= x < y; got ({x}, {y}).")
    canon = canonicalize_basis_family(family)
    pair = (x, y)
    if (canon, pair) not in CBS_PARAMETERS:
        available = [p for (fam, p) in CBS_PARAMETERS if fam == canon]
        raise ValueError(
            f"No verified CBS parameters for family={canon!r}, pair={pair}. "
            f"Available pairs for {canon!r}: {available}."
        )
    return CBS_PARAMETERS[(canon, pair)]
