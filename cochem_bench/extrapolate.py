"""CoChem-BENCH Extrapolation Algorithms.

Two-point Complete Basis Set (CBS) extrapolations:
- Exponential SCF extrapolation: E_SCF(X) = E_SCF(inf) + A * exp(-alpha * sqrt(X))
- Inverse-power correlation extrapolation: E_corr(X) = E_corr(inf) + B * X^(-beta)
  (Helgaker inverse-cubic beta = 3.0, Neese-Valeev optimal beta, and F12 beta = 7.0).
"""

from __future__ import annotations

import math
from typing import Literal

from pydantic import BaseModel

from .constants import canonicalize_basis_family, get_cbs_parameters


class CBSExtrapolationResult(BaseModel):
    """Structured result of a two-point CBS extrapolation."""

    basis_family: str
    cardinal_x: int
    cardinal_y: int
    scheme: str
    alpha: float
    beta: float
    e_scf_x_eh: float
    e_scf_y_eh: float
    e_corr_x_eh: float
    e_corr_y_eh: float
    e_scf_cbs_eh: float
    e_corr_cbs_eh: float
    e_total_cbs_eh: float
    delta_scf_eh: float
    delta_corr_eh: float
    provenance_tag: str = "[D]"


def extrapolate_scf_exponential(
    e_scf_x: float,
    e_scf_y: float,
    x: int,
    y: int,
    alpha: float,
) -> float:
    """Two-point exponential SCF CBS extrapolation.

    E_SCF(X) = E_SCF(inf) + A * exp(-alpha * sqrt(X))
    => E_SCF(inf) = E_SCF(Y) + (E_SCF(Y) - E_SCF(X)) / (exp(alpha * (sqrt(Y) - sqrt(X))) - 1)
    """
    if y <= x or x < 2:
        raise ValueError(f"Cardinal numbers must satisfy 2 <= x < y; got x={x}, y={y}")
    if alpha <= 0.0:
        raise ValueError(f"Exponential parameter alpha must be positive; got {alpha}")
    denom = math.expm1(alpha * (math.sqrt(float(y)) - math.sqrt(float(x))))
    return float(e_scf_y + (e_scf_y - e_scf_x) / denom)


def extrapolate_corr_power(
    e_corr_x: float,
    e_corr_y: float,
    x: int,
    y: int,
    beta: float = 3.0,
) -> float:
    """Two-point inverse-power correlation CBS extrapolation (Helgaker X^-beta).

    E_corr(X) = E_corr(inf) + B * X^(-beta)
    => E_corr(inf) = (Y^beta * E_corr(Y) - X^beta * E_corr(X)) / (Y^beta - X^beta)
    """
    if y <= x or x < 2:
        raise ValueError(f"Cardinal numbers must satisfy 2 <= x < y; got x={x}, y={y}")
    if beta <= 0.0:
        raise ValueError(f"Power exponent beta must be positive; got {beta}")
    cx = float(x) ** beta
    cy = float(y) ** beta
    return float((cy * e_corr_y - cx * e_corr_x) / (cy - cx))


def cbs_extrapolate(
    e_scf_x: float,
    e_scf_y: float,
    e_corr_x: float | None = None,
    e_corr_y: float | None = None,
    *,
    e_tot_x: float | None = None,
    e_tot_y: float | None = None,
    x: int = 3,
    y: int = 4,
    basis_family: str = "cc-pVnZ",
    scheme: Literal["helgaker", "neese_valeev", "f12"] = "helgaker",
    alpha: float | None = None,
    beta: float | None = None,
) -> CBSExtrapolationResult:
    """Perform verified two-point SCF + correlation CBS extrapolation.

    Either (e_corr_x, e_corr_y) or (e_tot_x, e_tot_y) must be supplied.
    If total energies are supplied instead of correlation energies,
    E_corr(n) = E_tot(n) - E_SCF(n) is computed directly.
    """
    if e_corr_x is None or e_corr_y is None:
        if e_tot_x is not None and e_tot_y is not None:
            e_corr_x = float(e_tot_x) - float(e_scf_x)
            e_corr_y = float(e_tot_y) - float(e_scf_y)
        else:
            raise ValueError(
                "Must supply either (e_corr_x, e_corr_y) or (e_tot_x, e_tot_y)."
            )

    canon_family = canonicalize_basis_family(basis_family)
    params = get_cbs_parameters(canon_family, x, y)

    eff_alpha = float(alpha) if alpha is not None else params.alpha
    if beta is not None:
        eff_beta = float(beta)
    elif scheme == "helgaker":
        eff_beta = params.beta_helgaker  # 3.0 exact inverse-cubic
    elif scheme == "neese_valeev":
        eff_beta = params.beta_neese_valeev
    elif scheme == "f12":
        eff_beta = 7.0
    else:
        raise ValueError(f"Unsupported CBS scheme {scheme!r}; use 'helgaker', 'neese_valeev', or 'f12'.")

    e_scf_cbs = extrapolate_scf_exponential(float(e_scf_x), float(e_scf_y), x, y, eff_alpha)
    e_corr_cbs = extrapolate_corr_power(float(e_corr_x), float(e_corr_y), x, y, eff_beta)
    e_tot_cbs = e_scf_cbs + e_corr_cbs

    return CBSExtrapolationResult(
        basis_family=canon_family,
        cardinal_x=x,
        cardinal_y=y,
        scheme=scheme,
        alpha=eff_alpha,
        beta=eff_beta,
        e_scf_x_eh=float(e_scf_x),
        e_scf_y_eh=float(e_scf_y),
        e_corr_x_eh=float(e_corr_x),
        e_corr_y_eh=float(e_corr_y),
        e_scf_cbs_eh=e_scf_cbs,
        e_corr_cbs_eh=e_corr_cbs,
        e_total_cbs_eh=e_tot_cbs,
        delta_scf_eh=e_scf_cbs - float(e_scf_y),
        delta_corr_eh=e_corr_cbs - float(e_corr_y),
    )


def calculate_halkier(corr_x3: float, corr_x4: float) -> float:
    """Halkier 1/X^3 two-point correlation energy CBS extrapolation for T/Q zeta.

    E_corr(CBS) = (4^3 * E(4) - 3^3 * E(3)) / (4^3 - 3^3) = (64*E(4) - 27*E(3)) / 37
    """
    return extrapolate_corr_power(corr_x3, corr_x4, x=3, y=4, beta=3.0)


def calculate_f12_cbs(corr_x3: float, corr_x4: float, beta: float = 7.0) -> float:
    """F12 1/X^beta CBS correlation energy extrapolation (default beta = 7.0)."""
    return extrapolate_corr_power(corr_x3, corr_x4, x=3, y=4, beta=beta)


def calculate_karton_martin(hf_x3: float, hf_x4: float, family: str = "aug-cc-pVXZ") -> float:
    """Karton-Martin exponential two-point HF CBS extrapolation for T/Q zeta.

    E_HF(X) = E_HF(CBS) + A * exp(-alpha * sqrt(X))
    """
    alpha = 5.46 if "aug" in family.lower() else 5.70
    return extrapolate_scf_exponential(hf_x3, hf_x4, x=3, y=4, alpha=alpha)


def compute_total_cbs_energy(
    hf_x3: float,
    hf_x4: float,
    corr_x3: float,
    corr_x4: float,
    is_f12: bool = True,
    basis_family: str = "aug-cc-pVXZ",
) -> dict[str, float]:
    """Compute total CBS energy combining HF and correlation components (T/Q)."""
    hf_cbs = calculate_karton_martin(hf_x3, hf_x4, family=basis_family)
    if is_f12:
        corr_cbs = calculate_f12_cbs(corr_x3, corr_x4, beta=7.0)
    else:
        corr_cbs = calculate_halkier(corr_x3, corr_x4)

    total_cbs = hf_cbs + corr_cbs
    return {
        "hf_cbs": hf_cbs,
        "corr_cbs": corr_cbs,
        "total_cbs": total_cbs,
    }
