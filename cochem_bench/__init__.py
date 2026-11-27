"""CoChem-BENCH: Complete Basis Set (CBS) Extrapolation and Ab Initio Thermochemistry."""

from __future__ import annotations

__version__ = "0.2.0"

from .constants import (
    CBS_PARAMETERS,
    HARTREE_TO_CM_INV,
    HARTREE_TO_EV,
    HARTREE_TO_KCAL_MOL,
    HARTREE_TO_KJ_MOL,
)
from .extrapolate import (
    CBSExtrapolationResult,
    calculate_f12_cbs,
    calculate_halkier,
    calculate_karton_martin,
    cbs_extrapolate,
    compute_total_cbs_energy,
)
from .orca_parser import OrcaEnergyResult, parse_orca_energy, parse_orca_output_text
from .orca_runner import (
    OrcaExecutionError,
    OrcaRunResult,
    resolve_orca_binary,
    run_orca_calculation,
    run_orca_input,
)
from .orca_writer import (
    generate_canonical_ccsd_t,
    generate_counterpoise_input,
    generate_dlpno_ccsd_f12,
    generate_dlpno_ccsd_f12_with_metadata,
    resolve_relativistic_basis,
)
from .provider import BenchProvider, bench_provider, dispatch_operation, execute_handoff
from .thermochemistry import (
    CompositeEnergyResult,
    RelativeEnergyEntry,
    RelativeEnergyTable,
    compute_composite_energy,
    compute_relative_energies,
)

__all__ = [
    "__version__",
    "CBS_PARAMETERS",
    "HARTREE_TO_KCAL_MOL",
    "HARTREE_TO_KJ_MOL",
    "HARTREE_TO_CM_INV",
    "HARTREE_TO_EV",
    "CBSExtrapolationResult",
    "calculate_halkier",
    "calculate_f12_cbs",
    "calculate_karton_martin",
    "compute_total_cbs_energy",
    "cbs_extrapolate",
    "OrcaEnergyResult",
    "parse_orca_energy",
    "parse_orca_output_text",
    "OrcaExecutionError",
    "OrcaRunResult",
    "resolve_orca_binary",
    "run_orca_input",
    "run_orca_calculation",
    "generate_dlpno_ccsd_f12",
    "generate_dlpno_ccsd_f12_with_metadata",
    "generate_canonical_ccsd_t",
    "generate_counterpoise_input",
    "resolve_relativistic_basis",
    "CompositeEnergyResult",
    "RelativeEnergyEntry",
    "RelativeEnergyTable",
    "compute_composite_energy",
    "compute_relative_energies",
    "BenchProvider",
    "bench_provider",
    "dispatch_operation",
    "execute_handoff",
]
