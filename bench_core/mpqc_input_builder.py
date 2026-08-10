#!/usr/bin/env python3
"""
CoChem-BENCH: MPQC Input Builder Wrapper Alias.
Re-exports input generation routines from mpqc_writer for alias path compatibility.
"""

from bench_core.mpqc_writer import (
    generate_mpqc_ccsd_f12,
    generate_counterpoise_input,
    compute_counterpoise_corrected_energy
)

__all__ = [
    "generate_mpqc_ccsd_f12",
    "generate_counterpoise_input",
    "compute_counterpoise_corrected_energy"
]
