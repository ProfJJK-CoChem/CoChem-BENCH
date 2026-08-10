#!/usr/bin/env python3
"""
CoChem-BENCH: Output Parser Wrapper Alias.
Re-exports output log parsing and diagnostic evaluation routines from verifier for alias path compatibility.
"""

from bench_core.verifier import (
    parse_multireference_diagnostics,
    verify_multireference_gate,
    parse_orca_output_file,
    verify_openmpi_shm_checksum
)

__all__ = [
    "parse_multireference_diagnostics",
    "verify_multireference_gate",
    "parse_orca_output_file",
    "verify_openmpi_shm_checksum"
]
