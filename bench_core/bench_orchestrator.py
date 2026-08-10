#!/usr/bin/env python3
"""
CoChem-BENCH: Bench Orchestrator Wrapper Alias.
Re-exports SHM verification and thermochemical functions for alias path compatibility.
"""

from bench_core.verifier import verify_openmpi_shm_checksum
from bench_core.dispatcher import (
    validate_geometry,
    compute_thermochemical_cbs,
    serialize_cbs_to_hdf5
)

__all__ = [
    "verify_openmpi_shm_checksum",
    "validate_geometry",
    "compute_thermochemical_cbs",
    "serialize_cbs_to_hdf5"
]
