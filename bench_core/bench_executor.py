#!/usr/bin/env python3
"""
CoChem-BENCH: Bench Executor Wrapper Alias.
Re-exports memory estimation and storage strategy functions from ram_estimator for alias path compatibility.
"""

from bench_core.ram_estimator import (
    calculate_dynamic_maxcore,
    estimate_ccsd_memory,
    set_storage_strategy,
    check_scratch_disk_space
)

__all__ = [
    "calculate_dynamic_maxcore",
    "estimate_ccsd_memory",
    "set_storage_strategy",
    "check_scratch_disk_space"
]
