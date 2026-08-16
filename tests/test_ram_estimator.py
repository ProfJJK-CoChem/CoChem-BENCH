#!/usr/bin/env python3
"""
PyTest suite for CoChem-BENCH RAM Estimator
Resolves BENCH-02: Validates analytical memory prediction, storage strategy switching, and disk checks.
"""

import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from cochem_bench.ram_estimator import estimate_ccsd_memory, set_storage_strategy, check_scratch_disk_space

def test_ram_estimation_analytical() -> None:
    # BENCH-03: Small system memory test
    mem_small = estimate_ccsd_memory(n_basis=50, n_elec=10, n_procs=1)
    assert isinstance(mem_small, float)
    assert mem_small > 0.0
    
    # Large system memory test (2000 basis functions)
    mem_large = estimate_ccsd_memory(n_basis=2000, n_elec=50, n_procs=16)
    assert mem_large > 100.0  # Should require >100 GB

def test_storage_strategy_fallback() -> None:
    # BENCH-12: Test fallbacks and TightPNO assignment
    strat_incore = set_storage_strategy(memory_gb=10.0, node_max_gb=64.0, n_procs=8)
    assert strat_incore["strategy"] == "INCORE"
    assert strat_incore["swap_fallback_triggered"] is False

    strat_disk = set_storage_strategy(memory_gb=70.0, node_max_gb=64.0, n_procs=8)
    assert strat_disk["strategy"] == "RIJCOSX_DISK"
    assert strat_disk["pno_cutoff"] == "TightPNO"
    assert strat_disk["swap_fallback_triggered"] is True

def test_check_scratch_disk_space() -> None:
    # BENCH-11: Check scratch disk space validation
    has_space = check_scratch_disk_space(required_gb=0.001)
    assert has_space is True
