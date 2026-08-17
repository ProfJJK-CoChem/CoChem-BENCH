#!/usr/bin/env python3
"""
PyTest suite for CoChem-BENCH RAM Estimator
Resolves BENCH-02: Validates analytical memory prediction, storage strategy switching, and disk checks.
"""

import pytest

from cochem_bench.ram_estimator import estimate_ccsd_memory, set_storage_strategy

def test_ram_estimation_analytical() -> None:
    # BENCH-03: Small system memory test
    mem_small: float = estimate_ccsd_memory(n_basis=50, n_elec=10, n_procs=1)
    assert isinstance(mem_small, float)
    assert mem_small > 0.0
    
    # Large system memory test (2000 basis functions)
    mem_large: float = estimate_ccsd_memory(n_basis=2000, n_elec=50, n_procs=16)
    assert mem_large > 100.0  # Should require >100 GB

def test_storage_strategy_fallback() -> None:
    # BENCH-12: Test fallbacks and TightPNO assignment
    strat_incore: dict[str, str | bool] = set_storage_strategy(memory_gb=10.0, node_max_gb=64.0, n_procs=8)
    assert strat_incore["strategy"] == "INCORE"
    assert strat_incore["swap_fallback_triggered"] is False

    strat_disk: dict[str, str | bool] = set_storage_strategy(memory_gb=70.0, node_max_gb=64.0, n_procs=8)
    assert strat_disk["strategy"] == "RIJCOSX_DISK"
    assert strat_disk["pno_cutoff"] == "TightPNO"
    assert strat_disk["swap_fallback_triggered"] is True
