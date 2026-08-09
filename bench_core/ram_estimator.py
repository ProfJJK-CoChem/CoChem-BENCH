#!/usr/bin/env python3
"""
CoChem-BENCH: Analytical RAM and Hardware Estimator
Provides analytical memory predictions, scratch disk space validation, and hardware fallbacks.
"""

import os
import shutil
import logging

logger = logging.getLogger("CoChem-BENCH.ram_estimator")

def estimate_ccsd_memory(n_basis: int, n_elec: int, n_aux: int = None, n_procs: int = 1, is_open_shell: bool = False) -> float:
    """
    Resolves BENCH-03: Analytical CCSD(T) memory estimation.
    Accounts for occupied (O) vs virtual (V) orbitals, N_aux auxiliary basis functions in RI-CCSD(T),
    spin-orbitals vs spatial orbitals, 3-center integral storage buffers, and per-core MPI thread footprint.
    Returns estimated total peak memory in GB.
    """
    n_occ = n_elec // 2 if not is_open_shell else n_elec
    n_virt = max(0, n_basis - n_occ)
    
    if n_aux is None:
        n_aux = 3 * n_basis  # Standard ratio for auxiliary fitting basis sets

    # Floating point numbers (double precision = 8 bytes)
    BYTES_PER_FLOAT = 8.0

    # 1. T2 Amplitudes tensor size: O^2 * V^2 doubles (or 4 * O^2 * V^2 for spin-orbitals)
    spin_factor = 4.0 if is_open_shell else 1.0
    t2_amplitudes_bytes = spin_factor * (n_occ**2) * (n_virt**2) * BYTES_PER_FLOAT

    # 2. 3-center integrals storage: (ij|A) -> O * V * N_aux
    three_center_bytes = n_occ * n_virt * n_aux * BYTES_PER_FLOAT

    # 3. Intermediate W and T1 matrices + (T) triples buffer
    triples_buffer_bytes = n_occ**3 * n_virt * BYTES_PER_FLOAT

    # Total per-thread memory in bytes
    per_thread_bytes = t2_amplitudes_bytes + 2.0 * three_center_bytes + triples_buffer_bytes

    # Total parallel memory (accounting for MPI per-process overhead ~500MB)
    total_memory_bytes = (per_thread_bytes * n_procs) + (500 * 1024 * 1024 * n_procs)
    
    memory_gb = total_memory_bytes / (1024.0**3)
    return memory_gb

def check_scratch_disk_space(required_gb: float, scratch_path: str = None) -> bool:
    """
    Resolves BENCH-11: Verifies available disk space on scratch directory ($ORCA_TMPDIR or system temp).
    Returns True if available space >= required_gb, else False.
    """
    if scratch_path is None:
        scratch_path = os.environ.get("ORCA_TMPDIR", os.environ.get("COCHEM_ARTIFACT_DIR", "."))
        
    p = os.path.abspath(scratch_path)
    if not os.path.exists(p):
        p = os.path.dirname(p)
        
    try:
        total, used, free = shutil.disk_usage(p)
        free_gb = free / (1024.0**3)
        if free_gb < required_gb:
            logger.warning(f"Scratch disk space warning: {free_gb:.2f} GB free, {required_gb:.2f} GB required at {p}")
            return False
        return True
    except Exception as e:
        logger.error(f"Error checking scratch disk space: {e}")
        return True

def set_storage_strategy(memory_gb: float, node_max_gb: float, n_procs: int = 1) -> dict:
    """
    Resolves BENCH-12: Dynamic Hardware-Aware Fallbacks.
    Returns execution strategy dictionary with PNO thresholds, memory limits, and storage type.
    """
    max_alloc_per_core = int((node_max_gb * 1024 * 0.85) / max(1, n_procs))
    
    if memory_gb < (0.90 * node_max_gb):
        return {
            "strategy": "INCORE",
            "method_prefix": "CCSD(T)",
            "maxcore_mb": max_alloc_per_core,
            "pno_cutoff": None,
            "swap_fallback_triggered": False
        }
    else:
        # BENCH-12 fix: Mandates TightPNO (TCutPNO=1e-7) when degrading from canonical to DLPNO
        logger.warning(f"Memory limit exceeded ({memory_gb:.2f} GB required vs {node_max_gb:.2f} GB node limit). Degrading to TightPNO DLPNO-CCSD(T).")
        return {
            "strategy": "RIJCOSX_DISK",
            "method_prefix": "DLPNO-CCSD(T)",
            "maxcore_mb": max_alloc_per_core,
            "pno_cutoff": "TightPNO",  # TCutPNO=1e-7
            "pno_tcut": 1e-7,
            "swap_fallback_triggered": True
        }
