#!/usr/bin/env python3
"""
CoChem-BENCH: Analytical RAM and Hardware Estimator
Provides analytical memory predictions, scratch disk space validation, and hardware fallbacks.
"""

import os
import shutil
import logging

logger = logging.getLogger("CoChem-BENCH.ram_estimator")

def estimate_ccsd_memory(
    n_basis: int, 
    n_elec: int, 
    n_aux: int = None, 
    n_procs: int = 1, 
    is_open_shell: bool = False,
    n_atoms: int = 2,
    engine: str = "ORCA",
    calculation_type: str = "ENERGY",
    gpu_available: bool = False,
    gpu_mem_gb: float = 16.0
) -> float:
    """
    Analytical CCSD(T) RAM & GPU memory estimator (§8.2, §9).
    Accounts for occupied/virtual orbitals, RI-CCSD(T) auxiliary basis, CFOUR analytic 2nd derivatives (36N^2 arithmetic),
    and gpu4pyscf FP64 VRAM limits.
    Returns total estimated CPU RAM in GB.
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
    triples_buffer_bytes = (n_occ**3) * n_virt * BYTES_PER_FLOAT

    # Total per-thread memory in bytes
    per_thread_bytes = t2_amplitudes_bytes + 2.0 * three_center_bytes + triples_buffer_bytes

    # Total parallel memory (accounting for MPI per-process overhead ~500MB)
    total_memory_bytes = (per_thread_bytes * n_procs) + (500 * 1024 * 1024 * n_procs)

    # CFOUR Analytic 2nd Derivatives (Hessian) 36N^2 arithmetic overhead
    if engine.upper() == "CFOUR" and calculation_type.upper() == "HESSIAN":
        cfour_36n2_bytes = 36.0 * (n_atoms**2) * (n_occ**2) * (n_virt**2) * BYTES_PER_FLOAT
        total_memory_bytes += cfour_36n2_bytes

    memory_gb = total_memory_bytes / (1024.0**3)
    return memory_gb

def estimate_ccsd_memory_detailed(
    n_basis: int, 
    n_elec: int, 
    n_aux: int = None, 
    n_procs: int = 1, 
    is_open_shell: bool = False,
    n_atoms: int = 2,
    engine: str = "ORCA",
    calculation_type: str = "ENERGY",
    gpu_available: bool = False,
    gpu_mem_gb: float = 16.0
) -> dict:
    """
    Detailed CCSD(T) memory estimation returning CPU, GPU VRAM, CFOUR overhead, and storage strategy.
    """
    cpu_gb = estimate_ccsd_memory(
        n_basis=n_basis, n_elec=n_elec, n_aux=n_aux, n_procs=n_procs,
        is_open_shell=is_open_shell, n_atoms=n_atoms, engine=engine,
        calculation_type=calculation_type, gpu_available=gpu_available, gpu_mem_gb=gpu_mem_gb
    )
    
    n_occ = n_elec // 2 if not is_open_shell else n_elec
    n_virt = max(0, n_basis - n_occ)
    BYTES_PER_FLOAT = 8.0

    if engine.upper() == "CFOUR" and calculation_type.upper() == "HESSIAN":
        cfour_36n2_bytes = 36.0 * (n_atoms**2) * (n_occ**2) * (n_virt**2) * BYTES_PER_FLOAT
        cfour_36n2_memory_gb = cfour_36n2_bytes / (1024.0**3)
    else:
        cfour_36n2_memory_gb = 0.0

    # gpu4pyscf FP64 VRAM limit calculation (N_basis^4 * 8 bytes)
    gpu_vram_gb = ((n_basis**4) * 8.0) / (1024.0**3)
    gpu_overflow = gpu_available and (gpu_vram_gb > gpu_mem_gb)

    recommendation = "INCORE"
    if cpu_gb > 32.0:
        recommendation = "RIJCOSX_DISK"
    if gpu_overflow:
        recommendation = "CPU_FALLBACK"

    return {
        "cpu_memory_gb": cpu_gb,
        "gpu_vram_gb": gpu_vram_gb,
        "gpu_overflow": gpu_overflow,
        "cfour_36n2_memory_gb": cfour_36n2_memory_gb,
        "engine": engine,
        "calculation_type": calculation_type,
        "recommendation": recommendation
    }

def check_scratch_disk_space(required_gb: float, scratch_path: str = None) -> bool:
    """
    Resolves BENCH-11 & Suggestion 11: Verifies available disk space on scratch directory ($ORCA_TMPDIR or system temp).
    Returns True if available space >= required_gb, else False (and False on exception).
    """
    if scratch_path is None:
        scratch_path = os.environ.get("ORCA_TMPDIR", os.environ.get("COCHEM_ARTIFACT_DIR", "."))
        
    p = os.path.abspath(scratch_path)
    if not os.path.exists(p):
        logger.warning(f"Scratch path {p} does not exist.")
        return False
        
    try:
        total, used, free = shutil.disk_usage(p)
        free_gb = free / (1024.0**3)
        if free_gb < required_gb:
            logger.warning(f"Scratch disk space warning: {free_gb:.2f} GB free, {required_gb:.2f} GB required at {p}")
            return False
        return True
    except Exception as e:
        logger.warning(f"Unable to verify scratch disk capacity at {p} due to exception: {e}. Defaulting scratch space check to False for safety.")
        return False

def calculate_dynamic_maxcore(nprocs: int = 1, node_max_gb: float = None, safety_factor: float = 0.85) -> int:
    """
    Suggestion 15: Calculates dynamic %maxcore allocation per core in MB based on available system RAM.
    Ensures total RAM requested (maxcore * nprocs) never exceeds physical RAM.
    """
    if node_max_gb is None:
        try:
            import psutil
            node_max_gb = psutil.virtual_memory().total / (1024.0**3)
        except Exception:
            node_max_gb = 16.0
    
    nprocs = max(1, nprocs)
    max_alloc_per_core = int((node_max_gb * 1024 * safety_factor) / nprocs)
    max_physical_per_core = int((node_max_gb * 1024) / nprocs)
    
    floor_mb = min(256, max_physical_per_core)
    floor_mb = max(1, floor_mb)
    
    res = max(floor_mb, max_alloc_per_core)
    if res * nprocs > node_max_gb * 1024:
        res = max_physical_per_core
    return max(1, res)

def set_storage_strategy(memory_gb: float, node_max_gb: float, n_procs: int = 1) -> dict:
    """
    Resolves BENCH-12: Dynamic Hardware-Aware Fallbacks.
    Returns execution strategy dictionary with PNO thresholds, memory limits, and storage type.
    """
    max_alloc_per_core = calculate_dynamic_maxcore(nprocs=n_procs, node_max_gb=node_max_gb)
    
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
