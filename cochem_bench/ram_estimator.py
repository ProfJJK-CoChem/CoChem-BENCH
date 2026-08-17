"""
CoChem-BENCH RAM & Resource Estimator
Estimates CCSD memory, storage strategy switching, GPU VRAM requirements, and scratch validation.
"""
import shutil
from pathlib import Path
from typing import Dict, Any, Optional

def estimate_ccsd_memory(
    n_basis: int,
    n_elec: int,
    n_procs: int = 1,
    n_atoms: int = 1,
    engine: str = "ORCA",
    calculation_type: str = "ENERGY",
    is_open_shell: bool = False
) -> float:
    """
    Estimates memory required for CCSD / DLPNO-CCSD calculation in GB.
    """
    det = estimate_ccsd_memory_detailed(
        n_basis=n_basis,
        n_elec=n_elec,
        n_procs=n_procs,
        n_atoms=n_atoms,
        engine=engine,
        calculation_type=calculation_type,
        is_open_shell=is_open_shell
    )
    return det["cpu_memory_gb"]

def estimate_ccsd_memory_detailed(
    n_basis: int,
    n_elec: int,
    n_procs: int = 1,
    n_atoms: int = 1,
    engine: str = "ORCA",
    calculation_type: str = "ENERGY",
    is_open_shell: bool = False,
    gpu_available: bool = False,
    gpu_mem_gb: float = 16.0
) -> Dict[str, Any]:
    """
    Detailed memory breakdown including CPU, GPU VRAM, and engine-specific overhead.
    """
    n_occ = max(1, n_elec // 2)
    n_vir = max(1, n_basis - n_occ)
    spin_factor = 2.0 if is_open_shell else 1.0

    # Base T2 amplitudes: N_occ^2 * N_vir^2 * 8 bytes
    t2_bytes = (n_occ ** 2) * (n_vir ** 2) * 8.0 * spin_factor
    base_cpu_gb = (t2_bytes / (1024.0 ** 3)) * 4.0 * max(1, n_procs)

    # CFOUR Hessian 36N^2 displacement coordinate factor
    cfour_36n2_memory_gb = 0.0
    if engine.upper() == "CFOUR" and calculation_type.upper() == "HESSIAN":
        cfour_bytes = 36.0 * (n_basis ** 2) * n_atoms * 8.0
        cfour_36n2_memory_gb = cfour_bytes / (1024.0 ** 3)
        base_cpu_gb += cfour_36n2_memory_gb
    elif calculation_type.upper() == "HESSIAN":
        # Other engines Hessian scaling
        base_cpu_gb += (n_atoms * 0.05)

    # GPU4PySCF full 4-center integral VRAM estimation
    gpu_vram_gb = ((n_basis ** 4) * 8.0) / (1024.0 ** 3)
    gpu_overflow = gpu_vram_gb > gpu_mem_gb
    recommendation = "CPU_FALLBACK" if gpu_overflow else "GPU_ACCELERATED"

    return {
        "cpu_memory_gb": float(base_cpu_gb),
        "gpu_vram_gb": float(gpu_vram_gb),
        "gpu_overflow": bool(gpu_overflow),
        "cfour_36n2_memory_gb": float(cfour_36n2_memory_gb),
        "recommendation": recommendation
    }

def set_storage_strategy(memory_gb: float, node_max_gb: float = 64.0, n_procs: int = 8) -> Dict[str, Any]:
    """
    Selects INCORE vs RIJCOSX_DISK storage strategy based on memory ceiling.
    """
    if memory_gb <= node_max_gb:
        return {
            "strategy": "INCORE",
            "pno_cutoff": "NormalPNO",
            "swap_fallback_triggered": False
        }
    else:
        return {
            "strategy": "RIJCOSX_DISK",
            "pno_cutoff": "TightPNO",
            "swap_fallback_triggered": True
        }

def calculate_dynamic_maxcore(nprocs: int = 4, node_max_gb: float = 16.0) -> int:
    """
    Calculates per-process maxcore in MB leaving headroom for OS and MPI runtime.
    """
    available_mb = node_max_gb * 1024.0 * 0.75
    per_proc_mb = int(available_mb / max(1, nprocs))
    return max(1000, per_proc_mb)

def check_scratch_disk_space(required_gb: float, scratch_path: Optional[str] = None) -> bool:
    """
    Checks if scratch disk space is sufficient. Safe exception handling for invalid paths.
    """
    try:
        target = Path(scratch_path) if scratch_path else Path.cwd()
        if not target.exists():
            return False
        usage = shutil.disk_usage(str(target))
        free_gb = usage.free / (1024.0 ** 3)
        return free_gb >= required_gb
    except Exception:
        return False
