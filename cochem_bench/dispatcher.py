"""
CoChem-BENCH Dispatcher & Thermochemical Serializer
Validates geometries, aggregates thermochemical CBS corrections, and manages SWMR HDF5 state.
"""
import math
from pathlib import Path
from typing import List, Tuple, Dict, Any, Union
import h5py

def validate_geometry(coords: List[Union[Tuple[str, float, float, float], List[Any]]]) -> bool:
    """
    Validates molecular geometry and detects severe interatomic clashes (< 0.60 Å).
    """
    n = len(coords)
    for i in range(n):
        c1 = coords[i]
        p1 = (float(c1[1]), float(c1[2]), float(c1[3]))
        for j in range(i + 1, n):
            c2 = coords[j]
            p2 = (float(c2[1]), float(c2[2]), float(c2[3]))
            dist = math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)
            if dist < 0.60:
                raise ValueError(f"Atomic clash detected between atom {i} ({c1[0]}) and {j} ({c2[0]}): distance {dist:.3f} Å < 0.60 Å")
    return True

def compute_thermochemical_cbs(
    electronic_cbs: float,
    zpe: float = 0.0,
    thermal_h_corr: float = 0.0,
    thermal_g_corr: float = 0.0
) -> Dict[str, float]:
    """
    Combines electronic CBS energy with ZPE and thermal enthalpy/free energy corrections.
    """
    h_cbs = electronic_cbs + zpe + thermal_h_corr
    g_cbs = electronic_cbs + zpe + thermal_g_corr
    return {
        "electronic_cbs": float(electronic_cbs),
        "zpe": float(zpe),
        "thermal_h_corr": float(thermal_h_corr),
        "thermal_g_corr": float(thermal_g_corr),
        "h_cbs": float(h_cbs),
        "g_cbs": float(g_cbs)
    }

def serialize_cbs_to_hdf5(filepath: Union[str, Path], data: Dict[str, Any]) -> None:
    """
    Serializes thermochemical benchmark state to SWMR-compliant HDF5 file.
    """
    p = Path(filepath)
    p.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(p, "w", libver="latest") as f:
        grp = f.create_group("thermochemistry_cbs")
        for k, v in data.items():
            if isinstance(v, (int, float)):
                grp.attrs[k] = v
            elif isinstance(v, str):
                grp.attrs[k] = v
