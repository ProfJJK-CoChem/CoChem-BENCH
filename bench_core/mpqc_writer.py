#!/usr/bin/env python3
"""
CoChem-BENCH: MPQC Input Writer
Generates publication-ready MPQC execution input files (.json) supporting CCSD(T)-F12.
"""

import json
import logging
from bench_core.basis_manager import ATOMIC_NUMBERS

logger = logging.getLogger("CoChem-BENCH.mpqc_writer")

def generate_mpqc_ccsd_f12(
    coords: list, 
    basis: str = "cc-pVTZ-F12", 
    charge: int = 0, 
    mult: int = 1
) -> str:
    """
    Constructs the complete MPQC json input string.
    """
    atoms = []
    for item in coords:
        sym = item[0]
        z = ATOMIC_NUMBERS.get(sym.strip().upper(), 1)
        atoms.append({"Z": z, "r": [item[1], item[2], item[3]]})
        
    input_dict = {
        "molecule": {
            "atoms": atoms,
            "charge": charge,
            "multiplicity": mult
        },
        "basis": basis,
        "method": "CCSD(T)-F12"
    }
    return json.dumps(input_dict, indent=2)

def generate_counterpoise_input(
    frag_a_coords: list, 
    frag_b_coords: list, 
    basis: str = "cc-pVTZ-F12", 
    charge: int = 0, 
    mult: int = 1
) -> dict:
    """
    Counterpoise input generation for MPQC JSON format.
    Generates inputs for Complex (AB), Monomer A with ghost B (A:B), and Monomer B with ghost A (B:A).
    """
    def format_coords(coords, ghost_coords):
        atoms = []
        for sym, x, y, z in coords:
            z_num = ATOMIC_NUMBERS.get(sym.strip().upper(), 1)
            atoms.append({"Z": z_num, "r": [x, y, z]})
        for sym, x, y, z in ghost_coords:
            z_num = ATOMIC_NUMBERS.get(sym.strip().upper(), 1)
            atoms.append({"Z": z_num, "r": [x, y, z], "ghost": True})
            
        input_dict = {
            "molecule": {
                "atoms": atoms,
                "charge": charge,
                "multiplicity": mult
            },
            "basis": basis,
            "method": "CCSD(T)-F12"
        }
        return json.dumps(input_dict, indent=2)

    complex_input = format_coords(frag_a_coords + frag_b_coords, [])
    monomer_a_cp = format_coords(frag_a_coords, [(c[0], c[1], c[2], c[3]) for c in frag_b_coords])
    monomer_b_cp = format_coords(frag_b_coords, [(c[0], c[1], c[2], c[3]) for c in frag_a_coords])

    return {
        "complex_input": complex_input,
        "monomer_a_input": monomer_a_cp,
        "monomer_b_input": monomer_b_cp
    }

def compute_counterpoise_corrected_energy(
    e_complex: float,
    e_monomer_a_cp: float,
    e_monomer_b_cp: float,
    e_monomer_a_raw: float = None,
    e_monomer_b_raw: float = None
) -> dict:
    delta_e_cp = e_complex - e_monomer_a_cp - e_monomer_b_cp
    hartree_to_kcal = 627.509474
    
    res = {
        "delta_e_cp_hartree": delta_e_cp,
        "delta_e_cp_kcal": delta_e_cp * hartree_to_kcal
    }
    
    if e_monomer_a_raw is not None and e_monomer_b_raw is not None:
        delta_e_raw = e_complex - e_monomer_a_raw - e_monomer_b_raw
        bsse_a = e_monomer_a_raw - e_monomer_a_cp
        bsse_b = e_monomer_b_raw - e_monomer_b_cp
        e_bsse = bsse_a + bsse_b
        res["delta_e_raw_hartree"] = delta_e_raw
        res["delta_e_raw_kcal"] = delta_e_raw * hartree_to_kcal
        res["e_bsse_hartree"] = e_bsse
        res["e_bsse_kcal"] = e_bsse * hartree_to_kcal

    return res
