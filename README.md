# CoChem-BENCH

**CoChem-BENCH** is the Basis Set Limit and Coupled-Cluster benchmark engine for the extended CoChem suite.

It is responsible for:
- Executing ultra-high accuracy explicitly correlated DLPNO-CCSD(T)-F12 computations.
- Performing formal $1/X^3$ Complete Basis Set (CBS) extrapolations by aggregating multi-tier SCF and correlation calculations.
- Enforcing core-valence basis sets (e.g., aug-cc-pCVTZ) when called by downstream modules like CoChem-SHIFT.
- Dynamically profiling memory (RAM) and swapping to Cholesky decomposition or disk-based RIJCOSX storage to prevent memory overflow on massive nodes.

## Usage
Please refer to the authoritative `CoChem_Master_User_Manual.md` located in the `CoChem-BASE` repository for full execution instructions across the entire pipeline.