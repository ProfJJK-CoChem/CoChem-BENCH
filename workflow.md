# CoChem-BENCH: Software Engineering Specification
**Target Phase:** Python Implementation

This document serves as the exact coding blueprint for the next LLM agent to construct the `CoChem-BENCH` repository.

## 1. Directory & File Architecture
```text
CoChem-BENCH/
├── bench_core/
│   ├── __init__.py
│   ├── dispatcher.py      # Entry point for BASE payload ingestion
│   ├── ram_estimator.py   # Analytical memory predictor for CCSD(T)
│   ├── extrapolate.py     # Karton-Martin and Halkier algorithms
│   └── orca_writer.py     # Generates the strict input files (.inp)
├── tests/
│   ├── test_extrapolate.py
│   └── test_ram_estimator.py
├── requirements.txt       # h5py, numpy, scipy
└── README.md
```

## 2. File-by-File Blueprint

### `bench_core/ram_estimator.py`
- **Purpose:** Predicts exact peak RAM usage to decide between In-Core vs Out-of-Core execution.
- **Functions:**
  - `def estimate_ccsd_memory(n_basis: int, n_elec: int) -> float:`
    - *Returns:* Memory in GB. Uses formula $O(O^2 V^4)$ estimation.
  - `def set_storage_strategy(memory_gb: float, node_max_gb: float) -> str:`
    - *Returns:* `"INCORE"` if `< 90%` of node limit, else `"RIJCOSX_DISK"`.

### `bench_core/extrapolate.py`
- **Purpose:** Performs Complete Basis Set (CBS) limits from X=3 and X=4 energies.
- **Functions:**
  - `def calculate_halkier(corr_x3: float, corr_x4: float) -> float:`
    - *Returns:* Limit using $1/3^3$ and $1/4^3$.
  - `def calculate_karton_martin(hf_x3: float, hf_x4: float) -> float:`
    - *Returns:* HF SCF CBS limit.

### `bench_core/orca_writer.py`
- **Purpose:** Constructs the actual ORCA execution text.
- **Functions:**
  - `def generate_dlpno_ccsd_f12(coords: list, basis: str, strategy: str) -> str:`
    - *Returns:* Complete ORCA input string enforcing DefGrid4 and AutoAux.

## 3. Execution Data Flow (The Payload Trace)
1. **Payload Ingest:** `dispatcher.py` receives JSON via `sys.argv` containing the path to `cochem_state.h5`.
2. **Matrix Read:** Reads the relaxed `xyz` coordinates from `/topos/mace_ensemble`.
3. **Hardware Profiling:** Calls `ram_estimator.py` to evaluate the active `NODE` memory limit.
4. **Execution Generation:** Calls `orca_writer.py` to generate the `.inp` files.
5. **Output Parsing:** Once ORCA finishes (handled asynchronously by NODE), BENCH parses the `.out` file, extracts the T1 diagnostic, and throws an `AssertionError` if $T1 > 0.02$.
6. **Serialization:** Writes the final CBS energy to `/base/thermo_limits/cbs_energy`.

## 4. PyTest Roadmap
- **Test 1 (`test_extrapolate.py`):** Provide known exact X=3 and X=4 energies from literature. Assert that `calculate_halkier` returns the exact literature CBS limit within $10^{-6}$ Hartrees.
- **Test 2 (`test_ram_estimator.py`):** Assert that a massive system (e.g., 2000 basis functions) forces `set_storage_strategy` to return `"RIJCOSX_DISK"`.
