# CoChem-BENCH: Architectural Changes (2026-08-07)

## 1. Explicit Correlation (F12) CBS Limits
**Target File:** `bench_core/extrapolate.py`
**Required Architectural Change:**
- BENCH must calculate true Complete Basis Set (CBS) limits by enforcing F12 explicit correlation for CCSD(T) pathways (e.g., DLPNO-CCSD(T)-F12/cc-pVTZ-F12). The $1/X^3$ Halkier extrapolation must exclusively use cc-pVTZ and cc-pVQZ equivalents.

## 2. Core-Valence Integration
**Target File:** `bench_core/basis_manager.py`
**Required Architectural Change:**
- BENCH must mandate `aug-cc-pCVTZ` (Core-Valence) basis sets automatically if it detects down-stream calls from CoChem-SHIFT (NMR). Core electron correlation is non-negotiable for magnetic shielding accuracy.

## 3. Dynamic Hardware-Aware Fallbacks
**Target File:** `bench_core/ram_estimator.py`
**Required Architectural Change:**
- The engine must profile RAM via NODE. If canonical CCSD(T) threatens to hit swap memory on a massive node, BENCH must intercept the command and silently degrade to `TightPNO DLPNO-CCSD(T)` to save the job from crashing, appending the swap event to the SCRIBE methodology.

## 4. Academic Integrity Enforcements
**Target File:** `bench_core/verifier.py`
**Required Architectural Change:**
- Disable geometrical Counterpoise (gCP) on CBS limit calculations.
- Force Grimme D4 dispersion citations into SCRIBE.
