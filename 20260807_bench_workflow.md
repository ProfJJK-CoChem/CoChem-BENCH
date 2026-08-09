# CoChem-BENCH: Execution Workflow (2026-08-07)

## Phase 1: Pre-Flight & Extrapolation Strategy
1. **Geometry Lock:** BENCH receives fully relaxed geometries from TOPOS (via r2SCAN-3c/MACE-OFF24m). It performs a strict standard-orientation serialization.
2. **RAM Prediction:** The RAM estimator predicts the amplitude storage requirements. If $>90\%$ of node memory is required, the algorithm switches to Cholesky Decomposition or disk-based RIJCOSX storage.
3. **CBS Dispatch:** BENCH dispatches independent Single Points to NODE: Hartree-Fock SCF energy, and the CCSD(T) correlation energy using X=3 and X=4 basis sets.

## Phase 2: Execution & Monitoring
1. **Dynamic TFLOPS:** The Jupyter UI generates a dynamic Plotly TFLOPS progress bar indicating the exact completion time of the canonical/DLPNO equations.
2. **T1 Diagnostic Guard:** Upon completion, the T1 diagnostic is checked. If $T1 > 0.02$, BENCH halts and flashes a red UI warning that the system is multireference, appending a disclaimer to SCRIBE.

## Phase 3: Reporting
1. **CBS Integration:** The Karton-Martin HF extrapolation and Halkier correlation extrapolations are summed.
2. **Provenance:** BENCH cryptographically locks the ORCA outputs and streams the final energies and specific BLAS/LAPACK metadata to `cochem_state.h5`.
