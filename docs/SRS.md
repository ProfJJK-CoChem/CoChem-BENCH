---
schema: cochem-srs/5.0
repo: CoChem-BENCH
status: ACTIVE_CORE (CBS/F12 and ORCA parser with (T) triples recovery operational; Boltzmann ensemble and fragment escalation pending)
pipeline_baseline: 5.0.0
source_commit: c3283eb0f23c328ca862177896766b93aa0ef237
decision_ref: swarm-1126-8efd
---

# Software Requirements Specification: CoChem-BENCH (Pipeline 5.0.0)

## 1. Purpose, Scope & System Classification

`CoChem-BENCH` (`ProfJJK-CoChem/CoChem-BENCH`, Python distribution `cochem-bench` v0.2.0) is classified as an **`ACTIVE_CORE`** (`wave1_core_9`) repository in the CoChem Pipeline 5.0.0 ecosystem. It provides deterministic Complete Basis Set (CBS) and F12 extrapolation algebra, composite thermochemical auditing (core-valence, scalar relativistic, ZPE, BSSE), multireference diagnostics (`T1`/`D1`), and real ORCA 6.x output parsing.

The installable package resides under `cochem_bench/` (22 production `.py` files, 0 syntax errors, 0 `NotImplementedError` stubs, 0 broken imports) with CLI entrypoint `cochem-bench = "cochem_bench.cli:main"` exposing 5 subcommands (`info`, `cbs`, `parse-orca`, `audit`, `rel-energy`). Core modules include `cochem_bench/extrapolate.py`, `cochem_bench/orca_parser.py`, `cochem_bench/orca_writer.py`, `cochem_bench/orca_runner.py`, `cochem_bench/basis_manager.py`, `cochem_bench/thermochemistry.py`, `cochem_bench/verifier.py`, `cochem_bench/ram_estimator.py`, `cochem_bench/dispatcher.py`, and `cochem_bench/provider.py`.

## 2. Pipeline 5.0.0 Architecture Alignment & Resolved Conflicts

`CoChem-BENCH` aligns with signed Portal ruling `swarm-1126-8efd` and resolves legacy v4.x conflicts:
- **`CONTRA-01` (Compute Topology & `HC-1`):** Heavy coupled-cluster and composite ORCA runs execute exclusively on local WSL2 (`host-guard`), GCP Spot CPU (`c2d`/`c4d`/`c3d`), or Colab GPU via the Pipeline 5.0.0 runner, never on GitHub Actions compute runners. `.github/workflows/zero_mock_ci.yml` executes zero-cron lightweight validation only.
- **`CONTRA-05` (`DLPNOExclusionFirewall` & Method Matrix §4.4):** Under Pipeline 5.0.0, `DLPNOExclusionFirewall` prohibits new automated `DLPNO-CCSD(T)` submissions in favor of canonical `CCSD(T)` / `CCSD(T)-F12b` and `revDSD-PBEP86-D4`. Legacy 2026-08-07 directives calling for automatic degradation from canonical `CCSD(T)` to `TightPNO DLPNO-CCSD(T)` are superseded; any memory-limited fallback requires explicit human/Portal approval rather than silent DLPNO downgrade.
- **`CONTRA-06` (Unified Contracts):** Cross-repo benchmark outputs target `cochem-job-handoff/1`, `cochem-quantities/1` (`OWN-117`), and `cochem-portal-dossier/1.1` DRAFT schemas in `CoChem-BASE`.
- **Research Notebook Provenance (`NB-1`):** Commit `e856bc29` explicitly labeled hand-typed TZ/QZ inputs in `notebooks/CoChem-BENCH_Research.ipynb` as illustrative algebraic demonstrations (not experimental/DFT reference energies) with zero stored outputs.

## 3. Authoritative Scientific Engines, Methods & Method Matrix Constraints

- **CBS & F12 Extrapolation (`cochem_bench/extrapolate.py`, `cochem_bench/constants.py`):** Implements two-point exponential SCF extrapolation and inverse-power Helgaker / Neese-Valeev / F12 correlation extrapolation with exact component accounting.
- **ORCA 6.x Input Writer, Runner & Parser (`cochem_bench/orca_writer.py`, `cochem_bench/orca_runner.py`, `cochem_bench/orca_parser.py`, `cochem_bench/basis_manager.py`):** Generates Boys-Bernardi BSSE ghost-atom inputs, core-valence (`aug-cc-pCVTZ`) and scalar relativistic (`DKH2`/`X2C`/`ZORA`) directives, enforces `DefGrid3` (disabling `gCP` on CBS tiers), and parses authentic ORCA 6.x logs including `(T)` perturbative triples recovery (`Final correlation energy` + `triples_correction_eh`, verified at `-0.007648001 Eh` on ORCA 6.1.1 `CCSD(T)/cc-pVTZ`, commit `c3283eb` / P1-02).
- **Diagnostics & Thermochemical Integrity (`cochem_bench/verifier.py`, `cochem_bench/thermochemistry.py`):** Implements `verify_multireference_gate`, `audit_diagnostics` (`T1`/`D1` thresholds), and `audit_frozen_core_bias`. Missing composite terms (core-valence, relativistic, ZPE) are emitted honestly as `PENDING` rather than fabricated or silently zeroed.
- **Method Matrix Constraints:** `cochem_bench/orca_writer.py` (line 384) currently emits `! TightOpt` shorthand and is tracked in Section 6 for migration to an explicit 5-threshold `%geom` block per Method Matrix §4.4.

## 4. Cross-Repo Interfaces, Contracts & Canonical Schemas

- **CoChem-BASE Module Provider:** `pyproject.toml` registers `[project.entry-points."cochem.modules"]` `bench = "cochem_bench.provider:bench_provider"` implemented in `cochem_bench/provider.py` and validated by `cochem_bench/core/payload_validation.py`.
- **Canonical Pipeline 5.0.0 Schemas:** Benchmark energy tables and uncertainty budgets align with `cochem-quantities/1` (`QUANTITIES.json`, `OWN-117`) and `cochem-job-handoff/1` (`HANDOFF.json`) DRAFT schemas hosted in `CoChem-BASE`.
- **Ecosystem Resource Queries:** `cochem_bench/ram_estimator.py` exposes memory estimation (`node_max_gb` parameter) for upstream orchestration without hard-depending on quarantined `CoChem-NODE`.

## 5. Verification, Zero-Mock CI Gates & Acceptance Criteria

- **Zero-Mock CI Gate:** `tests/test_zero_mock_and_no_sim_gate.py` enforces strict AST scanning of `.py` and `.ipynb` files (`.anti_spoof_amnesty.json` strictly `[]`), passing `zero_mock_gate` on PR #3 via `.github/workflows/zero_mock_ci.yml`.
- **Unit & Scientific Tests:** 7 domain test modules verify CBS algebra, ORCA I/O, diagnostics, and memory estimation: `tests/test_bench_core.py`, `tests/test_cbs_scientific.py`, `tests/test_extrapolate.py`, `tests/test_m3_verification.py`, `tests/test_orca_parser_and_provider.py`, `tests/test_orca_writer_and_runner.py`, and `tests/test_ram_estimator.py`.
- **Branch Protection:** `main` is protected requiring authentic `blind-s3-p10` parity before PR #3 merges.

## 6. Verified Implementation Gaps, Known Holds & Non-Goals

The following verified gaps from the 2026-11-26 Opus-verified audit remain open:
1. **Macroscopic Boltzmann Ensemble & Fragment Escalation Absent:** Macroscopic Boltzmann ensemble synthesis and Weak/Strong Fragment-Based Escalation Protocols for multi-molecular clusters are not yet implemented anywhere in `cochem_bench/`.
2. **Static RAM Default & Unwired SCRIBE Logging:** `cochem_bench/ram_estimator.py` (line 80) defaults statically to `node_max_gb=64.0` with no live hardware/orchestrator integration (consistent with `CoChem-NODE` quarantine), and automated `CoChem-SCRIBE` methodology citation logging is not wired.
3. **Method Matrix `! TightOpt` Drift:** `cochem_bench/orca_writer.py` (line 384) emits `! TightOpt` shorthand and requires replacement with explicit 5-threshold `%geom` convergence criteria.
4. **Non-Goals:** `CoChem-BENCH` does not run conformer enumeration (`CoChem-TOPOS`), rovibrational VPT2 spectroscopy (`CoChem-TORQ`), or automated DLPNO fallback submissions (prohibited by `DLPNOExclusionFirewall`).

## 7. Legacy Specification Disposition (Non-Normative Archive)

Prior to `swarm-1126-8efd`, `CoChem-BENCH` had no formal SRS. The following historical repository files are non-normative archives superseded by `docs/SRS.md`:
- `20260807_bench_architectural_changes.md` (references retired `bench_core/*` paths replaced by `cochem_bench/*` and retired auto-DLPNO degradation superseded by `DLPNOExclusionFirewall`).
- `20260807_bench_workflow.md` and `workflow.md` (historical pre-v5 workflow notes).
- `docs/USER_GUIDE.md` remains the operational user guide subordinate to this normative `docs/SRS.md`.
