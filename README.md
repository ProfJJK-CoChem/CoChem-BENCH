# CoChem-BENCH: Complete Basis Set (CBS) Extrapolation & Ab Initio Thermochemistry

**Author/PI:** Dr. Joshua John Klaassen (Tennessee Tech University)  
**ORCiD:** [https://orcid.org/0009-0007-1506-4401](https://orcid.org/0009-0007-1506-4401)  
**GitHub Organization:** [https://github.com/ProfJJK-CoChem](https://github.com/ProfJJK-CoChem)  
**License:** Apache License 2.0

---

## 1. Module Overview

`CoChem-BENCH` (`cochem-bench` v0.2.0) provides high-accuracy *ab initio* Complete Basis Set (CBS) extrapolation, additive composite thermochemistry (core-valence $\Delta E_{\mathrm{CV}}$, scalar relativistic $\Delta E_{\mathrm{DKH2/X2C}}$, ZPE/thermal), relative energy ranking tables (`Eh`, `kcal/mol`, `kJ/mol`, `cm^-1` using exact CODATA constants), real ORCA 6.x output parsing, and multi-reference / convergence diagnostic auditing for CHEM 410 and CHEM 490 workflows.

## 2. Core Capabilities

- **Two-Point CBS Extrapolation (`cochem_bench.extrapolate`)**:
  - Exponential SCF extrapolation: $E_{\mathrm{SCF}}(X) = E_{\mathrm{SCF}}^\infty + A \exp(-\alpha \sqrt{X})$
  - Inverse-power correlation extrapolation: $E_{\mathrm{corr}}(X) = E_{\mathrm{corr}}^\infty + B X^{-\beta}$ (exact Helgaker $\beta = 3$, Neese-Valeev 2011 optimal $(\alpha, \beta)$, and F12 $\beta = 7$)
  - Verified parameter tables for `cc-pVnZ`, `aug-cc-pVnZ`, `def2`, and `ano-pVnZ`.
- **Additive Composite Thermochemistry & Relative Energies (`cochem_bench.thermochemistry`)**:
  - $\Delta E_{\mathrm{CV}}$, $\Delta E_{\mathrm{DKH2/X2C}}$, $\Delta \mathrm{ZPE}$, and thermal H/G corrections.
  - Strict integrity: missing corrections are explicitly tracked as `PENDING` (never silently substituted with `0.0`).
- **Real ORCA 6.x Parser (`cochem_bench.orca_parser`)**:
  - Parses `FINAL SINGLE POINT ENERGY`, SCF reference energy, correlation energy, `T1`/`D1` diagnostics, SCF/optimization convergence, `SlowConv` flags, and normal termination status directly from real ORCA `.out` files.
- **`CoChem-BASE` Provider (`cochem_bench.provider:bench_provider`)**:
  - Conforms to `cochem.module-provider/1` / `cochem.module-handoff/1` via the `cochem.modules` entry point `bench`.
  - Exposes operations: `cbs_extrapolate`, `relative_energies`, `parse_orca_energy`, `audit_diagnostics`.

## 3. Installation & CLI Usage

```bash
pip install --no-build-isolation -e .
cochem-bench info
cochem-bench parse-orca path/to/orca.out
cochem-bench cbs --e-scf-x -76.057133912 --e-scf-y -76.064825000 --e-corr-x -0.267422659 --e-corr-y -0.286950000 -x 3 -y 4 --basis-family cc-pVnZ
```

## 4. Documentation & User Guide

- **[User Guide (`docs/USER_GUIDE.md`)](docs/USER_GUIDE.md)** — Worked CLI and Python examples covering `cochem-bench parse-orca`, Helgaker/Neese–Valeev/F12 CBS extrapolation, composite thermochemistry & relative energy ranking, diagnostics auditing, and CHEM 490 professor-return parsing.

## 5. License

Licensed under the **Apache License, Version 2.0**. See [LICENSE](LICENSE) for full terms.
