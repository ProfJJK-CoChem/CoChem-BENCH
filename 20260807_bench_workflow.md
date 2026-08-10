# CoChem-BENCH: Execution Workflow & Working Set Benchmark Protocol (v4 Method Matrix §17)

**Status**: ARCHIVED & SUPERSEDED by v4 Method Matrix Standard (§17)  
**Effective Date**: 2026-08-09  

---

## 1. Archive Notice
The legacy 4-stage execution workflow (2026-08-07) has been archived. All benchmark validation, error quantification, and reference dataset verifications in `CoChem-BENCH` are governed by the **v4 6-System Working Set** (§17) and **Split-Conformal Prediction Protocol** (§17.5).

---

## 2. Mandatory v4 6-System Working Set Protocol (§17)

Every benchmark campaign must evaluate and report performance against the 6 reference systems:

| # | System | Measured Rotational Constants (MHz) | Primary Test Objective |
|---|---|---|---|
| 1 | **Ar–ketene ($\text{H}_2\text{CCO}\cdots\text{Ar}$)**, $A_1$ state | $A = 10447.9248$, $B = 1918.0138$, $C = 1606.7642$ | Rare gas + tunnelling + measured dipole components ($\mu_a, \mu_b$). |
| 2 | **Ar–oxazole ($\text{C}_3\text{H}_3\text{ON}\cdot\text{Ar}$)** | $A = 5012.89486$, $B = 1398.42815$, $C = 1388.95284$ | $^{14}\text{N}$ quadrupole tensor ($\chi_{aa}, \chi_{bb}, \chi_{cc}$) + near-prolate asymmetry ($\kappa = -0.99$). |
| 3 | **H₂CO···H³⁵Cl ($\text{H}_2\text{CO}\cdots\text{H}^{35}\text{Cl}$)** | $B = 2687.856$, $C = 2527.412$, $A$ fixed at 42 GHz | Hydrogen bond + large Cl quadrupole + prior-disclosure rule testing. |
| 4 | **Water dimer $(\text{H}_2\text{O})_2 / (\text{D}_2\text{O})_2$** | $(\text{H}_2\text{O})_2$: $A = 227580.432$; $(\text{D}_2\text{O})_2$: $A = 120327.492$ | Tunnelling stress test (9 orders of magnitude splittings) + H/D isotope ratios. |
| 5 | **NH₃···HCOOH ($\text{NH}_3\cdots\text{HCOOH}$)** | Experimental $V_3 = 195.18$ $\text{cm}^{-1}$ | Internal rotation barrier $V_3$ and A–E splitting test. |
| 6 | **C₆H₆–HCN ($\text{C}_6\text{H}_6\cdots\text{HCN}$) vs Ar₃–HCN ($\text{Ar}_3\cdots\text{HCN}$)** | $R_{cm} = 3.96$ Å vs $3.47$ Å | Matched topology pair binding discrimination. |

---

## 3. Key Invariants & Provenance Discipline (§12.5)
1. **Convergence Block**: Geometry optimizations must enforce the 5-threshold `%geom` block (`TolE 1e-7`, `TolRMSG 3e-6`, `TolMaxG 1e-5`, `TolRMSD 5e-5`, `TolMaxD 1e-4`).
2. **BSSE Rules**: CP-OPT is prohibited on CBS extrapolation rows and enforced on non-augmented triple-zeta sets.
3. **Frozen-Core Auditor**: Rows claiming $\le 0.5\%$ error in $B_e$ must evaluate core-valence correlation (`cc-pCVTZ`).
4. **Provenance Tagging**: Every value carries `[M]` (measured), `[D]` (derived), or `[E]` (estimated) tags. Rule 7 compliance is strictly enforced.
