# CoChem-BENCH User Guide

This guide provides worked CLI and Python examples for `CoChem-BENCH` (`cochem-bench` v0.2.0):

1. [Parsing ORCA Outputs (`cochem-bench parse-orca`)](#1-parsing-orca-outputs-cochem-bench-parse-orca)
2. [Two-Point CBS Helgaker / Neese–Valeev / F12 Extrapolation](#2-two-point-cbs-helgaker--neesevaleev--f12-extrapolation)
3. [Additive Composite Thermochemistry & Relative Energies](#3-additive-composite-thermochemistry--relative-energies)
4. [Diagnostics & Integrity Audits (`cochem-bench audit`)](#4-diagnostics--integrity-audits-cochem-bench-audit)
5. [CHEM 490 Professor-Return Parsing Workflow](#5-chem-490-professor-return-parsing-workflow)

---

## 1. Parsing ORCA Outputs (`cochem-bench parse-orca`)

`cochem_bench.orca_parser` extracts verified energies, convergence state, `T1`/`D1` coupled-cluster diagnostics, and `SlowConv` flags directly from real ORCA 6.x `.out` logs.

### CLI Example

```bash
# Inspect module metadata and CODATA conversion constants
cochem-bench info

# Parse a converged ORCA 6.1.1 output file to structured JSON
cochem-bench parse-orca path/to/h2o_ccsd_pvtz.out
```

### Python Example

```python
from pathlib import Path
from cochem_bench.orca_parser import parse_orca_energy

res = parse_orca_energy(Path("path/to/h2o_ccsd_pvtz.out"))
print("ORCA version:", res.orca_version)
print("Final SPE (Eh):", res.final_single_point_energy_eh)
print("SCF Energy (Eh):", res.scf_energy_eh)
print("Corr Energy (Eh):", res.correlation_energy_eh)
print("T1 diagnostic:", res.t1_diagnostic)
print("Converged:", res.scf_converged, res.opt_converged, res.terminated_normally)
```

---

## 2. Two-Point CBS Helgaker / Neese–Valeev / F12 Extrapolation

`cochem_bench.extrapolate` implements two-point Complete Basis Set (CBS) extrapolations:

- **SCF (exponential):** $E_{\mathrm{SCF}}(X) = E_{\mathrm{SCF}}^\infty + A \exp(-\alpha \sqrt{X})$
- **Correlation (inverse-power):** $E_{\mathrm{corr}}(X) = E_{\mathrm{corr}}^\infty + B X^{-\beta}$
  - `scheme="helgaker"`: exact inverse-cubic $\beta = 3.0$
  - `scheme="neese_valeev"`: Neese–Valeev (2011) basis-pair optimal $(\alpha, \beta)$
  - `scheme="f12"`: explicitly correlated $\beta = 7.0$
- Supported basis families: `cc-pVnZ`, `aug-cc-pVnZ`, `def2`, `ano-pVnZ`.

### CLI Example

```bash
cochem-bench cbs \
  --e-scf-x -76.057133912 \
  --e-scf-y -76.064825000 \
  --e-corr-x -0.267422659 \
  --e-corr-y -0.286950000 \
  -x 3 -y 4 \
  --basis-family cc-pVnZ \
  --scheme helgaker
```

*(You may also pass `--e-tot-x` and `--e-tot-y` instead of `--e-corr-x`/`--e-corr-y`; `CoChem-BENCH` automatically computes $E_{\mathrm{corr}} = E_{\mathrm{tot}} - E_{\mathrm{SCF}}$.)*

### Python Example

```python
from cochem_bench.extrapolate import cbs_extrapolate

cbs = cbs_extrapolate(
    e_scf_x=-76.057170,
    e_scf_y=-76.064825,
    e_corr_x=-0.275031,
    e_corr_y=-0.290322,
    x=3,
    y=4,
    basis_family="cc-pVnZ",
    scheme="helgaker",
)
print(f"E_SCF(CBS)  = {cbs.e_scf_cbs_eh:.9f} Eh (alpha={cbs.alpha})")
print(f"E_corr(CBS) = {cbs.e_corr_cbs_eh:.9f} Eh (beta={cbs.beta})")
print(f"E_tot(CBS)  = {cbs.e_total_cbs_eh:.9f} Eh  [tag: {cbs.provenance_tag}]")
```

---

## 3. Additive Composite Thermochemistry & Relative Energies

`cochem_bench.thermochemistry` assembles additive composite energies ($E_{\mathrm{CBS}} + \Delta E_{\mathrm{CV}} + \Delta E_{\mathrm{DKH2/X2C}} + \Delta\mathrm{ZPE}$) and ranks conformers/isomers in `Eh`, `kcal/mol`, `kJ/mol`, and `cm^-1` using exact CODATA constants.

**Strict Integrity Policy:** Missing corrections are **never** silently substituted with `0.0`. Any species missing a required correction is returned with `status="PENDING"`, `pending_corrections=[...]`, and `None` relative energies.

### CLI Example

Prepare `conformers.json`:

```json
{
  "species": [
    {
      "name": "S1",
      "e_cbs_eh": -344.364232,
      "delta_cv_eh": -0.020100,
      "delta_zpe_eh": 0.114500,
      "required_corrections": ["delta_cv", "delta_zpe"]
    },
    {
      "name": "S2",
      "e_cbs_eh": -344.360000,
      "delta_cv_eh": -0.019950,
      "delta_zpe_eh": 0.114300,
      "required_corrections": ["delta_cv", "delta_zpe"]
    },
    {
      "name": "S3_awaiting_cv",
      "e_cbs_eh": -344.358500,
      "delta_zpe_eh": 0.114100,
      "required_corrections": ["delta_cv", "delta_zpe"]
    }
  ]
}
```

Run:

```bash
cochem-bench rel-energy conformers.json --reference S1
```

### Python Example

```python
from cochem_bench.thermochemistry import compute_composite_energy, compute_relative_energies

s1 = compute_composite_energy(
    name="S1",
    e_cbs_eh=-344.364232,
    e_ae_cv_eh=-344.520100,
    e_fc_cv_eh=-344.500000,
    delta_zpe_eh=0.114500,
    required_corrections=("delta_cv", "delta_zpe"),
)
s2 = compute_composite_energy(
    name="S2",
    e_cbs_eh=-344.360000,
    delta_cv_eh=-0.019950,
    delta_zpe_eh=0.114300,
    required_corrections=("delta_cv", "delta_zpe"),
)
s3 = compute_composite_energy(
    name="S3_pending",
    e_cbs_eh=-344.358500,
    delta_zpe_eh=0.114100,
    required_corrections=("delta_cv", "delta_zpe"),
)

table = compute_relative_energies([s1, s2, s3], reference="S1")
for entry in table.entries:
    print(entry.rank, entry.name, entry.status, entry.delta_e_kcal_mol, entry.pending_corrections)
```

---

## 4. Diagnostics & Integrity Audits (`cochem-bench audit`)

`cochem_bench.verifier.audit_diagnostics` checks:

- **Termination & convergence:** `ORCA TERMINATED NORMALLY`, `SCF CONVERGED`, geometry optimization convergence, and absence of `SlowConv`.
- **Multi-reference character:** Closed-shell thresholds `T1 <= 0.020`, `D1 <= 0.050`; open-shell (`spin_multiplicity > 1`) thresholds `T1 <= 0.045`, `D1 <= 0.070`. Missing diagnostics (e.g., canonical CCSD without D1) are reported honestly as `"NOT_REPORTED"`.
- **Frozen-core bias gate:** Flags accuracy claims `<= 0.5%` if core-valence correlation was not evaluated (`-0.81%` systematic rotational-constant bias).
- **Spin contamination:** Flags open-shell deviations `|S^2_actual - S^2_expected| / S^2_expected > 10%`.

### CLI Example

```bash
# Closed-shell diagnostic audit
cochem-bench audit path/to/h2o_ccsd_pvtz.out --multiplicity 1

# Audit with a tight accuracy claim and core-valence check
cochem-bench audit path/to/h2o_ccsd_pvtz.out --multiplicity 1 --claimed-error-pct 0.3 --has-core-valence
```

### Python Example

```python
from cochem_bench.verifier import audit_diagnostics

report = audit_diagnostics(
    "path/to/h2o_ccsd_pvtz.out",
    spin_multiplicity=1,
    claimed_error_pct=0.3,
    has_core_valence=False,
)
print("Audit status:", report["status"], "| Findings:", report["findings"])
```

---

## 5. CHEM 490 Professor-Return Parsing Workflow

In `CHEM 490` directed research, students stage lightweight inputs on their own machines and receive converged heavy ORCA 6.x `.out` files back from the professor's workstation or cloud runner. `CoChem-BENCH` turns a directory of returned `.out` files into a provenance-audited CBS & relative-energy package in four steps:

```python
from pathlib import Path
from cochem_bench.extrapolate import cbs_extrapolate
from cochem_bench.orca_parser import parse_orca_energy
from cochem_bench.thermochemistry import compute_composite_energy, compute_relative_energies
from cochem_bench.verifier import audit_diagnostics, export_provenance_json

return_dir = Path("chem490_professor_return")
conformer_specs = [
    ("Conformer_A", "confA_tz.out", "confA_qz.out"),
    ("Conformer_B", "confB_tz.out", "confB_qz.out"),
]

composite_results = []
for name, tz_file, qz_file in conformer_specs:
    tz_path = return_dir / tz_file
    qz_path = return_dir / qz_file

    # 1. Parse returned ORCA outputs directly from disk
    tz = parse_orca_energy(tz_path)
    qz = parse_orca_energy(qz_path)

    # 2. Gate every returned output with convergence & T1/D1 audits
    tz_audit = audit_diagnostics(tz, spin_multiplicity=1)
    qz_audit = audit_diagnostics(qz, spin_multiplicity=1)
    if not (tz_audit["passed"] and qz_audit["passed"]):
        raise RuntimeError(f"{name} failed diagnostic audit: {tz_audit['findings'] + qz_audit['findings']}")

    # 3. Extrapolate T/Q energies to CBS (Helgaker beta=3)
    cbs = cbs_extrapolate(
        e_scf_x=tz.scf_energy_eh,
        e_scf_y=qz.scf_energy_eh,
        e_corr_x=tz.correlation_energy_eh,
        e_corr_y=qz.correlation_energy_eh,
        x=3,
        y=4,
        basis_family="cc-pVnZ",
        scheme="helgaker",
    )

    # 4. Assemble composite energy (uncomputed CV/DKH2/ZPE remain explicitly PENDING)
    composite_results.append(
        compute_composite_energy(name=name, e_cbs_eh=cbs.e_total_cbs_eh)
    )

# 5. Rank conformers and export Standing-Rule-7 provenance JSON
table = compute_relative_energies(composite_results)
export_provenance_json(
    input_dict={"return_dir": str(return_dir), "basis_pair": "cc-pVTZ/cc-pVQZ"},
    results_dict=table.model_dump(mode="json"),
    output_path=return_dir / "chem490_bench_provenance.json",
    orca_output_path=return_dir / conformer_specs[0][2],
)
```

For `CoChem-BASE` pipeline orchestration, the same operations (`cbs_extrapolate`, `relative_energies`, `parse_orca_energy`, `audit_diagnostics`) can be dispatched via `cochem_bench.provider:execute_handoff(handoff_path, output_directory)`.
