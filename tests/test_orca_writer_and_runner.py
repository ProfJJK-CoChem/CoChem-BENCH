#!/usr/bin/env python3
"""Zero-mock tests for CoChem-BENCH ORCA 6.1.1 writer and authentic ORCA 6 runner."""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

import pytest

from cochem_bench.basis_manager import (
    PERIODIC_TABLE_ATOMIC_NUMBERS,
    get_atomic_number,
    is_3d_transition_metal,
    is_heavy_element,
    requires_scalar_relativistic,
)
from cochem_bench.orca_runner import (
    OrcaExecutionError,
    resolve_orca_binary,
    run_orca_calculation,
    run_orca_input,
)
from cochem_bench.orca_writer import (
    _resolve_relativistic_basis,
    generate_canonical_ccsd_t,
    generate_counterpoise_input,
    generate_dlpno_ccsd_f12,
    generate_dlpno_ccsd_f12_with_metadata,
    resolve_relativistic_basis,
)


def _header_tokens(deck: str) -> list[str]:
    header = [line for line in deck.splitlines() if line.startswith("!")][0]
    return header.split()


def test_complete_periodic_table_coverage_through_oganesson() -> None:
    """Periodic table lookup must cover every IUPAC element Z=1..118 without truncation at Xe."""
    assert len(PERIODIC_TABLE_ATOMIC_NUMBERS) == 118
    assert get_atomic_number("H") == 1
    assert get_atomic_number("Sc") == 21
    assert get_atomic_number("Zn") == 30
    assert get_atomic_number("Br") == 35
    assert get_atomic_number("Kr") == 36
    assert get_atomic_number("Xe") == 54
    assert get_atomic_number("Cs") == 55
    assert get_atomic_number("Au") == 79
    assert get_atomic_number("Rn") == 86
    assert get_atomic_number("U") == 92
    assert get_atomic_number("Og") == 118
    # Ghost atom colon suffix normalized cleanly
    assert get_atomic_number("Au:") == 79

    assert is_3d_transition_metal("Fe") is True
    assert is_3d_transition_metal("Cu") is True
    assert is_3d_transition_metal("Ca") is False
    assert is_3d_transition_metal("Ga") is False

    assert is_heavy_element("Br") is False
    assert is_heavy_element("Kr") is True
    assert is_heavy_element("Cs") is True
    assert is_heavy_element("Au") is True
    assert is_heavy_element("U") is True
    assert is_heavy_element("Og") is True

    assert requires_scalar_relativistic("O") is False
    assert requires_scalar_relativistic("Fe") is True
    assert requires_scalar_relativistic("Br") is True
    assert requires_scalar_relativistic("Kr") is True
    assert requires_scalar_relativistic("Og") is True


def test_orca6_writer_eliminates_invalid_bare_cabs_and_relativistic_keywords() -> None:
    """No generated deck may contain bare CABS or invalid 'Relativistic DKH2'/'Relativistic ZORA'."""
    h2 = [("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.74)]
    hbr = [("Br", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.414)]
    auh = [("Au", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.524)]
    feh = [("Fe", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.60)]
    og2 = [("Og", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.80)]

    # 1. Light F12 SP
    deck_f12 = generate_dlpno_ccsd_f12(h2, basis="cc-pVDZ-F12", is_opt=False)
    tokens_f12 = _header_tokens(deck_f12)
    assert "CABS" not in tokens_f12
    assert "cc-pVDZ-F12-CABS" in tokens_f12
    assert "TightPNO" in tokens_f12
    assert "AutoAux" in tokens_f12
    assert "NumGrad" not in tokens_f12

    # 2. Light aug-cc-pVTZ F12 Opt -> includes TightPNO, cc-pVTZ-F12-CABS, aug-cc-pVTZ/C, NumGrad
    deck_opt = generate_dlpno_ccsd_f12(h2, basis="aug-cc-pVTZ", is_opt=True)
    tokens_opt = _header_tokens(deck_opt)
    assert "CABS" not in tokens_opt
    assert "cc-pVTZ-F12-CABS" in tokens_opt
    assert "aug-cc-pVTZ/C" in tokens_opt
    assert "TightPNO" in tokens_opt
    assert "Opt" in tokens_opt
    assert "NumGrad" in tokens_opt

    # 3. Heavy-element DKH2 (HBr, AuH, OgH) and 3d TM (FeH)
    for coords, basis in [(hbr, "def2-SVP"), (auh, "def2-TZVP"), (feh, "def2-TZVP"), (og2, "def2-TZVP")]:
        deck_rel = generate_dlpno_ccsd_f12(coords, basis=basis, rel_mode="auto")
        tokens_rel = _header_tokens(deck_rel)
        assert "Relativistic" not in deck_rel
        assert "CABS" not in tokens_rel
        assert "DKH2" in tokens_rel
        assert "TightPNO" in tokens_rel
        assert "DefGrid3" in tokens_rel
        assert "TightSCF" in tokens_rel
        assert 'CABS "AutoAux"' in deck_rel

    # 4. ZORA explicit mode on HBr
    deck_zora = generate_dlpno_ccsd_f12(hbr, basis="def2-SVP", rel_mode="ZORA")
    tokens_zora = _header_tokens(deck_zora)
    assert "Relativistic" not in deck_zora
    assert "ZORA" in tokens_zora
    assert "ZORA-def2-SVP" in tokens_zora
    assert 'CABS "AutoAux"' in deck_zora

    # 5. Canonical CCSD(T) SP and Opt
    deck_cc_sp = generate_canonical_ccsd_t(h2, basis="cc-pVDZ", is_opt=False)
    tokens_cc_sp = _header_tokens(deck_cc_sp)
    assert "CCSD(T)" in tokens_cc_sp
    assert "TightPNO" not in tokens_cc_sp
    assert "CABS" not in deck_cc_sp
    assert "NumGrad" not in tokens_cc_sp

    deck_cc_opt = generate_canonical_ccsd_t(h2, basis="cc-pVDZ", is_opt=True)
    tokens_cc_opt = _header_tokens(deck_cc_opt)
    assert "CCSD(T)" in tokens_cc_opt
    assert "Opt" in tokens_cc_opt
    assert "NumGrad" in tokens_cc_opt


def test_orca_runner_fails_closed_when_binary_missing(tmp_path: Path) -> None:
    """Missing ORCA executable must fail closed with FileNotFoundError (zero synthesis)."""
    missing = tmp_path / "does_not_exist_orca"
    with pytest.raises(FileNotFoundError):
        resolve_orca_binary(missing)

    prev_orca_bin = os.environ.get("ORCA_BIN")
    try:
        os.environ["ORCA_BIN"] = str(missing)
        with pytest.raises(FileNotFoundError):
            resolve_orca_binary()
    finally:
        if prev_orca_bin is None:
            os.environ.pop("ORCA_BIN", None)
        else:
            os.environ["ORCA_BIN"] = prev_orca_bin

    deck = generate_dlpno_ccsd_f12([("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.74)], basis="cc-pVDZ-F12")
    with pytest.raises(FileNotFoundError):
        run_orca_input(deck, tmp_path / "missing_run", orca_bin=missing)


def test_orca_runner_rejects_non_orca6_banner_and_nonzero_rc(tmp_path: Path) -> None:
    """Runner must enforce rc=0, normal termination, and authentic Program Version 6.x banner."""
    # Real executable subprocess (python) returning non-zero exit code
    if os.name == "posix":
        bad_rc = tmp_path / "bad_rc.sh"
        bad_rc.write_text("#!/bin/sh\necho 'failed'\nexit 4\n", encoding="utf-8")
        bad_rc.chmod(0o755)
        bad_banner = tmp_path / "bad_banner.sh"
        bad_banner.write_text(
            "#!/bin/sh\necho 'Program Version 5.0.4'\necho 'FINAL SINGLE POINT ENERGY -1.0'\n"
            "echo '****ORCA TERMINATED NORMALLY****'\nexit 0\n",
            encoding="utf-8",
        )
        bad_banner.chmod(0o755)
    else:
        bad_rc = tmp_path / "bad_rc.cmd"
        bad_rc.write_text("@echo off\r\necho failed\r\nexit /b 4\r\n", encoding="utf-8")
        bad_banner = tmp_path / "bad_banner.cmd"
        bad_banner.write_text(
            "@echo off\r\necho Program Version 5.0.4\r\necho FINAL SINGLE POINT ENERGY -1.0\r\n"
            "echo ****ORCA TERMINATED NORMALLY****\r\nexit /b 0\r\n",
            encoding="utf-8",
        )

    deck = generate_dlpno_ccsd_f12([("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.74)], basis="cc-pVDZ-F12")
    with pytest.raises(OrcaExecutionError, match="rc=4"):
        run_orca_input(deck, tmp_path / "run_bad_rc", orca_bin=bad_rc)
    with pytest.raises(OrcaExecutionError, match="Program Version 6"):
        run_orca_input(deck, tmp_path / "run_bad_banner", orca_bin=bad_banner)


def _live_orca_binary_or_skip() -> Path:
    try:
        return resolve_orca_binary()
    except FileNotFoundError:
        pytest.skip("Live ORCA 6 binary not installed on this host")


def test_live_orca6_h2_dlpno_f12_sp_and_opt_and_canonical_ccsdt(tmp_path: Path) -> None:
    """Live ORCA 6.1.1 execution of BENCH DLPNO-CCSD(T)-F12 SP, Opt+NumGrad, and canonical CCSD(T)."""
    orca_bin = _live_orca_binary_or_skip()
    h2_sp = [("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.74)]
    h2_opt = [("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.76)]

    # 1. DLPNO-CCSD(T)-F12 SP
    res_sp = run_orca_calculation(
        h2_sp,
        tmp_path / "h2_f12_sp",
        basis="cc-pVDZ-F12",
        is_opt=False,
        job_name="h2_f12_sp",
        orca_bin=orca_bin,
    )
    assert res_sp.returncode == 0
    assert res_sp.orca_version.startswith("6.")
    assert -1.20 < res_sp.final_single_point_energy_eh < -1.15
    assert res_sp.t1_diagnostic is not None and 0.0 < res_sp.t1_diagnostic < 0.02
    assert res_sp.audit["passed"] is True
    assert res_sp.provenance["orca_version"].startswith("ORCA 6.")
    assert res_sp.output_sha256 == hashlib.sha256(Path(res_sp.output_path).read_bytes()).hexdigest()

    # 2. DLPNO-CCSD(T)-F12 Opt with NumGrad
    res_opt = run_orca_calculation(
        h2_opt,
        tmp_path / "h2_f12_opt",
        basis="cc-pVDZ-F12",
        is_opt=True,
        in_hess="XTB2",
        job_name="h2_f12_opt",
        orca_bin=orca_bin,
    )
    assert res_opt.returncode == 0
    assert res_opt.parsed.opt_converged is True
    assert res_opt.t1_diagnostic is not None and 0.0 < res_opt.t1_diagnostic < 0.02
    assert res_opt.final_single_point_energy_eh <= res_sp.final_single_point_energy_eh + 1e-4

    # 3. Canonical CCSD(T) Opt with NumGrad
    res_cc = run_orca_input(
        generate_canonical_ccsd_t(h2_opt, basis="cc-pVDZ", is_opt=True, in_hess="Lindh"),
        tmp_path / "h2_canon_opt",
        job_name="h2_canon_opt",
        orca_bin=orca_bin,
        require_opt_converged=True,
    )
    assert res_cc.returncode == 0
    assert res_cc.parsed.opt_converged is True
    assert res_cc.t1_diagnostic is not None and 0.0 < res_cc.t1_diagnostic < 0.02
    assert -1.18 < res_cc.final_single_point_energy_eh < -1.14


def test_live_orca6_hbr_heavy_element_dkh2_and_zora(tmp_path: Path) -> None:
    """Live ORCA 6.1.1 execution of BENCH heavy-element scalar-relativistic DKH2 & ZORA decks on HBr."""
    orca_bin = _live_orca_binary_or_skip()
    hbr = [("Br", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.414)]

    # DKH2 (auto-selected on HBr)
    res_dkh2 = run_orca_calculation(
        hbr,
        tmp_path / "hbr_dkh2",
        basis="def2-SVP",
        rel_mode="auto",
        is_opt=False,
        job_name="hbr_dkh2",
        orca_bin=orca_bin,
    )
    assert res_dkh2.returncode == 0
    assert res_dkh2.orca_version.startswith("6.")
    assert res_dkh2.final_single_point_energy_eh < -2600.0
    assert res_dkh2.t1_diagnostic is not None and 0.0 < res_dkh2.t1_diagnostic < 0.02
    assert res_dkh2.audit["passed"] is True

    # ZORA (explicit on HBr)
    res_zora = run_orca_calculation(
        hbr,
        tmp_path / "hbr_zora",
        basis="def2-SVP",
        rel_mode="ZORA",
        is_opt=False,
        job_name="hbr_zora",
        orca_bin=orca_bin,
    )
    assert res_zora.returncode == 0
    assert res_zora.orca_version.startswith("6.")
    assert res_zora.final_single_point_energy_eh < -2600.0
    assert res_zora.t1_diagnostic is not None and 0.0 < res_zora.t1_diagnostic < 0.02
    assert res_zora.audit["passed"] is True
    assert res_zora.provenance["input"]["effective_basis"] == "ZORA-def2-SVP"
    assert res_zora.provenance["input"]["basis_swapped"] is True


def test_relativistic_basis_rejects_unrecognized_bases_and_mismatches() -> None:
    """M3 regression: unrecognized bases like 6-31G* must fail closed under DKH2/ZORA instead of inventing 6-31G*-DK."""
    hbr = [("Br", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.414)]
    auh = [("Au", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.524)]

    for bad_basis in ("6-31G*", "6-311++G(d,p)", "STO-3G", "LANL2DZ", "6-31G*-DK"):
        for rel_mode in ("DKH2", "ZORA", "auto"):
            with pytest.raises(ValueError, match="Unsupported basis"):
                generate_dlpno_ccsd_f12(hbr, basis=bad_basis, rel_mode=rel_mode)
            with pytest.raises(ValueError, match="Unsupported basis"):
                generate_dlpno_ccsd_f12(auh, basis=bad_basis, rel_mode=rel_mode)
        with pytest.raises(ValueError, match="Unsupported basis"):
            _resolve_relativistic_basis(bad_basis, "DKH2", [35, 1])
        with pytest.raises(ValueError, match="Unsupported basis"):
            _resolve_relativistic_basis(bad_basis, "ZORA", [35, 1])

    # Cross-Hamiltonian mismatches must also fail closed
    with pytest.raises(ValueError, match="incompatible with DKH2"):
        resolve_relativistic_basis("ZORA-def2-TZVP", "DKH2", [35, 1])
    with pytest.raises(ValueError, match="incompatible with ZORA"):
        resolve_relativistic_basis("DKH-def2-TZVP", "ZORA", [35, 1])
    with pytest.raises(ValueError, match="incompatible with ZORA"):
        resolve_relativistic_basis("cc-pVTZ-DK", "ZORA", [35, 1])


def test_relativistic_basis_swap_emits_warning_provenance_comment_and_metadata() -> None:
    """M3 regression: legitimate DKH2/ZORA swaps emit UserWarning, .inp provenance comment, and structured metadata."""
    hbr = [("Br", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.414)]
    auh = [("Au", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.524)]

    cases = [
        (hbr, "def2-TZVP", "DKH2", "DKH-def2-TZVP"),
        (hbr, "ma-def2-TZVP", "DKH2", "ma-DKH-def2-TZVP"),
        (hbr, "cc-pVTZ", "DKH2", "cc-pVTZ-DK"),
        (hbr, "aug-cc-pVTZ", "DKH2", "aug-cc-pVTZ-DK"),
        (hbr, "cc-pwCVTZ", "DKH2", "cc-pwCVTZ-DK"),
        (hbr, "cc-pVTZ", "ZORA", "ZORA-def2-TZVP"),
        (hbr, "def2-TZVP", "ZORA", "ZORA-def2-TZVP"),
        (auh, "def2-TZVP", "DKH2", "x2c-TZVPall"),
        (auh, "cc-pVTZ", "ZORA", "x2c-TZVPall"),
    ]
    for coords, req_basis, rel_mode, expected_eff in cases:
        with pytest.warns(UserWarning, match=expected_eff):
            deck, meta = generate_dlpno_ccsd_f12_with_metadata(
                coords, basis=req_basis, rel_mode=rel_mode
            )
        tokens = _header_tokens(deck)
        assert expected_eff in tokens
        assert meta["requested_basis"] == req_basis
        assert meta["effective_basis"] == expected_eff
        assert meta["basis_swapped"] is True
        assert len(meta["warnings"]) == 1
        assert meta["provenance_comment"] is not None
        assert meta["provenance_comment"].startswith("# CoChem-BENCH basis provenance:")
        assert req_basis in meta["provenance_comment"]
        assert expected_eff in meta["provenance_comment"]
        assert meta["provenance_comment"] in deck.splitlines()

    # Already-relativistic bases must NOT swap or emit warnings
    import warnings as _warnings

    for rel_basis, rel_mode in [
        ("DKH-def2-TZVP", "DKH2"),
        ("cc-pVTZ-DK", "DKH2"),
        ("SARC-DKH-TZVP", "DKH2"),
        ("x2c-TZVPall", "DKH2"),
        ("ZORA-def2-TZVP", "ZORA"),
        ("SARC-ZORA-TZVP", "ZORA"),
        ("x2c-TZVPall", "ZORA"),
    ]:
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            deck, meta = generate_dlpno_ccsd_f12_with_metadata(
                hbr, basis=rel_basis, rel_mode=rel_mode
            )
        assert caught == []
        assert meta["basis_swapped"] is False
        assert meta["effective_basis"] == rel_basis
        assert meta["warnings"] == []
        assert meta["provenance_comment"] is None
        assert "# CoChem-BENCH basis provenance:" not in deck

    # Counterpoise wrapper also exposes structured basis provenance metadata
    with pytest.warns(UserWarning, match="DKH-def2-SVP"):
        cp_res = generate_counterpoise_input(
            [("Br", 0.0, 0.0, 0.0)],
            [("H", 0.0, 0.0, 1.414)],
            basis="def2-SVP",
            rel_mode="auto",
        )
    assert cp_res["basis_provenance"]["effective_basis"] == "DKH-def2-SVP"
    assert cp_res["basis_provenance"]["basis_swapped"] is True


def test_build_orca_env_prepends_discovered_openmpi_bin_and_lib(tmp_path: Path) -> None:
    """DEF-S1-04 regression: _build_orca_env must prepend OpenMPI bin/ to PATH and lib/ to LD_LIBRARY_PATH."""
    from cochem_bench.orca_runner import _build_orca_env

    engines = tmp_path / "engines"
    orca_dir = engines / "orca-native-v1" / "orca_6_1_1_linux"
    mpi_root = engines / "openmpi4.1.8-native-v1"
    mpi_bin = mpi_root / "bin"
    mpi_lib = mpi_root / "lib"
    custom_mpi_bin = tmp_path / "custom_mpi" / "bin"
    custom_mpi_lib = tmp_path / "custom_mpi" / "lib"
    for folder in (orca_dir, mpi_bin, mpi_lib, custom_mpi_bin, custom_mpi_lib):
        folder.mkdir(parents=True)
    orca_bin = orca_dir / "orca"
    orca_bin.write_text("#!/bin/sh\n", encoding="utf-8")
    custom_mpirun = custom_mpi_bin / "mpirun"
    custom_mpirun.write_text("#!/bin/sh\n", encoding="utf-8")

    deck_pal2 = "! HF STO-3G\n%pal nprocs 2 end\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"

    # 1. Sibling openmpi*/bin and openmpi*/lib discovered automatically near ORCA binary
    saved_env = {
        k: os.environ.get(k)
        for k in ("COCHEM_ORCA_MPIRUN_BIN", "ORCA_MPI_BIN", "COCHEM_MPIRUN_BIN", "COCHEM_MPIEXEC_BIN", "COCHEM_ORCA_LD_LIBRARY_PATH")
    }
    try:
        for k in saved_env:
            os.environ.pop(k, None)
        env = _build_orca_env(orca_bin, deck_pal2)
        path_entries = env["PATH"].split(os.pathsep)
        assert str(orca_dir) == path_entries[0]
        assert str(mpi_bin) in path_entries
        assert env["OMP_NUM_THREADS"] == "1"
        if os.name == "posix":
            ld_entries = env["LD_LIBRARY_PATH"].split(":")
            assert str(mpi_lib) in ld_entries

        # 2. Explicit COCHEM_ORCA_MPIRUN_BIN takes priority ahead of sibling directories
        os.environ["COCHEM_ORCA_MPIRUN_BIN"] = str(custom_mpirun)
        env2 = _build_orca_env(orca_bin, deck_pal2)
        path_entries2 = env2["PATH"].split(os.pathsep)
        assert path_entries2[0] == str(orca_dir)
        assert path_entries2[1] == str(custom_mpi_bin)
        assert str(mpi_bin) in path_entries2
        if os.name == "posix":
            ld_entries2 = env2["LD_LIBRARY_PATH"].split(":")
            assert str(custom_mpi_lib) in ld_entries2
    finally:
        for k, v in saved_env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def test_live_orca6_parallel_nprocs2_execution(tmp_path: Path) -> None:
    """DEF-S1-04 live check: parallel nprocs=2 calculation succeeds with OpenMPI auto-discovery."""
    if os.name != "posix":
        pytest.skip("Parallel OpenMPI ORCA execution check targets Linux/WSL2 POSIX hosts")
    orca_bin = _live_orca_binary_or_skip()
    # HF has 4 correlated valence orbitals (10 pairs), satisfying ORCA_MDCI nprocs <= pairs.
    hf = [("F", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.917)]
    res = run_orca_calculation(
        hf,
        tmp_path / "hf_pal2",
        basis="cc-pVDZ-F12",
        is_opt=False,
        nprocs=2,
        job_name="hf_pal2",
        orca_bin=orca_bin,
    )
    assert res.returncode == 0
    assert res.orca_version.startswith("6.")
    assert -100.60 < res.final_single_point_energy_eh < -100.30
    assert res.audit["passed"] is True
