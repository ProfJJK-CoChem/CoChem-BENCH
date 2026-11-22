"""Command-line interface for CoChem-BENCH (`cochem-bench`)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from .extrapolate import cbs_extrapolate
from .orca_parser import parse_orca_energy
from .provider import bench_provider
from .thermochemistry import compute_relative_energies
from .verifier import audit_diagnostics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cochem-bench",
        description="CoChem-BENCH: CBS Extrapolation, Composite Thermochemistry, and ORCA Parsing.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("info", help="Print bench_provider metadata as JSON.")

    p_cbs = sub.add_parser("cbs", help="Run two-point Helgaker/Neese-Valeev CBS extrapolation.")
    p_cbs.add_argument("--e-scf-x", type=float, required=True, help="SCF energy at cardinal X (Eh).")
    p_cbs.add_argument("--e-scf-y", type=float, required=True, help="SCF energy at cardinal Y (Eh).")
    p_cbs.add_argument("--e-corr-x", type=float, default=None, help="Correlation energy at X (Eh).")
    p_cbs.add_argument("--e-corr-y", type=float, default=None, help="Correlation energy at Y (Eh).")
    p_cbs.add_argument("--e-tot-x", type=float, default=None, help="Total energy at X (Eh).")
    p_cbs.add_argument("--e-tot-y", type=float, default=None, help="Total energy at Y (Eh).")
    p_cbs.add_argument("-x", "--cardinal-x", type=int, default=3)
    p_cbs.add_argument("-y", "--cardinal-y", type=int, default=4)
    p_cbs.add_argument("--basis-family", default="cc-pVnZ")
    p_cbs.add_argument("--scheme", choices=["helgaker", "neese_valeev", "f12"], default="helgaker")

    p_orca = sub.add_parser("parse-orca", help="Parse an ORCA output file.")
    p_orca.add_argument("path", type=Path, help="Path to ORCA .out file.")

    p_audit = sub.add_parser("audit", help="Run diagnostic audit on an ORCA output file.")
    p_audit.add_argument("path", type=Path, help="Path to ORCA .out file.")
    p_audit.add_argument("--multiplicity", type=int, default=1)
    p_audit.add_argument("--claimed-error-pct", type=float, default=None)
    p_audit.add_argument("--has-core-valence", action="store_true")

    p_rel = sub.add_parser("rel-energy", help="Compute relative energy ranking table from JSON.")
    p_rel.add_argument("input_json", type=Path, help="Path to JSON file containing species list.")
    p_rel.add_argument("--reference", default=None)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "info":
        print(json.dumps(bench_provider(), indent=2))
        return 0

    if args.command == "cbs":
        res = cbs_extrapolate(
            e_scf_x=args.e_scf_x,
            e_scf_y=args.e_scf_y,
            e_corr_x=args.e_corr_x,
            e_corr_y=args.e_corr_y,
            e_tot_x=args.e_tot_x,
            e_tot_y=args.e_tot_y,
            x=args.cardinal_x,
            y=args.cardinal_y,
            basis_family=args.basis_family,
            scheme=args.scheme,
        )
        print(res.model_dump_json(indent=2))
        return 0

    if args.command == "parse-orca":
        res = parse_orca_energy(args.path)
        print(res.model_dump_json(indent=2))
        return 0

    if args.command == "audit":
        res = audit_diagnostics(
            args.path,
            spin_multiplicity=args.multiplicity,
            claimed_error_pct=args.claimed_error_pct,
            has_core_valence=args.has_core_valence,
        )
        print(json.dumps(res, indent=2))
        return 0

    if args.command == "rel-energy":
        raw = json.loads(args.input_json.read_text(encoding="utf-8"))
        species = raw["species"] if isinstance(raw, dict) and "species" in raw else raw
        table = compute_relative_energies(species, reference=args.reference)
        print(table.model_dump_json(indent=2))
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
