"""Command line entry point: water-bench <command>."""
from __future__ import annotations

import argparse
from pathlib import Path

from . import box, forcefield, gadi, gromacs, imd, protocol, report, trc
from .analysis import density_hov, dielectric, diffusion, rotation


def _run_dir(model: str, engine: str, root: Path) -> Path:
    return root / model / engine


def cmd_build(args) -> None:
    """Write topology, box and inputs for both engines without running anything."""
    for model in args.models:
        gromos_dir = _run_dir(model, "gromos", args.root)
        gadi.prepare(model, gromos_dir)
        print(f"{model}: GROMOS inputs in {gromos_dir}")

        gmx_dir = _run_dir(model, "gromacs", args.root)
        gromacs.write_all(gmx_dir, model)
        water_box = box.build(gmx_dir)
        (gmx_dir / "run.slurm").write_text(
            gromacs.slurm_script(model, str(gmx_dir))
        )
        print(f"{model}: GROMACS inputs in {gmx_dir} (box edge {water_box.edge:.4f} nm)")


def cmd_run(args) -> None:
    """Run the GROMOS ladder for one model on Gadi."""
    run_dir = _run_dir(args.model, "gromos", args.root)
    gadi.run_model(args.model, run_dir, dry_run=args.dry_run)


def cmd_analyse(args) -> None:
    """Compute every property for one model and engine."""
    run_dir = _run_dir(args.model, args.engine, args.root)
    results = analyse_run(args.model, args.engine, run_dir)
    for key, value in sorted(results.values.items()):
        print(f"{key:14s} {value:.6g}")


def analyse_run(model: str, engine: str, run_dir: Path) -> report.Results:
    """Run all four analyses over one finished run directory."""
    results = report.Results(model=model, engine=engine)

    energies = sorted(run_dir.glob("md_*.tre*"))
    if energies:
        thermo = density_hov.analyse(energies, model, run_dir / "analysis")
        results.values["density"] = thermo.density.mean
        results.uncertainties["density"] = thermo.density.error
        results.values["hov"] = thermo.hov
        results.uncertainties["hov"] = thermo.hov_error
        if thermo.hov_polarisation_corrected is not None:
            results.values["hov_polarisation_corrected"] = thermo.hov_polarisation_corrected

    trajectories = sorted(run_dir.glob("md_*.trc*"))
    if trajectories:
        charges = forcefield.solvent_charges(run_dir / f"{model}.top")

        lags, msd, edge = diffusion.mean_squared_displacement(trc.read_all(trajectories))
        d = diffusion.diffusion_from_msd(lags, msd, edge)
        results.values["diffusion"] = d.d_corrected
        results.values["diffusion_pbc"] = d.d_pbc
        results.uncertainties["diffusion"] = d.d_error

        for vector, key in (("HH", "tau2_HH"), ("OH", "tau2_OH")):
            rot = rotation.analyse(trc.read_all(trajectories), vector=vector)
            results.values[key] = rot.tau2
        dipole_rot = rotation.analyse(trc.read_all(trajectories), vector="dipole")
        results.values["tau1_dipole"] = dipole_rot.tau1

        eps = dielectric.dielectric_constant(trc.read_all(trajectories), charges)
        results.values["dielectric"] = eps.epsilon
        results.uncertainties["dielectric"] = eps.epsilon_error

    return results


def cmd_report(args) -> None:
    """Build the comparison table across every model and engine that has results."""
    collected = []
    for model in args.models:
        for engine in args.engines:
            run_dir = _run_dir(model, engine, args.root)
            if not run_dir.exists():
                continue
            collected.append(analyse_run(model, engine, run_dir))
    if not collected:
        raise SystemExit(f"no finished runs found under {args.root}")
    path = report.write_report(collected, args.output)
    print(f"wrote {path}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="water-bench", description=__doc__
    )
    parser.add_argument("--root", type=Path, default=protocol.RUN_ROOT,
                        help="directory holding <model>/<engine> run directories")
    sub = parser.add_subparsers(dest="command", required=True)

    build = sub.add_parser("build", help="write inputs for both engines")
    build.add_argument("--models", nargs="+", default=sorted(protocol.MODELS))
    build.set_defaults(func=cmd_build)

    run = sub.add_parser("run", help="run the GROMOS ladder on Gadi")
    run.add_argument("--model", required=True, choices=sorted(protocol.MODELS))
    run.add_argument("--dry-run", action="store_true",
                     help="print the shim invocations instead of submitting")
    run.set_defaults(func=cmd_run)

    analyse = sub.add_parser("analyse", help="compute properties for one run")
    analyse.add_argument("--model", required=True, choices=sorted(protocol.MODELS))
    analyse.add_argument("--engine", default="gromos", choices=("gromos", "gromacs"))
    analyse.set_defaults(func=cmd_analyse)

    rep = sub.add_parser("report", help="build the comparison table")
    rep.add_argument("--models", nargs="+", default=sorted(protocol.MODELS))
    rep.add_argument("--engines", nargs="+", default=("gromos", "gromacs"))
    rep.add_argument("--output", type=Path, default=Path("results"))
    rep.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
