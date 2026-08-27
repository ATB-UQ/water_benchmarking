"""Command line entry point: water-bench <command>."""
from __future__ import annotations

import argparse
from pathlib import Path

from . import box, forcefield, gadi, gromacs, imd, protocol, report, trc
from .analysis import aggregate, density_hov, errors


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


def cmd_submit_gromacs(args) -> None:
    """Stage and submit the GROMACS ladder for one model on Setonix."""
    run_dir = _run_dir(args.model, "gromacs", args.root)
    gromacs.write_all(run_dir, args.model)
    box.build(run_dir)
    job_id = gromacs.submit(args.model, run_dir, dry_run=args.dry_run)
    if job_id:
        print(f"{args.model}: Setonix job {job_id}")


def cmd_analyse(args) -> None:
    """Compute every property for one model and engine."""
    run_dir = _run_dir(args.model, args.engine, args.root)
    if getattr(args, "collect", False) and args.engine == "gromacs":
        gromacs.collect(args.model, run_dir)
    results = analyse_run(args.model, args.engine, run_dir)
    for key, value in sorted(results.values.items()):
        print(f"{key:14s} {value:.6g}")


def analyse_run(model: str, engine: str, run_dir: Path) -> report.Results:
    """Compute every property for one finished run, whichever engine produced it.

    The two engines differ only in how the numbers are read off disk: GROMOS keeps
    energies in .tre and positions in .trc, GROMACS in .edr and .xtc.  Once the
    trajectory is in the shared g96 form, diffusion, rotation and the dielectric
    constant run through exactly the same code for both -- which is what makes an
    engine-to-engine difference mean something.
    """
    results = report.Results(model=model, engine=engine)
    charges = forcefield.EXPECTED_CHARGES[model]

    if engine == "gromos":
        energies = sorted(run_dir.glob("md_*.tre*"))
        if energies:
            thermo = density_hov.analyse(energies, model, run_dir / "analysis")
            _record_thermodynamics(results, thermo)
        trajectories = sorted(run_dir.glob("md_*.trc*"))
    else:
        energies = sorted(run_dir.glob("md_*.edr"))
        if energies:
            series = gromacs.energy_series(energies, ("Density", "Potential"))
            _record_gromacs_thermodynamics(results, series, model)
        trajectories = _converted_trajectories(run_dir)

    if trajectories:
        summary = aggregate.analyse_run(trajectories, charges)
        results.values["diffusion"] = summary.diffusion.d_corrected
        results.values["diffusion_pbc"] = summary.diffusion.d_pbc
        results.uncertainties["diffusion"] = summary.diffusion.d_error
        results.values["tau2_HH"] = summary.tau2["HH"]
        results.values["tau2_OH"] = summary.tau2["OH"]
        results.values["tau1_dipole"] = summary.tau1["dipole"]
        results.values["dielectric"] = summary.dielectric.epsilon
        results.uncertainties["dielectric"] = summary.dielectric.epsilon_error

    return results


def _record_thermodynamics(results: report.Results, thermo) -> None:
    results.values["density"] = thermo.density.mean
    results.uncertainties["density"] = thermo.density.error
    results.values["hov"] = thermo.hov
    results.uncertainties["hov"] = thermo.hov_error
    if thermo.hov_polarisation_corrected is not None:
        results.values["hov_polarisation_corrected"] = thermo.hov_polarisation_corrected


def _record_gromacs_thermodynamics(results: report.Results, series, model: str) -> None:
    """The same two numbers, from gmx energy instead of ene_ana."""
    density = errors.block_average(errors.drop_equilibration(series["Density"]))
    per_molecule = errors.drop_equilibration(series["Potential"]) / protocol.N_WATERS
    energy = errors.block_average(per_molecule)

    results.values["density"] = density.mean
    results.uncertainties["density"] = density.error
    hov = -energy.mean + density_hov.GAS_CONSTANT * protocol.TEMPERATURE
    results.values["hov"] = hov
    results.uncertainties["hov"] = energy.error
    if model == "spce":
        results.values["hov_polarisation_corrected"] = (
            hov - density_hov.SPCE_POLARISATION_CORRECTION
        )


def _converted_trajectories(run_dir: Path) -> list[Path]:
    """Convert each .xtc to g96 on demand, reusing anything already converted.

    The text form is an order of magnitude larger than the .xtc it came from, so
    these are treated as scratch: `water-bench analyse --clean` removes them again.
    """
    converted = []
    for xtc in sorted(run_dir.glob("md_*.xtc")):
        tpr = xtc.with_suffix(".tpr")
        g96 = xtc.with_suffix(".g96")
        if not tpr.exists():
            raise FileNotFoundError(f"{tpr} is needed to convert {xtc.name}")
        if not g96.exists():
            gromacs.to_g96(xtc, tpr, g96)
        converted.append(g96)
    return converted


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

    submit = sub.add_parser("submit-gromacs",
                            help="stage and sbatch the GROMACS ladder on Setonix")
    submit.add_argument("--model", required=True, choices=sorted(protocol.MODELS))
    submit.add_argument("--dry-run", action="store_true")
    submit.set_defaults(func=cmd_submit_gromacs)

    analyse = sub.add_parser("analyse", help="compute properties for one run")
    analyse.add_argument("--model", required=True, choices=sorted(protocol.MODELS))
    analyse.add_argument("--engine", default="gromos", choices=("gromos", "gromacs"))
    analyse.add_argument("--collect", action="store_true",
                         help="pull results from Setonix first (gromacs only)")
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
