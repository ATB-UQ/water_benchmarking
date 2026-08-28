"""Command line entry point: water-bench <command>."""
from __future__ import annotations

import argparse
from pathlib import Path

from . import box, gadi, gromacs, imd, protocol, report, trc
from .analysis import aggregate, density_hov, errors


def _run_dir(model: str, engine: str, root: Path) -> Path:
    return root / model / engine


def cmd_build(args) -> None:
    """Write topology, box and inputs for every engine a model is run with.

    Not every model is run with both.  OPC3 has no 54A7 building block, so there is
    no GROMOS topology to build for it; the registry says which engines apply and
    this loop asks rather than assuming.
    """
    for model in args.models:
        engines = protocol.engines_for(model)

        if "gromos" in engines:
            gromos_dir = _run_dir(model, "gromos", args.root)
            gadi.prepare(model, gromos_dir)
            print(f"{model}: GROMOS inputs in {gromos_dir}")

        if "gromacs" in engines:
            gmx_dir = _run_dir(model, "gromacs", args.root)
            gromacs.write_all(gmx_dir, model)
            water_box = box.build(gmx_dir)
            (gmx_dir / "run.slurm").write_text(
                gromacs.slurm_script(model, str(gmx_dir))
            )
            print(f"{model}: GROMACS inputs in {gmx_dir} (box edge {water_box.edge:.4f} nm)")

        skipped = [e for e in ("gromos", "gromacs") if e not in engines]
        if skipped:
            print(f"{model}: not run with {', '.join(skipped)}")


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
    # Not every value is a number: dielectric_relation records which relation was
    # used to turn the dipole fluctuation into eps, and formatting a str with
    # ":.6g" raises -- after the whole trajectory analysis has been done and with
    # nothing written to disk, so the work is simply lost.
    width = max((len(key) for key in results.values), default=14)
    for key, value in sorted(results.values.items()):
        rendered = value if isinstance(value, str) else f"{value:.6g}"
        print(f"{key:{width}s}  {rendered}")


def analyse_run(model: str, engine: str, run_dir: Path) -> report.Results:
    """Compute every property for one finished run, whichever engine produced it.

    The two engines differ only in how the numbers are read off disk: GROMOS keeps
    energies in .tre and positions in .trc, GROMACS in .edr and .xtc.  Once the
    trajectory is in the shared g96 form, diffusion, rotation and the dielectric
    constant run through exactly the same code for both -- which is what makes an
    engine-to-engine difference mean something.
    """
    results = report.Results(model=model, engine=engine)
    entry = protocol.model(model)
    charges = entry.charges
    # GROMOS production is ten independent replicates with fresh velocities, so
    # the head of each is discarded; the GROMACS segments are one continuous
    # trajectory already past its equilibration and lose nothing.
    discard = 0.1 if engine == "gromos" else 0.0

    if engine == "gromos":
        energies = sorted(run_dir.glob("md_*.tre*"))
        if energies:
            thermo = density_hov.analyse(energies, model, run_dir / "analysis", discard=discard)
            _record_thermodynamics(results, thermo)
            results.density_series = thermo.density_series
        trajectories = sorted(run_dir.glob("md_*.trc*"))
    else:
        energies = sorted(run_dir.glob("md_*.edr"))
        if energies:
            series = gromacs.energy_series(energies, ("Density", "Potential"))
            _record_gromacs_thermodynamics(results, series, model)
            results.density_series = series["Density"]
        trajectories = _gromacs_segments(run_dir)

    if trajectories:
        summary = aggregate.analyse_run(
            trajectories, charges, discard=discard, r_oh=entry.r_oh
        )
        results.summary = summary
        results.values["diffusion"] = summary.diffusion.d_corrected
        results.values["diffusion_pbc"] = summary.diffusion.d_pbc
        results.uncertainties["diffusion"] = summary.diffusion.d_error
        results.values["tau2_HH"] = summary.tau2["HH"]
        results.values["tau2_OH"] = summary.tau2["OH"]
        results.values["tau1_dipole"] = summary.tau1["dipole"]
        # The relation that turns the dipole fluctuation into eps depends on the
        # boundary the engine actually realises, and the two engines differ.
        # GROMOS honours the finite eps_rf: its fluctuation is 35% below the
        # conducting-boundary value, exactly as Neumann's relation predicts for
        # eps ~ 70, and that relation is the right one for it.  GROMACS's GPU
        # runs show no eps_rf dependence at all (eps_rf = 61 and infinity give
        # the same fluctuation), so for them the conducting-boundary relation
        # eps = 1 + y is the one that describes the box.  Read the other way
        # round, GROMOS gives 44 and GROMACS 140; read this way both give the
        # published values.  The choice is recorded beside the number.
        eps = summary.dielectric
        if engine == "gromos":
            results.values["dielectric"] = eps.epsilon_neumann
            results.values["dielectric_relation"] = f"Neumann, eps_rf = {protocol.EPSILON_RF:.0f}"
            # Propagate the y uncertainty through the Neumann relation's slope.
            results.uncertainties["dielectric"] = eps.epsilon_error * eps.neumann_sensitivity
        else:
            results.values["dielectric"] = eps.epsilon
            results.values["dielectric_relation"] = "conducting boundary, eps = 1 + y"
            results.uncertainties["dielectric"] = eps.epsilon_error
        results.values["dielectric_y"] = eps.y
        results.values["dielectric_conducting"] = eps.epsilon
        results.values["dielectric_neumann"] = eps.epsilon_neumann
        results.values["dielectric_neumann_sensitivity"] = eps.neumann_sensitivity

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
    polarisation = density_hov.POLARISATION_CORRECTION.get(model)
    if polarisation is not None:
        results.values["hov_polarisation_corrected"] = hov - polarisation


def _gromacs_segments(run_dir: Path) -> list:
    """One converting, self-cleaning frame source per .xtc segment.

    The g96 text form is roughly thirteen times the size of the .xtc it comes
    from, so each is made, streamed and deleted in turn: peak disk stays at one
    segment rather than the 28 GB a whole model would otherwise need.
    """
    sources = []
    for xtc in sorted(run_dir.glob("md_*.xtc")):
        tpr = xtc.with_suffix(".tpr")
        if not tpr.exists():
            raise FileNotFoundError(f"{tpr} is needed to convert {xtc.name}")

        def source(xtc=xtc, tpr=tpr):
            g96 = xtc.with_suffix(".g96")
            gromacs.to_g96(xtc, tpr, g96)
            try:
                yield from trc.read_frames(g96)
            finally:
                g96.unlink(missing_ok=True)

        sources.append(source)
    return sources


def cmd_diagnostics(args) -> None:
    """Analyse the dielectric-investigation runs and write their table."""
    import pickle

    from . import diagnostics

    run_dir = args.root / "diagnostics"
    rows = []
    for diagnostic in diagnostics.load_manifest():
        rows.append(diagnostics.analyse(diagnostic, run_dir))
        print(f"{diagnostic.tag}: eps = {diagnostic.values['dielectric']:.1f}", flush=True)
    protocol_row = None
    if args.protocol_results and Path(args.protocol_results).exists():
        protocol_row = next(r for r in pickle.load(open(args.protocol_results, "rb"))
                            if r.model == "spc")
    path = diagnostics.write(protocol_row, rows, args.output)
    print(f"wrote {path}")


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
    run.add_argument("--model", required=True,
                     choices=protocol.models_for_engine("gromos"))
    run.add_argument("--dry-run", action="store_true",
                     help="print the shim invocations instead of submitting")
    run.set_defaults(func=cmd_run)

    submit = sub.add_parser("submit-gromacs",
                            help="stage and sbatch the GROMACS ladder on Setonix")
    submit.add_argument("--model", required=True,
                        choices=protocol.models_for_engine("gromacs"))
    submit.add_argument("--dry-run", action="store_true")
    submit.set_defaults(func=cmd_submit_gromacs)

    analyse = sub.add_parser("analyse", help="compute properties for one run")
    analyse.add_argument("--model", required=True, choices=sorted(protocol.MODELS))
    analyse.add_argument("--engine", default="gromos", choices=("gromos", "gromacs"))
    analyse.add_argument("--collect", action="store_true",
                         help="pull results from Setonix first (gromacs only)")
    analyse.set_defaults(func=cmd_analyse)

    diag = sub.add_parser("diagnostics", help="analyse the dielectric-investigation runs")
    diag.add_argument("--output", type=Path, default=Path("results"))
    diag.add_argument("--protocol-results", type=Path, default=None,
                      help="pickled main-run Results, for the reference row")
    diag.set_defaults(func=cmd_diagnostics)

    rep = sub.add_parser("report", help="build the comparison table")
    rep.add_argument("--models", nargs="+", default=sorted(protocol.MODELS))
    rep.add_argument("--engines", nargs="+", default=("gromos", "gromacs"))
    rep.add_argument("--output", type=Path, default=Path("results"))
    rep.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
