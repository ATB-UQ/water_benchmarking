"""Analyse the dielectric-investigation runs with the main-run code, unchanged.

Each diagnostic changed one thing about the protocol -- thermostat, cutoff, box
size, electrostatics -- to find why the reaction-field dielectric constant read
twice the published value.  They were 1 ns each and only ever needed for eps,
but there is no reason to report only that: put through exactly the analysis
the main runs get, every one yields density, heat of vaporisation, diffusion and
rotational relaxation too, and the table then says how much each of those
depends on the cutoff, the box and the boundary condition.  That is a result in
its own right for a protocol that will be run on solvated systems of every size.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import experiment, gromacs, protocol, report, trc
from .analysis import aggregate, density_hov, errors

MANIFEST = Path(__file__).resolve().parents[2] / "diagnostics" / "manifest.json"


@dataclass
class Diagnostic:
    tag: str
    label: str
    n_molecules: int
    cutoff: float
    electrostatics: str
    changed: str
    values: dict = field(default_factory=dict)
    uncertainties: dict = field(default_factory=dict)
    summary: object = None


def load_manifest(path: Path = MANIFEST) -> list[Diagnostic]:
    return [Diagnostic(**entry) for entry in json.loads(path.read_text())]


def analyse(diagnostic: Diagnostic, run_dir: Path, model: str = "spc") -> Diagnostic:
    """All five properties for one diagnostic run."""
    n = diagnostic.n_molecules
    edr = run_dir / f"{diagnostic.tag}.edr"
    xtc = run_dir / f"{diagnostic.tag}.xtc"
    tpr = run_dir / f"{diagnostic.tag}.tpr"

    series = gromacs.energy_series([edr], ("Density", "Potential"))
    density = errors.block_average(errors.drop_equilibration(series["Density"]))
    per_molecule = errors.drop_equilibration(series["Potential"]) / n
    energy = errors.block_average(per_molecule)
    diagnostic.values["density"] = density.mean
    diagnostic.uncertainties["density"] = density.error
    diagnostic.values["hov"] = -energy.mean + density_hov.GAS_CONSTANT * protocol.TEMPERATURE
    diagnostic.uncertainties["hov"] = energy.error

    def frames():
        g96 = xtc.with_suffix(".g96")
        gromacs.to_g96(xtc, tpr, g96)
        try:
            yield from trc.read_frames(g96, n_molecules=n)
        finally:
            g96.unlink(missing_ok=True)

    # The 8x box would need eight times the memory for the rotational
    # correlation functions; 2048 molecules of it is the same statistics as the
    # main runs and is what gets sampled.
    summary = aggregate.analyse_run([frames], _charges(model), n_molecules=n)
    diagnostic.summary = summary
    diagnostic.values["diffusion"] = summary.diffusion.d_corrected
    diagnostic.values["diffusion_pbc"] = summary.diffusion.d_pbc
    diagnostic.uncertainties["diffusion"] = summary.diffusion.d_error
    diagnostic.values["tau2_HH"] = summary.tau2["HH"]
    diagnostic.values["tau2_OH"] = summary.tau2["OH"]
    diagnostic.values["tau1_dipole"] = summary.tau1["dipole"]
    diagnostic.values["dielectric"] = summary.dielectric.epsilon
    diagnostic.uncertainties["dielectric"] = summary.dielectric.epsilon_error
    diagnostic.values["dielectric_y"] = summary.dielectric.y
    diagnostic.values["dielectric_neumann"] = summary.dielectric.epsilon_neumann
    return diagnostic


def _charges(model: str):
    from .forcefield import EXPECTED_CHARGES
    return EXPECTED_CHARGES[model]


#: Columns of the diagnostics table.  y and the Neumann value belong here even
#: though they were dropped from the summary: this table is the demonstration
#: that the fluctuation is the same everywhere and only the relation moves.
COLUMNS = (
    ("density", "rho", "{:.1f}"),
    ("hov", "dH_vap", "{:.2f}"),
    ("diffusion", "D", "{:.2f}"),
    ("tau2_HH", "tau2(HH)", "{:.2f}"),
    ("dielectric_y", "y", "{:.1f}"),
    ("dielectric", "eps=1+y", "{:.1f}"),
    ("dielectric_neumann", "eps Neumann(61)", "{:.0f}"),
)


def table(protocol_row: report.Results | None, diagnostics: list[Diagnostic]) -> list[list[str]]:
    header = ["run", "changed from protocol", "R_c (nm)", "electrostatics", "N"] + \
             [label for _key, label, _fmt in COLUMNS]
    rows = [header]
    if protocol_row is not None:
        rows.append(["protocol (10 ns)", "-", f"{protocol.CUTOFF}", f"RF eps_rf = {protocol.EPSILON_RF:.0f}",
                     str(protocol.N_WATERS)] + [_cell(protocol_row.values, key, fmt) for key, _l, fmt in COLUMNS])
    for d in diagnostics:
        cells = [_cell(d.values, key, fmt) for key, _l, fmt in COLUMNS]
        if "61" not in d.electrostatics:
            # Neumann's finite-eps_rf relation has no meaning for a conducting
            # boundary or PME; printing it there would only invite comparison.
            cells[[k for k, _l, _f in COLUMNS].index("dielectric_neumann")] = "-"
        rows.append([d.label + " (1 ns)", d.changed, f"{d.cutoff}", d.electrostatics,
                     str(d.n_molecules)] + cells)
    rows.append(["experiment", "-", "-", "-", "-"] + [
        fmt.format(experiment.EXPERIMENT[key].value * report.SCALE.get(key, 1.0))
        if key in experiment.EXPERIMENT else "-" for key, _l, fmt in COLUMNS])
    return rows


def _cell(values: dict, key: str, fmt: str) -> str:
    value = values.get(key)
    if value is None:
        return "-"
    return fmt.format(value * report.SCALE.get(key, 1.0))


def write(protocol_row, diagnostics: list[Diagnostic], output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = table(protocol_row, diagnostics)
    report.write_csv(rows, output_dir / "diagnostics.csv")
    text = [
        "# Protocol diagnostics (SPC, 298.15 K, 1 atm)",
        "",
        "Each run changes one thing relative to the benchmark protocol; all are analysed with",
        "the same code as the main runs. Units: rho kg m^-3, dH_vap kJ mol^-1, D 1e-9 m^2 s^-1",
        "(Yeh-Hummer corrected), tau2 ps. y is the dimensionless box-dipole fluctuation;",
        "eps = 1 + y is the conducting-boundary relation, eps Neumann(61) the finite-eps_rf one",
        "(see README: the latter sits near its pole here and is not the number to quote).",
        "",
        report.to_markdown(rows),
        "",
    ]
    path = output_dir / "diagnostics.md"
    path.write_text("\n".join(text))
    return path
