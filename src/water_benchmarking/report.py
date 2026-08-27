"""Collect the per-model, per-engine results into one comparison table."""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from . import experiment, protocol

#: Display order and formatting for the benchmark's five properties.
PROPERTIES = (
    ("density", "Density", "kg m^-3", "{:.1f}"),
    ("hov", "dH_vap", "kJ mol^-1", "{:.2f}"),
    ("diffusion", "Self-diffusion D", "1e-9 m^2 s^-1", "{:.2f}"),
    ("tau2_HH", "Rot. corr. time tau2(HH)", "ps", "{:.2f}"),
    ("tau1_dipole", "Rot. corr. time tau1(mu)", "ps", "{:.2f}"),
    # eps = 1 + y.  The fluctuation y and the Neumann(eps_rf) value are still
    # computed and kept in Results.values for anyone who wants them; the README
    # explains why they are not the number to quote.
    ("dielectric", "Dielectric constant", "-", "{:.1f}"),
)

#: Diffusion is reported in units of 1e-9 m^2/s to keep the table readable.
SCALE = {"diffusion": 1e9}


@dataclass
class Results:
    """Everything measured for one (model, engine) pair."""

    model: str
    engine: str
    values: dict = field(default_factory=dict)
    uncertainties: dict = field(default_factory=dict)
    #: The trajectory analysis this came from, kept so the report can plot the
    #: curves behind the numbers -- a converged epsilon and an unconverged one
    #: look identical in a table.
    summary: object = None
    density_series: object = None

    @property
    def key(self) -> str:
        return f"{self.model.upper()}/{self.engine}"


def _format(key: str, value: float | None, error: float | None, template: str) -> str:
    if value is None:
        return "-"
    scaled = value * SCALE.get(key, 1.0)
    text = template.format(scaled)
    if error is not None and error == error:      # not NaN
        # Two significant figures for the error whatever the value's format:
        # a 0.003 kJ/mol uncertainty printed at two decimals reads as "+/- 0.00".
        text += f" +/- {error * SCALE.get(key, 1.0):.2g}"
    return text


def build_table(results: list[Results]) -> list[list[str]]:
    header = ["Property", "Unit"] + [r.key for r in results] + ["Experiment", "Source"]
    rows = [header]
    for key, label, unit, template in PROPERTIES:
        reference = experiment.EXPERIMENT.get(key)
        row = [label, unit]
        for result in results:
            row.append(
                _format(key, result.values.get(key), result.uncertainties.get(key), template)
            )
        if reference:
            row.append(template.format(reference.value * SCALE.get(key, 1.0)))
            row.append(reference.source)
        else:
            row += ["-", "-"]
        rows.append(row)
    return rows


def deviation_table(results: list[Results]) -> list[list[str]]:
    """Percentage deviation from experiment, and whether the model looks reproduced."""
    header = ["Property"] + [r.key for r in results]
    rows = [header]
    for key, label, _unit, _template in PROPERTIES:
        if key not in experiment.EXPERIMENT or key in experiment.NOT_DIRECTLY_COMPARABLE:
            continue
        row = [label]
        for result in results:
            value = result.values.get(key)
            if value is None:
                row.append("-")
                continue
            percent = experiment.deviation(value, key)
            known = experiment.within_literature(value, result.model, key)
            flag = "" if known is None else ("" if known else "  [!]")
            row.append(f"{percent:+.1f}%{flag}")
        rows.append(row)
    return rows


def _notes(results: list[Results]) -> list[str]:
    """Footnotes for flagged values whose cause is understood."""
    lines = []
    for result in results:
        for key, _label, _unit, _template in PROPERTIES:
            note = experiment.NOTES.get((result.model, key))
            value = result.values.get(key)
            if note and value is not None and experiment.within_literature(value, result.model, key) is False:
                lines.append(f"- {result.key}, {key}: {note}")
    return lines + [""] if lines else []


def write_csv(rows: list[list[str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        csv.writer(handle).writerows(rows)


def to_markdown(rows: list[list[str]]) -> str:
    widths = [max(len(str(row[i])) for row in rows) for i in range(len(rows[0]))]
    lines = []
    for index, row in enumerate(rows):
        lines.append("| " + " | ".join(str(c).ljust(widths[i]) for i, c in enumerate(row)) + " |")
        if index == 0:
            lines.append("|" + "|".join("-" * (w + 2) for w in widths) + "|")
    return "\n".join(lines)


def write_plots(results: list[Results], output_dir: Path) -> list[Path]:
    """Draw every figure for which the runs carry data."""
    from .analysis import plots

    runs = {r.key: r.summary for r in results if r.summary is not None}
    written = []
    if runs:
        written.append(plots.dielectric_convergence(runs, output_dir / "dielectric.png"))
        written.append(plots.mean_squared_displacement(runs, output_dir / "msd.png"))
        written.append(plots.rotational_correlation(runs, output_dir / "c2_HH.png"))
    series = {r.key: r.density_series for r in results if r.density_series is not None}
    if series:
        written.append(plots.density_series(series, output_dir / "density.png"))
    return written


def write_report(results: list[Results], output_dir: Path) -> Path:
    """Write summary.csv, summary.md and the figures; return the markdown path."""
    output_dir.mkdir(parents=True, exist_ok=True)
    table = build_table(results)
    write_csv(table, output_dir / "summary.csv")
    figures = write_plots(results, output_dir)

    text = [
        "# Water model benchmark",
        "",
        f"{protocol.N_WATERS} water molecules, {protocol.TEMPERATURE} K, 1 atm, "
        f"{protocol.PRODUCTION_SEGMENTS} x "
        f"{protocol.PRODUCTION_STEPS * protocol.TIMESTEP / 1000:.0f} ns production.",
        f"Single-range cutoff {protocol.CUTOFF} nm, reaction field "
        f"eps_rf = {protocol.EPSILON_RF:.0f}, timestep {protocol.TIMESTEP * 1000:.0f} fs.",
        "",
        "## Results",
        "",
        to_markdown(table),
        "",
        "## Deviation from experiment",
        "",
        "`[!]` marks a value outside the published range for that model, which points",
        "at the setup rather than at the model -- unless a note below says otherwise.",
        "tau1(mu) is omitted: the experimental Debye time is a collective quantity and",
        "not directly comparable to the single-molecule correlation time simulated.",
        "",
        to_markdown(deviation_table(results)),
        "",
        *_notes(results),
    ]
    if figures:
        text += ["## Figures", ""]
        text += [f"![{f.stem}]({f.name})" for f in figures]
        text += [""]
    path = output_dir / "summary.md"
    path.write_text("\n".join(text))
    return path
