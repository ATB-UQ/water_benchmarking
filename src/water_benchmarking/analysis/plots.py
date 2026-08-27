"""Figures for the benchmark.

Four plots, each chosen because it shows something a single number in the summary
table cannot: whether the quantity had converged, and whether the two engines
agree on the whole curve rather than only on its summary.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

import matplotlib
import numpy as np

matplotlib.use("Agg")           # no display on the server
import matplotlib.pyplot as plt  # noqa: E402

from .. import experiment  # noqa: E402

#: One colour per model, consistent across every figure; engines are distinguished
#: by line style so a model's two curves read as a pair.
MODEL_COLOURS = {"spc": "#0072B2", "spce": "#D55E00"}
ENGINE_STYLES = {"GROMOS": "-", "GROMACS": "--"}


def _style(key: str) -> dict:
    model, _, engine = key.partition("/")
    return {
        "color": MODEL_COLOURS.get(model.lower(), "#666666"),
        "linestyle": ENGINE_STYLES.get(engine, "-"),
        "label": key,
        "linewidth": 1.6,
    }


def _finish(axes, path: Path) -> Path:
    axes.legend(frameon=False, fontsize=9)
    axes.spines[["top", "right"]].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    axes.figure.tight_layout()
    axes.figure.savefig(path, dpi=150)
    plt.close(axes.figure)
    return path


def dielectric_convergence(runs: Mapping[str, object], output: Path) -> Path:
    """Running epsilon against time -- the plot that says whether to believe it.

    A flat tail means the estimate has converged; a curve still climbing at the end
    of the run means the number in the table is a lower bound, not a result.
    """
    _figure, axes = plt.subplots(figsize=(6, 4))
    for key, summary in runs.items():
        eps = summary.dielectric
        if len(eps.running):
            axes.plot(eps.times, eps.running, **_style(key))
    axes.axhline(experiment.EXPERIMENT["dielectric"].value, color="k",
                 linestyle=":", linewidth=1, label="experiment")
    axes.set_xlabel("simulated time (ps)")
    axes.set_ylabel(r"$\epsilon$ estimated from 0 to $t$")
    axes.set_title("Dielectric constant convergence")
    return _finish(axes, output)


def mean_squared_displacement(runs: Mapping[str, object], output: Path) -> Path:
    """MSD with the Einstein fitting window marked.

    The fit is only meaningful where the curve is straight; showing the window
    makes it obvious whether it was placed on a linear stretch.
    """
    _figure, axes = plt.subplots(figsize=(6, 4))
    low = high = None
    for key, summary in runs.items():
        d = summary.diffusion
        axes.plot(d.lag_times, d.msd, **_style(key))
        low, high = d.fit_range
    if low is not None:
        axes.axvspan(low, high, color="grey", alpha=0.15, label="Einstein fit window")
    axes.set_xlabel("lag time (ps)")
    axes.set_ylabel(r"MSD (nm$^2$)")
    axes.set_title("Mean squared displacement")
    return _finish(axes, output)


def rotational_correlation(runs: Mapping[str, object], output: Path,
                           vector: str = "HH") -> Path:
    """C2(t) on a log axis, where a single-exponential decay is a straight line."""
    from .rotation import NOISE_FLOOR

    _figure, axes = plt.subplots(figsize=(6, 4))
    for key, summary in runs.items():
        c2 = summary.c2.get(vector)
        if c2 is None:
            continue
        lags = summary.rotation_lags[: len(c2)]
        # NaN, not a mask: a masked point is skipped and the line bridges the gap,
        # which on a log axis draws a cliff wherever C2 crosses zero.
        shown = np.where(c2 > 0, c2, np.nan)
        axes.semilogy(lags, shown, **_style(key))
    axes.axhline(NOISE_FLOOR, color="grey", linestyle=":", linewidth=1,
                 label="noise floor (integration stops)")
    axes.set_ylim(NOISE_FLOOR / 5, 1.2)
    axes.set_xlabel("time (ps)")
    axes.set_ylabel(rf"$C_2(t)$, {vector} vector")
    axes.set_title(f"Rotational correlation, {vector}")
    return _finish(axes, output)


def density_series(series: Mapping[str, object], output: Path) -> Path:
    """Density against time, against the experimental value."""
    _figure, axes = plt.subplots(figsize=(6, 4))
    for key, values in series.items():
        axes.plot(range(len(values)), values, alpha=0.8, **_style(key))
    axes.axhline(experiment.EXPERIMENT["density"].value, color="k",
                 linestyle=":", linewidth=1, label="experiment")
    axes.set_xlabel("frame")
    axes.set_ylabel(r"density (kg m$^{-3}$)")
    axes.set_title("Density")
    return _finish(axes, output)
