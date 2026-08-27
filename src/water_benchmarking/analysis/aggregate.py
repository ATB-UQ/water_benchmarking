"""Per-segment trajectory analysis, aggregated across a run.

Everything the trajectory can tell us is extracted in a **single pass per
segment**, for two reasons.

*Memory.* A 1 ns segment at 0.1 ps sampling is 10^4 frames of 2048 molecules;
holding the positions of all ten segments at once would be ~15 GB. Only the
reduced quantities are kept -- centres of mass, unit vectors, the box dipole --
and the frames themselves are discarded as they stream past.

*Correctness.* MSD and the rotational correlation functions are computed within a
segment and then averaged over segments, never across a segment boundary. The
boundary is a restart: correlations do carry across it physically, but the
100 ps fitting window is a tenth of a segment, so nothing is lost by staying
inside one -- and a lag that straddled a boundary would silently mix a duplicated
frame into the average.

The dielectric constant is the exception: it needs the longest possible series of
the box dipole, which is three numbers a frame, so those are simply concatenated.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .. import protocol
from ..box import MASSES, MOLAR_MASS
from . import dielectric, diffusion, rotation

#: Vectors whose reorientation experiment can see.
VECTORS = ("HH", "OH", "dipole")

#: Molecules sampled for the rotational correlation functions.  All 2048 would
#: cost four times the memory for an error bar that is already far below the
#: systematic uncertainty -- 512 molecules over 10^4 frames is ample statistics.
ROTATION_MOLECULES = 512


@dataclass
class SegmentSummary:
    """What one segment contributes, after the frames are gone."""

    lags: np.ndarray
    msd: np.ndarray
    rotation_lags: np.ndarray
    c1: dict = field(default_factory=dict)
    c2: dict = field(default_factory=dict)
    dipoles: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    volumes: np.ndarray = field(default_factory=lambda: np.empty(0))
    times: np.ndarray = field(default_factory=lambda: np.empty(0))
    mean_edge: float = 0.0


def analyse_segment(
    frames: Iterable,
    charges: Sequence[float],
    # 50 ps, not 20: tau_1 of the dipole is 3-8 ps and is taken as the integral
    # of C_1, which at 20 ps still had a few percent of its area outstanding.
    max_lag: float = 50.0,
    rotation_molecules: int = ROTATION_MOLECULES,
    check_whole: bool = True,
) -> SegmentSummary:
    """Stream one segment, returning everything downstream needs from it."""
    from ..gromacs import assert_whole_molecules

    charges = np.asarray(charges, dtype=float)
    centres, vectors, dipoles, volumes, times = [], {v: [] for v in VECTORS}, [], [], []
    selection = None

    for index, frame in enumerate(frames):
        if check_whole and index == 0:
            # One check per segment: a conversion is either right or wrong for the
            # whole file, and the test is not free.
            assert_whole_molecules(frame)
        if selection is None:
            count = min(rotation_molecules, frame.positions.shape[0])
            # A fixed stride, not a random draw: reproducible, and it samples the
            # box uniformly rather than clustering.
            selection = np.arange(0, frame.positions.shape[0],
                                  max(1, frame.positions.shape[0] // count))[:count]

        centres.append(
            ((frame.positions * MASSES[None, :, None]).sum(axis=1) / MOLAR_MASS
             ).astype(np.float32)
        )
        unit = rotation.molecular_vectors(frame.positions[selection])
        for name in VECTORS:
            vectors[name].append(unit[name].astype(np.float32))
        dipoles.append(dielectric.total_dipole(frame.positions, charges))
        volumes.append(frame.volume)
        times.append(frame.time)

    if len(times) < 10:
        raise ValueError(f"segment has only {len(times)} frames")

    centres_arr = np.asarray(centres)
    times_arr = np.asarray(times)
    edges = np.cbrt(np.asarray(volumes))
    lags = times_arr - times_arr[0]

    path = diffusion.unwrap(centres_arr, edges)
    msd = diffusion._msd_fft(path)

    dt = float(np.median(np.diff(times_arr)))
    max_lag_frames = min(len(times_arr), int(max_lag / dt) + 1)
    c1, c2 = {}, {}
    for name in VECTORS:
        stack = np.asarray(vectors[name])
        c1[name], c2[name] = rotation.correlation_functions(stack, max_lag_frames)

    return SegmentSummary(
        lags=lags,
        msd=msd,
        rotation_lags=np.arange(max_lag_frames) * dt,
        c1=c1,
        c2=c2,
        dipoles=np.asarray(dipoles),
        volumes=np.asarray(volumes),
        times=times_arr,
        mean_edge=float(edges.mean()),
    )


def _stack_and_average(arrays: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Mean and standard error across segments, truncated to the shortest."""
    length = min(len(a) for a in arrays)
    stacked = np.vstack([a[:length] for a in arrays])
    mean = stacked.mean(axis=0)
    error = stacked.std(axis=0, ddof=1) / np.sqrt(len(arrays)) if len(arrays) > 1 else \
        np.zeros_like(mean)
    return mean, error


@dataclass
class RunSummary:
    diffusion: diffusion.DiffusionResult
    dielectric: dielectric.DielectricResult
    tau1: dict
    tau2: dict
    rotation_lags: np.ndarray
    c2: dict
    n_segments: int


def analyse_run(
    segments: Sequence,
    charges: Sequence[float],
    n_molecules: int = protocol.N_WATERS,
    temperature: float = protocol.TEMPERATURE,
) -> RunSummary:
    """Analyse every segment of one run and combine them.

    `segments` is a sequence of either paths or zero-argument callables returning
    an iterable of frames.  The callable form exists so the GROMACS path can
    convert a .xtc, stream it, and delete the 2.8 GB text file again before the
    next one is made -- converting all ten up front would need 28 GB per model.
    """
    from ..trc import read_frames

    summaries = []
    for segment in segments:
        source = segment if callable(segment) else (
            lambda path=segment: read_frames(path, n_molecules)
        )
        summaries.append(analyse_segment(source(), charges))
    if not summaries:
        raise ValueError("no trajectories to analyse")

    msd, _msd_error = _stack_and_average([s.msd for s in summaries])
    lags = summaries[0].lags[: len(msd)]
    edge = float(np.mean([s.mean_edge for s in summaries]))
    d = diffusion.diffusion_from_msd(lags, msd, edge, temperature=temperature)

    rotation_lags = summaries[0].rotation_lags
    tau1, tau2, c2_mean = {}, {}, {}
    for name in VECTORS:
        c1, _ = _stack_and_average([s.c1[name] for s in summaries])
        c2, _ = _stack_and_average([s.c2[name] for s in summaries])
        lag = rotation_lags[: len(c2)]
        c2_mean[name] = c2
        tau1[name] = float(np.trapezoid(c1, lag))
        tau2[name] = float(np.trapezoid(c2, lag))

    # The dielectric constant wants one long series, not ten short ones; the box
    # dipole was kept during the streaming pass so nothing is re-read here.
    eps = dielectric.from_dipoles(
        np.vstack([s.dipoles for s in summaries]),
        np.concatenate([s.volumes for s in summaries]),
        np.concatenate([s.times for s in summaries]),
        temperature=temperature,
        running_every=500,
    )

    return RunSummary(
        diffusion=d, dielectric=eps, tau1=tau1, tau2=tau2,
        rotation_lags=rotation_lags[: len(next(iter(c2_mean.values())))],
        c2=c2_mean, n_segments=len(summaries),
    )
