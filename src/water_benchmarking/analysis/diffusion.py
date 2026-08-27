"""Self-diffusion coefficient from the mean squared displacement.

Two things are easy to get wrong here and both are handled explicitly:

* **Unwrapping.**  Trajectory frames hold molecules wrapped into the box, so a
  molecule crossing the boundary appears to jump a whole box length.  Displacements
  are therefore taken between consecutive frames under the minimum image
  convention and accumulated, which reconstructs the continuous path.
* **Finite size.**  D from a periodic box is systematically low, because a molecule
  drags its own periodic images through the solvent.  The Yeh-Hummer correction
  D_inf = D_pbc + 2.837 k_B T / (6 pi eta L) removes the leading term; it is about
  +10% for a 4 nm box, so an uncorrected D cannot be compared with experiment.
  Both numbers are reported so the correction is never hidden.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .. import protocol

BOLTZMANN = 1.380649e-23          # J K^-1
YEH_HUMMER_XI = 2.837297          # the cubic-lattice self-interaction constant
WATER_VISCOSITY = 0.89e-3         # Pa s, experimental at 298.15 K


@dataclass
class DiffusionResult:
    d_pbc: float                  # m^2 s^-1, as simulated
    d_corrected: float            # m^2 s^-1, Yeh-Hummer corrected
    d_error: float
    lag_times: np.ndarray         # ps
    msd: np.ndarray               # nm^2
    fit_range: tuple[float, float]
    #: log-log slope of the MSD across the fit window; 1.00 is diffusive.  A value
    #: below ~0.95 means the window reaches into the sub-diffusive cage regime.
    linearity: float = float("nan")


def unwrap(centres: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Turn wrapped centres of mass into a continuous path.

    centres is (n_frames, n_molecules, 3); edges is (n_frames,).  The box breathes
    under NPT, so each step uses its own edge.
    """
    steps = np.diff(centres, axis=0)
    box = edges[1:, None, None]
    steps -= np.round(steps / box) * box
    path = np.empty_like(centres)
    path[0] = centres[0]
    path[1:] = centres[0] + np.cumsum(steps, axis=0)
    return path


def _msd_fft(path: np.ndarray) -> np.ndarray:
    """MSD averaged over all time origins, via the Fast Correlation Algorithm.

    Direct evaluation is O(n_frames^2) per molecule, which is hopeless for 10^4
    frames; this is O(n log n).  MSD(k) = S1(k) - 2 S2(k), where S1 comes from
    cumulative sums of squares and S2 is the position autocorrelation.
    """
    n_frames, n_molecules, _ = path.shape
    squares = (path**2).sum(axis=2)                      # (frames, molecules)

    total = squares.sum(axis=0)
    s1 = np.empty((n_frames, n_molecules))
    forward = np.concatenate([[np.zeros(n_molecules)], np.cumsum(squares, axis=0)])
    for lag in range(n_frames):
        # sum over origins of |r(t)|^2 + |r(t+lag)|^2
        head = forward[n_frames - lag] - forward[0]
        tail = total - (forward[lag] - forward[0])
        s1[lag] = (head + tail) / (n_frames - lag)

    size = 2 * n_frames
    fft = np.fft.rfft(path, n=size, axis=0)
    acf = np.fft.irfft(fft * np.conjugate(fft), n=size, axis=0)[:n_frames]
    s2 = acf.sum(axis=2) / (n_frames - np.arange(n_frames))[:, None]

    return (s1 - 2 * s2).mean(axis=1)


def mean_squared_displacement(frames: Iterable) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (lag_times_ps, msd_nm2, mean_edge_nm) for one trajectory."""
    from ..box import MASSES, MOLAR_MASS

    centres = []
    edges = []
    times = []
    for frame in frames:
        centres.append((frame.positions * MASSES[None, :, None]).sum(axis=1) / MOLAR_MASS)
        edges.append(frame.edge)
        times.append(frame.time)

    centres_arr = np.asarray(centres)
    edges_arr = np.asarray(edges)
    times_arr = np.asarray(times)
    if len(centres_arr) < 10:
        raise ValueError("need at least 10 frames for an MSD")

    path = unwrap(centres_arr, edges_arr)
    msd = _msd_fft(path)
    lags = times_arr - times_arr[0]
    return lags, msd, float(edges_arr.mean())


#: Fit window as a fraction of the longest lag.  The lower bound clears the
#: ballistic and cage-rattling regime by a wide margin; the upper bound stops
#: where the number of time origins per lag has fallen to half and the MSD
#: starts to get noisy.  For a 1 ns segment this is 100-500 ps.
FIT_WINDOW = (0.10, 0.50)


def diffusion_from_msd(
    lag_times: np.ndarray,
    msd: np.ndarray,
    edge: float,
    fit_from: float | None = None,
    fit_to: float | None = None,
    temperature: float = protocol.TEMPERATURE,
    viscosity: float = WATER_VISCOSITY,
) -> DiffusionResult:
    """Einstein fit of the MSD, plus the finite-size correction.

    The window defaults to a fixed fraction of the available lag range rather
    than fixed times, so it scales with the segment length and is not
    accidentally a tenth of the data.
    """
    longest = float(lag_times[-1])
    if fit_from is None:
        fit_from = FIT_WINDOW[0] * longest
    if fit_to is None:
        fit_to = FIT_WINDOW[1] * longest
    window = (lag_times >= fit_from) & (lag_times <= fit_to)
    if window.sum() < 3:
        raise ValueError(
            f"only {window.sum()} points between {fit_from} and {fit_to} ps; "
            "trajectory is too short or too coarsely sampled"
        )
    slope, _intercept = np.polyfit(lag_times[window], msd[window], 1)
    positive = window & (lag_times > 0) & (msd > 0)
    linearity = float(np.polyfit(np.log(lag_times[positive]), np.log(msd[positive]), 1)[0])
    # MSD in nm^2 vs ps; 6D = slope for three dimensions.
    d_pbc = slope / 6.0 * 1e-18 / 1e-12          # nm^2/ps -> m^2/s

    correction = (
        YEH_HUMMER_XI * BOLTZMANN * temperature
        / (6.0 * np.pi * viscosity * edge * 1e-9)
    )

    # Uncertainty from splitting the fit window in half: if the MSD is not yet
    # linear the two halves disagree, which is exactly the error that matters.
    midpoint = (fit_from + fit_to) / 2.0
    halves = []
    for low, high in ((fit_from, midpoint), (midpoint, fit_to)):
        sub = (lag_times >= low) & (lag_times <= high)
        if sub.sum() >= 3:
            halves.append(np.polyfit(lag_times[sub], msd[sub], 1)[0] / 6.0 * 1e-6)
    error = abs(halves[0] - halves[1]) / 2.0 if len(halves) == 2 else float("nan")

    return DiffusionResult(
        d_pbc=d_pbc,
        d_corrected=d_pbc + correction,
        d_error=error,
        lag_times=lag_times,
        msd=msd,
        fit_range=(fit_from, fit_to),
        linearity=linearity,
    )
