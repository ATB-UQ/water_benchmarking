"""Rotational relaxation of water, as measured by NMR and dielectric relaxation.

C_l(t) = <P_l(u(0) . u(t))> for the molecular vectors that experiment can see:

* the H-H vector and the O-H vector, whose l = 2 correlation times are what NMR
  relaxation of D2O and H2O report;
* the molecular dipole (the HOH bisector), whose l = 1 correlation time is the
  microscopic counterpart of the Debye relaxation time.

Both orders are computed by FFT.  C_1 is the autocorrelation of the three
components of u.  C_2 needs <(u(0).u(t))^2>, which is not an autocorrelation of u
but of the second-rank tensor Q = u (x) u -- expanding the square gives
sum_ab <Q_ab(0) Q_ab(t)>, so six independent autocorrelations (three diagonal,
three off-diagonal counted twice) reconstruct it exactly.  Doing it directly over
time-origin pairs would be O(n_frames^2) and is not affordable at 0.1 ps sampling.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

# Off-diagonal tensor components appear twice in the sum over a and b.
_TENSOR_WEIGHTS = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])


@dataclass
class RotationResult:
    lag_times: np.ndarray         # ps
    c1: np.ndarray
    c2: np.ndarray
    tau1: float                   # ps
    tau2: float                   # ps
    tau2_fit: float               # ps, from an exponential fit of the tail


def molecular_vectors(positions: np.ndarray) -> dict[str, np.ndarray]:
    """Unit vectors of one frame, keyed by name.

    positions is (n_molecules, 3, 3) in the order OW, HW1, HW2.  The molecules are
    rigid and whole, so no periodic correction is needed inside a molecule.
    """
    oxygen = positions[:, 0, :]
    h1 = positions[:, 1, :]
    h2 = positions[:, 2, :]

    vectors = {
        "OH": h1 - oxygen,
        "HH": h2 - h1,
        # The dipole of a symmetric 3-site water lies along the HOH bisector.
        "dipole": (h1 + h2) / 2.0 - oxygen,
    }
    return {
        name: vector / np.linalg.norm(vector, axis=1, keepdims=True)
        for name, vector in vectors.items()
    }


def _autocorrelation(series: np.ndarray) -> np.ndarray:
    """Autocorrelation along axis 0, averaged over all time origins.

    series is (n_frames, n_molecules, n_components); returns (n_frames,) summed
    over components and averaged over molecules.
    """
    n_frames = series.shape[0]
    size = 2 * n_frames
    fft = np.fft.rfft(series, n=size, axis=0)
    acf = np.fft.irfft(fft * np.conjugate(fft), n=size, axis=0)[:n_frames]
    normalisation = (n_frames - np.arange(n_frames))[:, None, None]
    return (acf / normalisation).mean(axis=1)


def correlation_functions(
    vectors: np.ndarray, max_lag_frames: int
) -> tuple[np.ndarray, np.ndarray]:
    """C_1 and C_2 for a stack of unit vectors (n_frames, n_molecules, 3)."""
    c1 = _autocorrelation(vectors).sum(axis=1)[:max_lag_frames]

    x, y, z = vectors[..., 0], vectors[..., 1], vectors[..., 2]
    tensor = np.stack([x * x, y * y, z * z, x * y, x * z, y * z], axis=-1)
    squared_dot = (_autocorrelation(tensor) * _TENSOR_WEIGHTS[None, :]).sum(axis=1)
    c2 = (3.0 * squared_dot[:max_lag_frames] - 1.0) / 2.0
    return c1, c2


#: Below this C(t) is indistinguishable from the statistical floor and must not
#: be integrated -- it is noise around zero, and on a log plot it is the
#: cliff-and-plateau garbage that looks like a broken analysis.
NOISE_FLOOR = 2e-3


def correlation_time(lags: np.ndarray, c: np.ndarray,
                     fit_from: float = 1.0, fit_to: float | None = None) -> tuple[float, float]:
    """tau = integral of C(t) up to the noise floor, plus the exponential tail.

    Returns (tau, tau_fit).  C(t) is integrated only while it is clearly above
    the floor; the remainder is added analytically as tau_fit * C(t_cut), using
    the decay time fitted on the last decade above the floor.  Integrating the
    noisy tail instead adds a random walk of area that does not average out
    over any practical run length.
    """
    above = c > NOISE_FLOOR
    cut = int(np.argmin(above)) if not above.all() else len(c)
    cut = max(cut, 3)
    t_cut = lags[cut - 1]

    if fit_to is None:
        fit_to = t_cut
    window = (lags >= fit_from) & (lags <= fit_to) & (c > NOISE_FLOOR)
    if window.sum() >= 3:
        slope = np.polyfit(lags[window], np.log(c[window]), 1)[0]
        tau_fit = float(-1.0 / slope) if slope < 0 else float("nan")
    else:
        tau_fit = float("nan")

    area = float(np.trapezoid(c[:cut], lags[:cut]))
    tail = tau_fit * float(c[cut - 1]) if tau_fit == tau_fit else 0.0
    return area + tail, tau_fit


def analyse(
    frames: Iterable,
    vector: str = "HH",
    max_lag: float = 20.0,
    fit_from: float = 1.0,
    fit_to: float = 10.0,
) -> RotationResult:
    """Correlation functions and correlation times for one molecular vector."""
    stack = []
    times = []
    for frame in frames:
        stack.append(molecular_vectors(frame.positions)[vector].astype(np.float32))
        times.append(frame.time)

    vectors = np.asarray(stack)
    times_arr = np.asarray(times)
    dt = float(np.median(np.diff(times_arr)))
    max_lag_frames = min(len(vectors), int(max_lag / dt) + 1)

    c1, c2 = correlation_functions(vectors, max_lag_frames)
    lags = np.arange(max_lag_frames) * dt

    tau1, _ = correlation_time(lags, c1, fit_from)
    tau2, tau2_fit = correlation_time(lags, c2, fit_from)

    return RotationResult(
        lag_times=lags, c1=c1, c2=c2, tau1=tau1, tau2=tau2, tau2_fit=tau2_fit
    )
