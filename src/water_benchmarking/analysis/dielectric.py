"""Static dielectric permittivity from total dipole moment fluctuations.

The simulation uses a reaction field, so the plain Kirkwood relation does not
apply: the surrounding continuum of permittivity eps_rf reacts back on the box.
The correct expression (Neumann 1983) is

    (eps - 1)(2 eps_rf + 1) / (2 eps_rf + eps) = <M^2> - <M>^2 / (3 eps_0 V k_B T)

which is solved for eps below.  Using the vacuum formula instead would understate
eps by roughly 15% at eps_rf = 61, so this is not a detail.

Convergence is slow -- <M^2> is a whole-box quantity, so one box gives one sample
per frame no matter how many molecules it holds -- which is why the running
estimate is reported alongside the final number.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from .. import protocol
from . import errors

# Conversion: charges in e, distances in nm, so M is in e nm.
ELEMENTARY_CHARGE = 1.602176634e-19       # C
VACUUM_PERMITTIVITY = 8.8541878128e-12    # F m^-1
BOLTZMANN = 1.380649e-23                  # J K^-1
NM = 1e-9


@dataclass
class DielectricResult:
    epsilon: float
    epsilon_error: float
    mean_square_fluctuation: float        # e^2 nm^2
    running: np.ndarray                   # epsilon estimated from the first n frames
    times: np.ndarray


def total_dipole(positions: np.ndarray, charges: np.ndarray) -> np.ndarray:
    """Box dipole moment M = sum_i q_i r_i, in e nm.

    No gathering is needed and none must be done: the molecules are rigid and
    written whole, and M is only well defined for a neutral system when each
    molecule is intact.  Wrapping atoms individually would scatter charge across
    the box and destroy M.
    """
    return (positions * charges[None, :, None]).sum(axis=(0, 1))


def _solve_epsilon(fluctuation_si: float, eps_rf: float) -> float:
    """Invert the reaction-field relation for eps.

    With y the right-hand side, (eps-1)(2 eps_rf+1) = y (2 eps_rf + eps), so
    eps (2 eps_rf + 1 - y) = 2 eps_rf + 1 + 2 y eps_rf.
    """
    y = fluctuation_si
    denominator = (2.0 * eps_rf + 1.0) - y
    if denominator <= 0:
        return float("inf")
    return (2.0 * eps_rf + 1.0 + 2.0 * y * eps_rf) / denominator


def _fluctuation_to_si(mean_sq: float, volume_nm3: float, temperature: float) -> float:
    """(<M^2> - <M>^2) / (3 eps_0 V k_B T), dimensionless."""
    numerator = mean_sq * (ELEMENTARY_CHARGE * NM) ** 2
    denominator = 3.0 * VACUUM_PERMITTIVITY * (volume_nm3 * NM**3) * BOLTZMANN * temperature
    return numerator / denominator


def dielectric_constant(
    frames: Iterable,
    charges: Sequence[float],
    temperature: float = protocol.TEMPERATURE,
    eps_rf: float = protocol.EPSILON_RF,
    running_every: int = 100,
) -> DielectricResult:
    """Static permittivity from a trajectory."""
    charges = np.asarray(charges, dtype=float)
    dipoles = []
    volumes = []
    times = []
    for frame in frames:
        dipoles.append(total_dipole(frame.positions, charges))
        volumes.append(frame.volume)
        times.append(frame.time)

    dipoles_arr = np.asarray(dipoles)
    volumes_arr = np.asarray(volumes)
    times_arr = np.asarray(times)
    if len(dipoles_arr) < 2:
        raise ValueError("need at least two frames for a dielectric constant")

    def epsilon_from(slice_end: int) -> float:
        m = dipoles_arr[:slice_end]
        mean_sq = float((m**2).sum(axis=1).mean() - (m.mean(axis=0) ** 2).sum())
        volume = float(volumes_arr[:slice_end].mean())
        return _solve_epsilon(_fluctuation_to_si(mean_sq, volume, temperature), eps_rf)

    running_points = range(running_every, len(dipoles_arr) + 1, running_every)
    running = np.array([epsilon_from(n) for n in running_points])
    running_times = times_arr[[n - 1 for n in running_points]] if len(running) else times_arr[:0]

    mean_sq_total = float((dipoles_arr**2).sum(axis=1).mean() - (dipoles_arr.mean(axis=0) ** 2).sum())
    epsilon = epsilon_from(len(dipoles_arr))

    # Error from block averaging on M^2: the fluctuation, not the mean dipole, is
    # what carries the statistical uncertainty.
    m_squared = (dipoles_arr**2).sum(axis=1)
    block = errors.block_average(m_squared)
    volume = float(volumes_arr.mean())
    shifted = _solve_epsilon(
        _fluctuation_to_si(mean_sq_total + block.error, volume, temperature), eps_rf
    )
    return DielectricResult(
        epsilon=epsilon,
        epsilon_error=abs(shifted - epsilon),
        mean_square_fluctuation=mean_sq_total,
        running=running,
        times=running_times,
    )
