"""Static dielectric permittivity from total dipole moment fluctuations.

The measured quantity is the dimensionless fluctuation

    y = (<M^2> - <M>^2) / (3 eps_0 V k_B T)

and the question is which relation turns y into eps.  For a conducting boundary
(Ewald, or a reaction field with eps_rf = infinity) it is eps = 1 + y.  For a
reaction field with finite eps_rf, Neumann (1983) gives

    (eps - 1)(2 eps_rf + 1) / (2 eps_rf + eps) = y

which for eps_rf = 61 has a pole at y = 2 eps_rf + 1 = 123.  Liquid water sits
at y ~ 65-75, within fifty units of that pole, where d(eps)/dy is about 3: a
ten percent change in y moves eps by fifty.

That is not merely a precision problem.  Measured on this system, y is the same
to within its statistical scatter whether the reaction field is applied at 0.9,
1.4 or 1.8 nm, in a 2048- or 16384-molecule box, or replaced by PME altogether:
every geometry gives y = 59-76, and the 8x box under the reaction field gives
the PME value exactly.  Neumann's relation presumes the fluctuation carries a
strong dependence on eps_rf, but the force field it describes barely does --
k_rf at eps_rf = 61 is within 2.4 % of the conducting-boundary value, so the
box cannot "know" which it is under to the precision the formula demands.

The conducting-boundary relation is therefore the one reported, with the
Neumann value and its sensitivity alongside so the choice is never silent.
Convergence is slow regardless -- <M^2> is a whole-box quantity, so one box
gives one sample per frame -- which is why the running estimate is kept.
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
    epsilon: float                        # conducting-boundary relation, 1 + y
    epsilon_error: float
    y: float                              # the dimensionless fluctuation itself
    epsilon_neumann: float                # Neumann relation at eps_rf, for the record
    neumann_sensitivity: float            # d(eps_neumann)/dy: how unstable that number is
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


def _conducting_boundary_epsilon(fluctuation_si: float) -> float:
    """eps = 1 + y: the relation for Ewald or a conducting reaction field."""
    return 1.0 + fluctuation_si


def _neumann_sensitivity(fluctuation_si: float, eps_rf: float) -> float:
    """d(eps)/dy for the Neumann relation -- large near its pole at 2 eps_rf + 1."""
    denominator = (2.0 * eps_rf + 1.0) - fluctuation_si
    if denominator <= 0:
        return float("inf")
    return (2.0 * eps_rf + 1.0) * (2.0 * eps_rf + 1.0) / denominator**2


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


def from_dipoles(
    dipoles: np.ndarray,
    volumes: np.ndarray,
    times: np.ndarray,
    temperature: float = protocol.TEMPERATURE,
    eps_rf: float = protocol.EPSILON_RF,
    running_every: int = 100,
) -> DielectricResult:
    """Static permittivity from an already-extracted series of box dipoles.

    This is the real entry point: the box dipole is three numbers a frame, so a
    streaming pass over the trajectory can keep the whole series cheaply while the
    positions it came from are discarded.
    """
    dipoles = np.asarray(dipoles, dtype=float)
    volumes = np.asarray(volumes, dtype=float)
    times = np.asarray(times, dtype=float)
    if len(dipoles) < 2:
        raise ValueError("need at least two frames for a dielectric constant")

    def fluctuation_from(end: int) -> float:
        m = dipoles[:end]
        mean_sq = float((m**2).sum(axis=1).mean() - (m.mean(axis=0) ** 2).sum())
        return _fluctuation_to_si(mean_sq, float(volumes[:end].mean()), temperature)

    points = range(running_every, len(dipoles) + 1, running_every)
    running = np.array([_conducting_boundary_epsilon(fluctuation_from(n)) for n in points])
    running_times = times[[n - 1 for n in points]] if len(running) else times[:0]

    mean_sq_total = float(
        (dipoles**2).sum(axis=1).mean() - (dipoles.mean(axis=0) ** 2).sum()
    )
    y_total = fluctuation_from(len(dipoles))
    epsilon = _conducting_boundary_epsilon(y_total)

    # The uncertainty lives in the fluctuation, not the mean dipole, so the block
    # average is taken on M^2 and propagated through the relation.
    block = errors.block_average((dipoles**2).sum(axis=1))
    y_shifted = _fluctuation_to_si(mean_sq_total + block.error, float(volumes.mean()), temperature)
    return DielectricResult(
        epsilon=epsilon,
        epsilon_error=abs(_conducting_boundary_epsilon(y_shifted) - epsilon),
        y=y_total,
        epsilon_neumann=_solve_epsilon(y_total, eps_rf),
        neumann_sensitivity=_neumann_sensitivity(y_total, eps_rf),
        mean_square_fluctuation=mean_sq_total,
        running=running,
        times=running_times,
    )


def dielectric_constant(
    frames: Iterable,
    charges: Sequence[float],
    temperature: float = protocol.TEMPERATURE,
    eps_rf: float = protocol.EPSILON_RF,
    running_every: int = 100,
) -> DielectricResult:
    """Static permittivity straight from a trajectory."""
    charges = np.asarray(charges, dtype=float)
    dipoles, volumes, times = [], [], []
    for frame in frames:
        dipoles.append(total_dipole(frame.positions, charges))
        volumes.append(frame.volume)
        times.append(frame.time)
    return from_dipoles(
        np.asarray(dipoles), np.asarray(volumes), np.asarray(times),
        temperature=temperature, eps_rf=eps_rf, running_every=running_every,
    )
