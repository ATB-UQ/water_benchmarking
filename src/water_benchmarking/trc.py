"""Streaming reader for GROMOS position trajectories.

One reader serves both engines: GROMACS .xtc is converted to .g96 with
`gmx trjconv`, which writes the same POSITIONRED/GENBOX blocks, so diffusion,
rotational relaxation and the dielectric constant are computed by identical code
for GROMOS and GROMACS.  Any difference in the results is then a difference in the
simulation, not in the analysis -- which is the only way an engine comparison
means anything.

Trajectories are large (0.1 ps sampling over 10 ns is ~10^5 frames), so frames are
yielded one at a time and never all held in memory.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from . import protocol


@dataclass
class Frame:
    """One trajectory frame: whole molecules, plus the box that framed them."""

    time: float                 # ps
    step: int
    positions: np.ndarray       # (n_molecules, 3, 3) nm, atom order OW HW1 HW2
    edge: float                 # nm, cubic box

    @property
    def volume(self) -> float:
        return self.edge**3


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_frames(path: Path, n_molecules: int = protocol.N_WATERS) -> Iterator[Frame]:
    """Yield every frame of one trajectory file."""
    n_atoms = n_molecules * 3
    with _open(path) as handle:
        block = None
        coords: list[list[float]] = []
        genbox: list[float] = []
        time = 0.0
        step = 0
        for line in handle:
            stripped = line.strip()
            if block is None:
                if stripped in ("TIMESTEP", "POSITIONRED", "POSITION", "GENBOX"):
                    block = stripped
                    if block == "GENBOX":
                        genbox = []
                continue
            if stripped == "END":
                if block == "GENBOX":
                    # GENBOX is ntb, then the three edges, then angles and origin.
                    edge = genbox[1]
                    positions = np.asarray(coords, dtype=float).reshape(n_molecules, 3, 3)
                    yield Frame(time=time, step=step, positions=positions, edge=edge)
                    coords = []
                block = None
                continue
            if stripped.startswith("#") or not stripped:
                continue
            fields = stripped.split()
            if block == "TIMESTEP":
                step, time = int(fields[0]), float(fields[1])
            elif block in ("POSITIONRED", "POSITION"):
                coords.append([float(x) for x in fields[-3:]])
            elif block == "GENBOX":
                genbox.extend(float(x) for x in fields)

        if coords and len(coords) == n_atoms:
            raise ValueError(f"{path}: trailing frame with no GENBOX block")


def read_all(paths: Iterable[Path], n_molecules: int = protocol.N_WATERS) -> Iterator[Frame]:
    """Yield frames from several segments in order.

    Consecutive segments repeat the configuration at the segment boundary (a
    segment starts from the previous one's final structure), so the duplicate is
    dropped -- otherwise every boundary contributes a zero-displacement step and
    biases the diffusion coefficient low.
    """
    previous_time = None
    for path in paths:
        for frame in read_frames(path, n_molecules):
            if previous_time is not None and frame.time <= previous_time:
                continue
            previous_time = frame.time
            yield frame
