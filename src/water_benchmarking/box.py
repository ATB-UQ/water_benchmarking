"""Cut a 2048-molecule cubic water box out of the equilibrated 5384-molecule box.

Cutting beats building from scratch: the source box is 5 ns of NPT equilibration
(Schmid 2009), so the local structure is already right and only the box edge and
the density need to relax.  We keep the 2048 molecules closest to the box centre
and then set the edge from the source number density, so the cut region is neither
compressed nor stretched relative to the liquid it came from.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence

import numpy as np

from . import protocol

ATOM_NAMES = ("OW", "HW1", "HW2")

# Seam repair: delete a molecule when an O-O contact is closer than this.  Water's
# first O-O peak is at 0.28 nm, so 0.24 removes real overlaps without thinning the
# first solvation shell.
SEAM_CLASH = 0.24
WINDOW_GROWTH = 0.02          # nm, how fast the cut window grows when short
LOCAL_DENSITY_RADIUS = 0.45   # nm, neighbour shell used to spot crowded molecules
# SPC/SPC-E masses, needed for the centre of mass and for the density assert.
MASSES = np.array([15.99940, 1.00800, 1.00800])
MOLAR_MASS = float(MASSES.sum())
AVOGADRO = 6.02214076e23


@dataclass
class WaterBox:
    """Coordinates of whole water molecules in a cubic periodic box."""

    positions: np.ndarray  # (n_molecules, 3, 3) nm, atom order OW HW1 HW2
    edge: float            # nm

    @property
    def n_molecules(self) -> int:
        return self.positions.shape[0]

    @property
    def density(self) -> float:
        """kg m^-3."""
        mass_kg = self.n_molecules * MOLAR_MASS / AVOGADRO * 1e-3
        volume_m3 = (self.edge * 1e-9) ** 3
        return mass_kg / volume_m3


def read_g96_positions(path: Path) -> tuple[np.ndarray, float]:
    """Read the POSITION/POSITIONRED and GENBOX blocks of a g96 file.

    Returns (positions[n_molecules, 3, 3], edge).  Assumes a cubic box of whole
    3-site waters, which is what every solvent box in the library is.
    """
    coords: List[Sequence[float]] = []
    box_edge = None
    with open(path) as handle:
        block = None
        genbox_values: List[float] = []
        for line in handle:
            stripped = line.strip()
            if block is None:
                if stripped in ("POSITION", "POSITIONRED", "GENBOX"):
                    block = stripped
                continue
            if stripped == "END":
                if block == "GENBOX":
                    # GENBOX: ntb line, then a b c, then the angles.
                    box_edge = genbox_values[1:4]
                block = None
                continue
            if stripped.startswith("#"):
                continue
            if block in ("POSITION", "POSITIONRED"):
                # POSITION carries 24 chars of labels; POSITIONRED is bare numbers.
                fields = stripped.split()
                coords.append([float(x) for x in fields[-3:]])
            elif block == "GENBOX":
                genbox_values.extend(float(x) for x in stripped.split())

    if box_edge is None:
        raise ValueError(f"{path}: no GENBOX block")
    a, b, c = box_edge
    if not (math.isclose(a, b, rel_tol=1e-6) and math.isclose(a, c, rel_tol=1e-6)):
        raise ValueError(f"{path}: box is not cubic ({a}, {b}, {c})")

    xyz = np.asarray(coords, dtype=float)
    if xyz.shape[0] % 3:
        raise ValueError(f"{path}: {xyz.shape[0]} atoms is not a whole number of waters")
    return xyz.reshape(-1, 3, 3), float(a)


def _centres_of_mass(positions: np.ndarray) -> np.ndarray:
    return (positions * MASSES[None, :, None]).sum(axis=1) / MOLAR_MASS


def _seam_repaired_window(positions, com, source_edge, edge, n_target):
    """Take the cubic window [0, edge)^3 of the periodic source and heal its seam.

    A sub-cube of a periodic box is not itself periodic: the two cut surfaces that
    become neighbours under the new boundary are uncorrelated, so they overlap.
    Every choice of window offset has this problem (a scan over 216 offsets left
    430 contacts below 0.25 nm), so the seam is repaired rather than avoided --
    the overlapping molecules are deleted, exactly as gmx solvate does.
    """
    from collections import Counter

    from scipy.spatial import cKDTree

    wrapped = np.mod(com, source_edge)
    inside = np.where(np.all(wrapped < edge, axis=1))[0]
    # Shift each molecule as a rigid unit by the wrap applied to its centre of mass.
    # Shift whole molecules by the wrap their centre of mass needed; applying
    # np.mod to the atoms would split molecules across the periodic boundary.
    shift = com[inside] - wrapped[inside]
    selected = positions[inside] - shift[:, None, :]

    oxygens = np.mod(selected[:, 0, :], edge)
    pairs = cKDTree(oxygens, boxsize=edge).query_pairs(SEAM_CLASH, output_type="ndarray")
    alive = np.ones(len(inside), dtype=bool)
    if len(pairs):
        # Drop the more promiscuous partner first, so one deletion can settle
        # several clashes and the repair stays as small as possible.
        clash_count = Counter(pairs.ravel().tolist())
        for i, j in sorted(pairs.tolist(), key=lambda p: -(clash_count[p[0]] + clash_count[p[1]])):
            if alive[i] and alive[j]:
                alive[j if clash_count[i] < clash_count[j] else i] = False
    return selected[alive]


def _trim_to(selected: np.ndarray, edge: float, n_target: int) -> np.ndarray:
    """Discard the surplus from the densest local environments.

    Removing from crowded regions evens the density out; removing from sparse ones
    would deepen the cavities the seam repair just made.
    """
    from scipy.spatial import cKDTree

    surplus = len(selected) - n_target
    if surplus <= 0:
        return selected
    com = _centres_of_mass(selected)
    wrapped = np.mod(com, edge)
    neighbours = cKDTree(wrapped, boxsize=edge).query_ball_point(wrapped, LOCAL_DENSITY_RADIUS)
    crowded_first = np.argsort([-len(n) for n in neighbours], kind="stable")
    return selected[np.sort(crowded_first[surplus:])]


def cut(source: Path, n_waters: int) -> WaterBox:
    """Cut a clash-free n_waters cubic box out of the equilibrated source box.

    Grow the cut window until the seam repair still leaves enough molecules, trim
    to exactly n_waters, then rescale the molecular centres back to the source
    number density so the box starts at the liquid density it came from.  Only the
    centres move: the rigid geometry SHAKE will enforce is never distorted.
    """
    positions, source_edge = read_g96_positions(source)
    if n_waters > positions.shape[0]:
        raise ValueError(f"asked for {n_waters} of {positions.shape[0]} molecules")
    com = _centres_of_mass(positions)

    number_density = positions.shape[0] / source_edge**3
    target_edge = (n_waters / number_density) ** (1.0 / 3.0)

    edge = target_edge
    while edge < source_edge:
        selected = _seam_repaired_window(positions, com, source_edge, edge, n_waters)
        if len(selected) >= n_waters:
            break
        edge += WINDOW_GROWTH
    else:
        raise RuntimeError(f"could not find {n_waters} clash-free molecules in {source}")

    selected = _trim_to(selected, edge, n_waters)

    scale = target_edge / edge
    centres = _centres_of_mass(selected)
    selected = selected + (centres * scale - centres)[:, None, :]

    # Wrap by centre of mass, never per atom: wrapping atoms independently would
    # tear a molecule in half across the boundary and SHAKE fails at step 0.
    centres = _centres_of_mass(selected)
    selected = selected - (np.floor(centres / target_edge) * target_edge)[:, None, :]
    return WaterBox(positions=selected, edge=target_edge)


def min_image_oxygen_distance(box: WaterBox) -> float:
    """Smallest O-O distance under periodic boundaries, in nm."""
    oxygens = box.positions[:, 0, :]
    delta = oxygens[:, None, :] - oxygens[None, :, :]
    delta -= np.round(delta / box.edge) * box.edge
    dist = np.sqrt((delta**2).sum(axis=-1))
    np.fill_diagonal(dist, np.inf)
    return float(dist.min())


def write_cnf(box: WaterBox, path: Path, title: str) -> None:
    """Write a GROMOS .cnf: POSITION plus a cubic GENBOX."""
    lines = ["TITLE", title, "END", "POSITION", "# first 24 chars ignored"]
    atom = 0
    for mol in range(box.n_molecules):
        for site, name in enumerate(ATOM_NAMES):
            atom += 1
            x, y, z = box.positions[mol, site]
            lines.append(
                f"{mol + 1:5d} SOLV  {name:<5s}{atom:6d}"
                f"{x:15.9f}{y:15.9f}{z:15.9f}"
            )
    # md++ requires LATTICESHIFTS whenever the box is periodic; a fresh box has
    # not shifted anything yet, so every atom carries the zero shift.
    lines.append("END")
    lines.append("LATTICESHIFTS")
    lines.extend(["    0    0    0"] * (box.n_molecules * len(ATOM_NAMES)))
    lines += [
        "END",
        "GENBOX",
        "    1",
        f"{box.edge:15.9f}{box.edge:15.9f}{box.edge:15.9f}",
        f"{90.0:15.9f}{90.0:15.9f}{90.0:15.9f}",
        f"{0.0:15.9f}{0.0:15.9f}{0.0:15.9f}",
        f"{0.0:15.9f}{0.0:15.9f}{0.0:15.9f}",
        "END",
    ]
    path.write_text("\n".join(lines) + "\n")


def write_gro(box: WaterBox, path: Path, title: str) -> None:
    """Write a GROMACS .gro of the same box (six decimals, SOL residues)."""
    lines = [title, f"{box.n_molecules * 3:5d}"]
    atom = 0
    for mol in range(box.n_molecules):
        for site, name in enumerate(ATOM_NAMES):
            atom += 1
            x, y, z = box.positions[mol, site]
            lines.append(
                f"{(mol + 1) % 100000:5d}SOL  {name:>5s}{atom % 100000:5d}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
            )
    lines.append(f"{box.edge:10.5f}{box.edge:10.5f}{box.edge:10.5f}")
    path.write_text("\n".join(lines) + "\n")


def build(output_dir: Path, n_waters: int = protocol.N_WATERS) -> WaterBox:
    """Cut the box, check it, and write both engine formats."""
    box = cut(protocol.SOURCE_WATER_BOX, n_waters)

    if box.n_molecules != n_waters:
        raise AssertionError(f"cut {box.n_molecules} molecules, wanted {n_waters}")
    half = box.edge / 2.0
    if half <= protocol.CUTOFF:
        raise AssertionError(
            f"edge {box.edge:.3f} nm gives half-box {half:.3f} nm, "
            f"which does not clear the {protocol.CUTOFF} nm cutoff"
        )
    closest = min_image_oxygen_distance(box)
    if closest < 0.22:
        raise AssertionError(f"closest O-O is {closest:.3f} nm; box has a clash")

    output_dir.mkdir(parents=True, exist_ok=True)
    title = (
        f"{n_waters} water molecules, cubic box a = {box.edge:.9f} nm\n"
        f"cut from {protocol.SOURCE_WATER_BOX.name} (Schmid 2009, 5 ns NPT)"
    )
    write_cnf(box, output_dir / f"water_{n_waters}.cnf", title)
    write_gro(box, output_dir / f"water_{n_waters}.gro", title.replace("\n", "; "))
    return box
