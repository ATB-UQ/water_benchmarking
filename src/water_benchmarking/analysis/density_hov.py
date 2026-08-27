"""Density and heat of vaporisation from the energy trajectory.

For a rigid, non-polarisable model the gas-phase molecule has no internal energy
and does not interact, so U_gas = 0 and no vacuum leg has to be simulated:

    dH_vap = -<U_pot>/N + RT

SPC/E is a special case.  Its charges are deliberately larger than a gas-phase
water's to mimic the polarisation of the liquid, and the energy cost of that
polarisation is not in the potential energy.  Berendsen's self-polarisation
correction (+5.22 kJ/mol subtracted from dH_vap) is what makes SPC/E's published
dH_vap comparable with experiment, so both numbers are reported.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .. import protocol
from . import errors

GAS_CONSTANT = 0.00831446         # kJ mol^-1 K^-1

#: Berendsen et al. 1987, J. Phys. Chem. 91:6269 -- the SPC/E self-polarisation term.
SPCE_POLARISATION_CORRECTION = 5.22  # kJ mol^-1


@dataclass
class ThermodynamicsResult:
    density: errors.Estimate           # kg m^-3
    potential_energy: errors.Estimate  # kJ mol^-1 per molecule
    hov: float                         # kJ mol^-1
    hov_error: float
    hov_polarisation_corrected: float | None
    pressure: errors.Estimate          # atm
    density_series: object = None      # per-frame series, for the convergence plot


def run_ene_ana(energy_files: Sequence[Path], properties: Sequence[str], work_dir: Path) -> dict:
    """Run gromos++ ene_ana and read back the per-frame series it writes.

    The library is chosen by parse test: each candidate is tried until ene_ana
    produces the requested series.  A wrong library does not read garbage, it
    fails outright, which is what makes the trial safe.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    failures = []
    for library in protocol.ENE_ANA_LIBS:
        command = [
            str(protocol.GROMOS_BIN / "ene_ana"),
            # ene_ana runs in work_dir (it writes <prop>.dat into the cwd), so
            # every input path has to be absolute.
            "@en_files", *[str(Path(f).resolve()) for f in energy_files],
            "@prop", *properties,
            "@library", str(library),
        ]
        result = subprocess.run(command, cwd=work_dir, capture_output=True, text=True)
        if result.returncode == 0 and all((work_dir / f"{p}.dat").exists() for p in properties):
            series = {}
            for prop in properties:
                data = np.loadtxt(work_dir / f"{prop}.dat", comments="#")
                series[prop] = data[:, 1] if data.ndim == 2 else data
            return series
        failures.append(f"{library.name}: {result.stderr.strip().splitlines()[-1] if result.stderr.strip() else 'no output'}")
        for prop in properties:
            (work_dir / f"{prop}.dat").unlink(missing_ok=True)
    raise RuntimeError(
        f"no ene_ana library parses {[Path(f).name for f in energy_files]}:\n  " + "\n  ".join(failures)
    )


def analyse(
    energy_files: Sequence[Path],
    model: str,
    work_dir: Path,
    n_molecules: int = protocol.N_WATERS,
    temperature: float = protocol.TEMPERATURE,
    discard: float = 0.1,
) -> ThermodynamicsResult:
    """Density, potential energy and heat of vaporisation for one model."""
    # One ene_ana call per file, and the leading fraction dropped from each: the
    # GROMOS production runs are independent replicates that each start from
    # freshly drawn velocities, so every one has its own short re-thermalisation.
    # Dropping the head of the concatenation would discard all of the first
    # replicate and none of the others.
    collected = {name: [] for name in ("densit", "totpot", "pressu")}
    for index, energy_file in enumerate(energy_files):
        series = run_ene_ana([energy_file], tuple(collected), work_dir / f"seg{index:02d}")
        for name, values in series.items():
            collected[name].append(errors.drop_equilibration(values, discard))
    series = {name: np.concatenate(values) for name, values in collected.items()}

    density = errors.block_average(series["densit"])
    per_molecule = series["totpot"] / n_molecules
    energy = errors.block_average(per_molecule)
    pressure = errors.block_average(series["pressu"])

    hov = -energy.mean + GAS_CONSTANT * temperature
    corrected = hov - SPCE_POLARISATION_CORRECTION if model == "spce" else None

    return ThermodynamicsResult(
        density_series=series["densit"],
        density=density,
        potential_energy=energy,
        hov=hov,
        hov_error=energy.error,
        hov_polarisation_corrected=corrected,
        pressure=pressure,
    )
