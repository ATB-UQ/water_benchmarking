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

import math
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

#: Gas-phase dipole and polarisability behind the self-polarisation correction
#: (Berendsen et al. 1987): a model whose dipole is enhanced over the gas-phase
#: value pays a self-energy to hold that enhancement, and it must be given back
#: before dH_vap is compared with experiment.
GAS_DIPOLE = 1.85            # D
POLARISABILITY = 1.44        # Angstrom^3
_DEBYE = 3.33564e-30         # C m
_EPS0 = 8.8541878128e-12
_AVOGADRO = 6.02214076e23


def dipole_moment(model: str) -> float:
    """The model's dipole in Debye, from its charges and geometry.

    mu = 2 q_H r_OH cos(theta/2).  Derived rather than tabulated so it cannot
    disagree with the parameters actually simulated; reproduces the values in
    Table I of Izadi & Onufriev 2016 (SPC/E 2.35, OPC3 2.43).
    """
    from .. import protocol

    entry = protocol.model(model)
    theta = 2 * math.asin(entry.r_hh / (2 * entry.r_oh))
    # e.Angstrom -> Debye; r is held in nm.
    return 2 * entry.charges[1] * (entry.r_oh * 10) * math.cos(theta / 2) * 4.803205


def self_polarisation_energy(model: str) -> float:
    """(mu - mu_gas)^2 / (2 alpha), in kJ mol^-1."""
    delta = (dipole_moment(model) - GAS_DIPOLE) * _DEBYE
    alpha_si = 4 * math.pi * _EPS0 * POLARISABILITY * 1e-30
    return delta ** 2 / (2 * alpha_si) * _AVOGADRO / 1000.0


#: Self-polarisation correction subtracted from dH_vap, per model, kJ mol^-1.
#: Opt-in per model rather than applied wherever the dipole is enhanced: it is
#: conventionally applied to the models parameterised with it in mind, and the
#: published SPC numbers this benchmark is checked against are uncorrected (SPC
#: would take 3.76 by the same formula).
#:
#: OPC3 is corrected because Izadi & Onufriev 2016 correct theirs -- their SPC/E
#: entry is 10.43 kcal/mol (43.64 kJ/mol), which is the corrected value, not the
#: ~49 kJ/mol a raw SPC/E run gives.  Both models then land within 0.4 kJ/mol of
#: the paper: this protocol gives 49.27 - 5.24 = 44.03 against their 43.64, and
#: 51.66 - 7.03 = 44.62 against their 44.89.
POLARISATION_CORRECTION = {
    "spce": SPCE_POLARISATION_CORRECTION,   # 5.22 published; the formula gives 5.24
    "opc3": 7.03,
}


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
    polarisation = POLARISATION_CORRECTION.get(model)
    corrected = hov - polarisation if polarisation is not None else None

    return ThermodynamicsResult(
        density_series=series["densit"],
        density=density,
        potential_energy=energy,
        hov=hov,
        hov_error=energy.error,
        hov_polarisation_corrected=corrected,
        pressure=pressure,
    )
