"""Build the pure-solvent GROMOS topologies for the water models under test.

make_top is called without @seq: an empty sequence is rejected outright
("Cannot find building block for"), but omitting it entirely yields exactly what
this benchmark needs -- SOLUTEATOM with NRP 0 and the requested solvent block.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from . import protocol


def build_topology(model: str, output_dir: Path) -> Path:
    """Write <model>.top for one water model and return its path."""
    entry = protocol.model(model)
    if entry.gromos_block is None:
        raise ValueError(
            f"{model!r} has no 54A7 solvent building block, so no GROMOS topology "
            f"can be built for it; it is run with {'/'.join(entry.engines)} only"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"{model}.top"

    result = subprocess.run(
        [
            str(protocol.GROMOS_BIN / "make_top"),
            "@build", str(protocol.FORCEFIELD_MTB),
            "@param", str(protocol.FORCEFIELD_IFP),
            "@solv", entry.gromos_block,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    destination.write_text(result.stdout)
    verify_topology(destination, model)
    return destination


def solvent_charges(topology: Path) -> list[float]:
    """Charges of the three solvent sites, read back out of a topology."""
    lines = topology.read_text().splitlines()
    start = lines.index("SOLVENTATOM")
    charges = []
    for line in lines[start + 1:]:
        if line.strip() == "END":
            break
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        if len(fields) == 5:
            charges.append(float(fields[-1]))
    return charges


#: Published charges, used to prove make_top picked the block we asked for -- and,
#: for the GROMACS path, to weight the box dipole (cli.analyse_run -> dielectric).
#: Held in the model registry so the two engines cannot disagree about a model.
EXPECTED_CHARGES = {name: entry.charges for name, entry in protocol.MODELS.items()}


def verify_topology(topology: Path, model: str) -> None:
    charges = solvent_charges(topology)
    expected = EXPECTED_CHARGES[model]
    if len(charges) != 3:
        raise AssertionError(f"{topology}: found {len(charges)} solvent sites, expected 3")
    for got, want in zip(charges, expected):
        if abs(got - want) > 1e-6:
            raise AssertionError(
                f"{topology}: solvent charges {charges} are not {model.upper()} {list(expected)}"
            )
    if "SOLUTEATOM" not in topology.read_text():
        raise AssertionError(f"{topology}: no SOLUTEATOM block")
