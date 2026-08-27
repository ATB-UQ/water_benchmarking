"""Experimental reference values for liquid water at 298.15 K and 1 atm.

Also carries the literature values for the two models, which serve a different
purpose: experiment says whether the model is good, the literature values say
whether *this run* reproduced the model.  A result far from both is a setup bug,
not a finding.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Reference:
    value: float
    unit: str
    source: str


EXPERIMENT = {
    "density": Reference(997.05, "kg m^-3", "Kell 1975, J. Chem. Eng. Data 20:97"),
    "hov": Reference(43.99, "kJ mol^-1", "Wagner & Pruss 2002, IAPWS-95"),
    "diffusion": Reference(2.30e-9, "m^2 s^-1", "Holz et al. 2000, PCCP 2:4740"),
    "dielectric": Reference(78.4, "-", "Fernandez et al. 1997, J. Phys. Chem. Ref. Data 26:1125"),
    "tau2_HH": Reference(2.0, "ps", "NMR relaxation; Ludwig 2001, Angew. Chem. 40:1808"),
    "tau2_OH": Reference(1.95, "ps", "NMR relaxation; Ludwig 2001"),
    "tau1_dipole": Reference(8.3, "ps", "Debye relaxation; Ronne et al. 1997, JCP 107:5319"),
}

#: Published values for these models, as a check that the run reproduces the model.
LITERATURE = {
    "spc": {
        "density": (972.0, 985.0),
        "hov": (41.0, 44.0),
        "diffusion": (3.9e-9, 4.3e-9),
        "dielectric": (62.0, 68.0),
        "tau2_HH": (0.9, 1.3),
    },
    "spce": {
        "density": (994.0, 1001.0),
        "hov": (46.0, 49.5),          # before the polarisation correction
        "diffusion": (2.4e-9, 2.8e-9),
        "dielectric": (68.0, 74.0),
        "tau2_HH": (1.7, 2.2),
    },
}


def deviation(simulated: float, key: str) -> float:
    """Percentage deviation of a simulated value from experiment."""
    reference = EXPERIMENT[key].value
    return 100.0 * (simulated - reference) / reference


def within_literature(simulated: float, model: str, key: str) -> bool | None:
    """True/False if a literature range exists for this model and property."""
    span = LITERATURE.get(model, {}).get(key)
    if span is None:
        return None
    return span[0] <= simulated <= span[1]
