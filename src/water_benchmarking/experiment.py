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
    # The single-molecule dipole correlation time is not the Debye time: tau_D is
    # collective and slower by roughly (2 eps + eps_inf) / (3 eps) ~ 1.5-2x.  The
    # reference is kept for context but excluded from the deviation table.
    "tau1_dipole": Reference(8.3, "ps", "Debye tau_D (collective); Ronne et al. 1997, JCP 107:5319"),
}

#: Properties whose experimental reference is a related but different quantity.
NOT_DIRECTLY_COMPARABLE = {"tau1_dipole"}

#: Footnotes for results that sit outside a model's published range for a reason
#: that is understood -- the range says "setup fault", the note says why not.
NOTES = {
    ("spc", "hov"): "published SPC dH_vap (41-44) is mostly at ~0.9 nm cutoffs; the 1.8 nm "
                    "cutoff here recovers more attractive energy, so a value at the top of "
                    "the range is the protocol, not an error",
    ("spc", "diffusion"): "the reported D carries the Yeh-Hummer finite-size correction "
                          "(+0.17); published SPC values are mostly uncorrected D_PBC, "
                          "for which this run gives 4.14",
    ("opc3", "dielectric"): "OPC3 was parameterised under Ewald/PME; this protocol is a "
                            "1.8 nm reaction field at eps_rf = 61, and eps is the property "
                            "most sensitive to the boundary condition (see results/"
                            "diagnostics.md) -- a shortfall here is the protocol, not a "
                            "setup fault",
    ("opc3", "diffusion"): "as for the dielectric constant: parameterised under PME, run "
                           "here under a reaction field, and D is the transport property "
                           "the cutoff moves most (the 0.9 nm control shifts it 7%)",
}

#: Published values for these models, as a check that the run reproduces the model.
LITERATURE = {
    "spc": {
        "density": (972.0, 985.0),
        "hov": (41.0, 44.0),
        "diffusion": (3.6e-9, 4.4e-9),   # 3.6 (Berendsen 1987) to 4.3 (van der Spoel 1998)
        "dielectric": (60.0, 70.0),      # 65 +/- 5 across RF and Ewald studies
        "tau2_HH": (0.9, 1.3),
    },
    "spce": {
        "density": (994.0, 1001.0),
        "hov": (46.0, 49.5),          # before the polarisation correction
        "diffusion": (2.4e-9, 2.8e-9),
        "dielectric": (68.0, 74.0),
        "tau2_HH": (1.7, 2.2),
    },
    # OPC3: Izadi & Onufriev 2016, J. Chem. Phys. 145:074501, Table III, at
    # 298.16 K and 1 bar under Ewald/PME:
    #     rho  = 0.996 +/- 0.001 g/cm^3      eps  = 78.4 +/- 1
    #     D    = 2.30 +/- 0.02 e-9 m^2/s     dHvap = 10.73 +/- 0.004 kcal/mol
    # dHvap converts to 44.89 kJ/mol, and it is a *self-polarisation corrected*
    # value: the paper's own SPC/E entry is 10.43 kcal/mol (43.64 kJ/mol), which is
    # the corrected number, not the ~49 kJ/mol a raw SPC/E run gives.  So the "hov"
    # range below is on the raw scale this benchmark reports -- 44.89 + OPC3's 7.03
    # correction ~ 51.9 -- and hov_polarisation_corrected is what compares with the
    # paper.  Ranges below are the paper's values with room for the reaction
    # field this protocol uses and OPC3 was not fitted under -- roughly the margin
    # the spce entry carries over its own published values.
    #
    # tau2_HH is absent because the paper does not report it: Table III lists ten
    # properties and no rotational correlation time.  It prints "[?]" in the
    # deviation table, which is the honest reading -- unchecked, not confirmed.
    "opc3": {
        "density": (992.0, 1000.0),      # 996 +/- 1
        "hov": (50.0, 53.5),             # raw; 44.89 + 7.03 polarisation ~ 51.9
        "diffusion": (2.1e-9, 2.6e-9),   # 2.30 +/- 0.02
        "dielectric": (73.0, 82.0),      # 78.4 +/- 1, the closest of any 3-point model
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
