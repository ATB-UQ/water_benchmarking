"""The one place the simulation protocol is defined.

Every constant here is either copied from the peptide validation protocol
(/ssd1_nas_md/protein_validation/runs/peptides_v3/...) or forced by the box size
choice.  The GROMOS .imd and the GROMACS .mdp are both rendered from these, which
is what makes the two engines comparable rather than merely similar.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Model:
    """One water model: what each engine needs to build it, and its geometry.

    ``gromos_block`` is the solvent building-block name in 54A7.mtb, passed to
    make_top as @solv.  It is None for a model 54A7 does not carry -- such a model
    is GROMACS-only, and ``engines`` says so.  Nothing else about the benchmark
    changes: the protocol, the box and every analysis are engine- and
    model-independent.
    """

    name: str
    charges: tuple           # (qO, qH, qH), e
    r_oh: float              # nm, the constrained O-H distance
    r_hh: float              # nm, the constrained H-H distance
    seed_index: int          # fixed forever: seeds of published runs must not move
    gromos_block: str | None = None
    engines: tuple = ("gromos", "gromacs")
    #: GROMACS include file for the [ moleculetype ].  A bare name is resolved
    #: against FORCEFIELD (the stock gromos54a7.ff); "opc3.itp" is shipped by this
    #: package and copied into the run directory -- see gromacs.write_topology.
    gromacs_itp: str = "spc.itp"

    def __post_init__(self):
        # A model with no 54A7 building block cannot be run with GROMOS, and a
        # registry that claims otherwise would fail deep inside make_top with a
        # "Cannot find building block" that names neither the model nor the reason.
        if self.gromos_block is None and "gromos" in self.engines:
            raise ValueError(
                f"{self.name}: no gromos_block, so 'gromos' cannot be one of its engines"
            )


# --- the models under test -----------------------------------------------------
# SPC and SPC/E differ only in charge (-0.82/+0.41 vs -0.8476/+0.4238); same
# geometry, same IAC, both in 54A7.mtb.  OPC3 (Izadi & Onufriev, JCP 145:074501,
# 2016) is not in 54A7 and has its own geometry and LJ parameters, so it is run
# with GROMACS only -- the engines have already been shown to agree to <0.1% on
# density and dH_vap, which is what the two-engine comparison was for.
MODELS = {
    "spc": Model(
        name="spc",
        charges=(-0.82, 0.41, 0.41),
        r_oh=0.1,
        r_hh=0.163299,
        seed_index=0,
        gromos_block="H2O",
        gromacs_itp="spc.itp",
    ),
    "spce": Model(
        name="spce",
        charges=(-0.8476, 0.4238, 0.4238),
        r_oh=0.1,
        r_hh=0.163299,
        seed_index=1,
        gromos_block="H2OE",
        gromacs_itp="spce.itp",
    ),
    "opc3": Model(
        name="opc3",
        charges=(-0.895170, 0.447585, 0.447585),
        r_oh=0.0978882,
        r_hh=0.1598507,
        seed_index=2,
        gromos_block=None,
        engines=("gromacs",),
        gromacs_itp="opc3.itp",
    ),
}


def model(name: str) -> Model:
    """The registry entry for one model, or a ValueError naming the alternatives."""
    try:
        return MODELS[name]
    except KeyError:
        raise ValueError(
            f"unknown model {name!r}; expected one of {sorted(MODELS)}"
        ) from None


def engines_for(name: str) -> tuple:
    """The engines this model is run with."""
    return model(name).engines


def models_for_engine(engine: str) -> list:
    """Every model that can be run with this engine, sorted."""
    return sorted(n for n, m in MODELS.items() if engine in m.engines)


# --- system --------------------------------------------------------------------
# 2048 is the smallest power of two whose box edge clears twice the 1.8 nm cutoff:
# 1024 waters at bulk density give a ~= 3.1 nm < 2 * 1.8, which breaks minimum image.
N_WATERS = 2048
ATOMS_PER_WATER = 3
N_ATOMS = N_WATERS * ATOMS_PER_WATER

TEMPERATURE = 298.15          # K
PRESSURE_GROMOS = 0.06102     # kJ mol^-1 nm^-3 == 1 atm
PRESSURE_BAR = 1.01325        # bar, the GROMACS spelling of the same thing

# --- interactions (peptide protocol, single-range) ------------------------------
CUTOFF = 1.8                  # nm; RCUTP == RCUTL == RCRF == ASHAPE
EPSILON_RF = 61.0             # the peptide-protocol value, used for BOTH models.
                              # SPC/E is often published with ~71; keeping 61 makes
                              # this a benchmark of the protocol ATB actually runs.
PAIRLIST_UPDATE = 5
PAIRLIST_GRID_SIZE = 0.4

# --- integration ---------------------------------------------------------------
TIMESTEP = 0.001              # ps
TAU_T = 0.1                   # ps
TAU_P = 0.5                   # ps
COMPRESSIBILITY_GROMOS = 4.575e-4   # kJ^-1 mol nm^3
COMPRESSIBILITY_BAR = 2.755e-5      # bar^-1, the same number in GROMACS units
NSCM = 1000                   # COM motion removal interval

# --- the run ladder ------------------------------------------------------------
# 0.1 ps sampling: tau_2 of water is ~2 ps, so the rotational correlation function
# needs frames far finer than the 50 ps the peptide runs wrote.
PRODUCTION_SEGMENTS = 10
PRODUCTION_STEPS = 1_000_000          # 1 ns per segment
PRODUCTION_WRITE = 100                # 0.1 ps
EQUILIBRATION_WRITE = 1000

EMIN_STEPS = 1000
EQ_STAGES = (
    # name,  steps,   temperature, pressure coupling, generate velocities
    ("eq1",  10_000,  50.0,        False, True),
    ("eq2",  20_000,  TEMPERATURE, False, False),
    ("eq3", 100_000,  TEMPERATURE, True,  False),
)

SEED_BASE = 770_000           # + model index, mirroring PEP_SEED_BASE

# --- paths ---------------------------------------------------------------------
GJW_LIB = Path(
    "/home/atb/ATB/gromos_job_wrapper/src/gromos_job_wrapper/lib"
)
SOURCE_WATER_BOX = GJW_LIB / "H2O_box.g96"      # 5384 SPC, cubic 5.4937 nm
FORCEFIELD_MTB = GJW_LIB / "54A7.mtb"
FORCEFIELD_IFP = GJW_LIB / "54A7.ifp"
# There is no single right ene_ana library.  Two md++ builds of the same
# nominal version write .tre files with different block layouts under the same
# ENEVERSION stamp (2023-04-15): the local /opt/gromos/1.6.0 md parses only with
# the library shipped beside it, and the gadi md_mpi build only with the copy in
# gromos_job_wrapper.  Each fails on the other with "Tried to read an integer for
# NUM_..." -- the stamp does not identify the layout, so the library is chosen
# per file by trying each until one parses.
ENE_ANA_LIBS = (
    Path("/opt/gromos/1.6.0/share/md++/ene_ana.md++.lib"),
    GJW_LIB / "ene_ana.md++.lib",
    GJW_LIB / "ene_ana_2015-06-23-A.md++.lib",
)
ENE_ANA_LIB = ENE_ANA_LIBS[0]   # kept for callers that only want a default

GROMOS_BIN = Path("/opt/gromos/1.6.0/bin")
GADI_MD_SHIM = Path("/home/atb/ATB/gromos_job_wrapper/deployment/gadi_md.sh")

RUN_ROOT = Path("/ssd1_nas_md/water_benchmarking")


def seed(model: str, replicate: int = 0) -> int:
    """Velocity seed for one model, optionally for one production replicate.

    Production is run as ten independent 1 ns replicates rather than one chained
    10 ns trajectory: they queue concurrently, so the wall time is one segment
    rather than ten.  That only works if each draws its own velocities -- started
    from the same equilibrated box with the same seed they would be ten copies of
    the same trajectory.

    The index is a fixed field of the model rather than its position in a sorted
    MODELS, so that adding a model cannot renumber the seeds of runs already
    published: sorted() would have put "opc3" first and moved both of them.
    """
    return SEED_BASE + 100 * MODELS[model].seed_index + replicate
