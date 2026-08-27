"""The one place the simulation protocol is defined.

Every constant here is either copied from the peptide validation protocol
(/ssd1_nas_md/protein_validation/runs/peptides_v3/...) or forced by the box size
choice.  The GROMOS .imd and the GROMACS .mdp are both rendered from these, which
is what makes the two engines comparable rather than merely similar.
"""
from pathlib import Path

# --- the two models under test -------------------------------------------------
# GROMOS solvent building-block names in 54A7.mtb.  The models differ only in
# charge: SPC -0.82/+0.41, SPC/E -0.8476/+0.4238.  Same geometry, same IAC.
MODELS = {"spc": "H2O", "spce": "H2OE"}

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
    """
    return SEED_BASE + 100 * sorted(MODELS).index(model) + replicate
