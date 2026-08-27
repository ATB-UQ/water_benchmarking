"""Render GROMOS .imd input files for a pure-water system.

The blocks below are the peptide validation protocol with the solute-only parts
removed: there is no solute, so no POSITIONRES, no ROTTRANS, and MULTIBATH couples
one bath instead of two.  Everything else -- cutoffs, reaction field, constraints,
thermostat and barostat -- is copied unchanged, which is the entire point: the
benchmark has to measure the protocol the peptide runs actually use.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from . import protocol


@dataclass
class Stage:
    """One .imd file: a named piece of the run ladder."""

    name: str
    steps: int
    temperature: float
    pressure_coupling: bool
    generate_velocities: bool
    write_every: int
    minimisation: bool = False
    initial_time: float = 0.0


def run_ladder() -> list[Stage]:
    """emin -> NVT at 50 K -> NVT at 298 K -> NPT -> 10 x 1 ns NPT production."""
    stages = [
        Stage("emin", protocol.EMIN_STEPS, protocol.TEMPERATURE, False, False,
              protocol.EQUILIBRATION_WRITE, minimisation=True)
    ]
    for name, steps, temperature, barostat, gen_vel in protocol.EQ_STAGES:
        stages.append(Stage(name, steps, temperature, barostat, gen_vel,
                            protocol.EQUILIBRATION_WRITE))
    for segment in range(1, protocol.PRODUCTION_SEGMENTS + 1):
        stages.append(
            Stage(
                f"md_{segment:02d}",
                protocol.PRODUCTION_STEPS,
                protocol.TEMPERATURE,
                True,
                False,
                protocol.PRODUCTION_WRITE,
                initial_time=(segment - 1) * protocol.PRODUCTION_STEPS * protocol.TIMESTEP,
            )
        )
    return stages


def _energy_minimisation() -> str:
    return """ENERGYMIN
#    NTEM    NCYC    DELE    DX0     DXM     NMIN    FLIM
        1       1    0.01    0.01    0.05     100     0.0
END
"""


def _multibath(temperature: float) -> str:
    # One bath: the system is a single molecular species, so there is no solute /
    # solvent split to make.  DOFSET 1 covers every atom.
    return f"""MULTIBATH
# ALGORITHM: 0 = weak coupling
#  ALGORITHM
          0
#  NBATHS
   1
# TEMP0(1 ... NBATHS)  TAU(1 ... NBATHS)
   {temperature:12.5f} {protocol.TAU_T:8.2f}
#  DOFSET
   1
# LAST(1 ... DOFSET)  COMBATH(1 ... DOFSET)  IRBATH(1 ... DOFSET)
   {protocol.N_ATOMS:8d} {1:8d} {1:8d}
END
"""


def _pressure_scale(active: bool) -> str:
    couple = 2 if active else 0
    return f"""PRESSURESCALE
#    COUPLE:    off(0), calc(1), scale(2)
#    SCALE:  off(0), iso(1), aniso(2), full(3), semianiso(4)
#    VIRIAL: none(0), atomic(1), group(2)
# COUPLE   SCALE    COMP    TAUP  VIRIAL
   {couple:6d} {1:7d} {protocol.COMPRESSIBILITY_GROMOS:9.7f} {protocol.TAU_P:7.2f} {2:7d}
# SEMIANISOTROPIC COUPLINGS(X, Y, Z)
       1     1     1
# PRES0(1...3,1...3)  1 atm = 0.06102 kJ mol^-1 nm^-3
 {protocol.PRESSURE_GROMOS:7.5f}       0       0
       0 {protocol.PRESSURE_GROMOS:7.5f}       0
       0       0 {protocol.PRESSURE_GROMOS:7.5f}
END
"""


def render(stage: Stage, model: str, title: str) -> str:
    """Build the complete .imd text for one stage."""
    # The pairlist is rebuilt every step during minimisation (the standard, non-grid
    # algorithm); grid pairlists assume the configuration changes only a little.
    pairlist_algorithm = 0 if stage.minimisation else 1
    generate = 1 if stage.generate_velocities else 0
    initial_temperature = stage.temperature if stage.generate_velocities else 0.0

    blocks = [
        f"TITLE\n{title}\nEND\n",
        f"""SYSTEM
#      NPM      NSM
         0      {protocol.N_WATERS}
END
""",
        f"""STEP
#   NSTLIM         T        DT
{stage.steps:10d} {stage.initial_time:9.1f} {protocol.TIMESTEP:9.3f}
END
""",
        f"""INITIALISE
#    NTIVEL   NTISHK   NTINHT    NTINHB    NTISHI  NTIRTC     NTICOM   NTISTI      IG     TEMPI
      {generate}       0        0         0         0       0          0        0       {protocol.seed(model)}     {initial_temperature:.2f}
END
""",
        f"""BOUNDCOND
#      NTB    NDFMIN
       1      {6 if not stage.minimisation else 0}
END
""",
        """FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
  1     1     1     1     1     1
# NEGR    NRE(1)
     1     %d
END
""" % protocol.N_ATOMS,
        # NTC 1, not the peptide protocol's NTC 2: md++ refuses solute SHAKE in a
        # solvent-only system ("solvent only simulation does not work with SHAKE
        # for solute").  Nothing is lost -- there is no solute to constrain, and
        # the rigid water geometry comes from SOLVENTCONSTR via NTCS below.
        """CONSTRAINT
#      NTC       NTCP   NTCP0(1)     NTCS      NTCS0(1)
         1          1    0.00001        1      0.00001
END
""",
        f"""PAIRLIST
#    algorithm    NSNB    RCUTP    RCUTL    SIZE    TYPE
     {pairlist_algorithm}  {protocol.PAIRLIST_UPDATE}  {protocol.CUTOFF}  {protocol.CUTOFF}     {protocol.PAIRLIST_GRID_SIZE}     0
END
""",
        f"""NONBONDED
# RCRF and ASHAPE track the long-range cutoff: the reaction field is the continuum
# beyond the cutoff, so a radius disagreeing with RCUTL is an inconsistent
# electrostatic treatment.
# NLRELE    APPAK      RCRF     EPSRF  NSLFEXCL
       1      0.0       {protocol.CUTOFF}          {protocol.EPSILON_RF:.0f}     1
# NSHAPE   ASHAPE    NA2CLC   TOLA2   EPSLS
       -1       {protocol.CUTOFF}        2   1e-10       0
# NKX    NKY   NKZ    KCUT
   10     10    10     100
# NGX   NGY   NGZ  NASORD  NFDORD   NALIAS  NSPORD
   32    32    32       3       2        3       4
# NQEVAL   FACCUR   NRDGRD   NWRGRD   NLRLJ    SLVDNS
  100000      1.6        0        0       0      33.3
END
""",
        f"""PRINTOUT
#     NTPR          NTPP
        {stage.write_every:d}          0
END
""",
        """COVALENTFORM
#    NTBBH    NTBAH     NTBDN
         0         0         0
END
""",
    ]

    if stage.minimisation:
        blocks.append(_energy_minimisation())
    else:
        blocks.append(f"""COMTRANSROT
#   NSCM
    {protocol.NSCM}
END
""")
        blocks.append(_multibath(stage.temperature))
        blocks.append(f"""WRITETRAJ
#    NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
      {stage.write_every:d}         0         0         0      {stage.write_every:d}         0         0
END
""")
        blocks.append(_pressure_scale(stage.pressure_coupling))

    return "".join(blocks)


def write_all(directory: Path, model: str) -> list[Path]:
    """Write every .imd of the ladder for one model."""
    directory.mkdir(parents=True, exist_ok=True)
    written = []
    for stage in run_ladder():
        title = (
            f"{protocol.N_WATERS} {model.upper()} water, {stage.name}\n"
            f"benchmark protocol: dt {protocol.TIMESTEP} ps, cutoff "
            f"{protocol.CUTOFF} nm, EPSRF {protocol.EPSILON_RF:.0f}"
        )
        path = directory / f"{stage.name}.imd"
        path.write_text(render(stage, model, title))
        written.append(path)
    return written
