"""The GROMACS mirror of the GROMOS runs, for Setonix GPU.

Same box, same models, same protocol -- the point is that any difference in the
measured properties is a difference between the engines, not between two
loosely-similar setups.  The settings below are the cross-engine equivalences
audited for the peptide campaign (/ssd1_nas_md/protein_validation/settings.md).

Two traps are handled explicitly and neither is obvious from the .mdp:

* `vdw-modifier` defaults to Potential-shift, which shifts the reported Lennard-Jones
  energy.  Forces are unaffected, so the trajectory is fine, but the potential
  energy is not -- and this benchmark computes the heat of vaporisation directly
  from the potential energy.  It must be `none` to match GROMOS.
* `-pme gpu` is a fatal error under a reaction field.  The GPU offload here is
  `-nb gpu -bonded gpu` and nothing else.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from . import imd, protocol

GROMACS_PREFIX = Path("/home/atb/opt/gromacs/2026.1")
FORCEFIELD = "gromos54a7.ff"

#: One GCD, not a whole node: at ~6k atoms a full GPU node runs at ~12% parallel
#: efficiency, so the two models are better run concurrently on separate GCDs.
SETONIX = {
    "account": "m72",
    "partition": "gpu",
    "gpus": 1,
    "cpus_per_task": 8,
    "modules": ("gromacs/2026.1",),
}


@dataclass
class MdpStage:
    name: str
    steps: int
    temperature: float
    pressure_coupling: bool
    generate_velocities: bool
    write_every: int
    minimisation: bool


def stages() -> list[MdpStage]:
    """The same ladder as the GROMOS run, stage for stage."""
    return [
        MdpStage(
            s.name, s.steps, s.temperature, s.pressure_coupling,
            s.generate_velocities, s.write_every, s.minimisation,
        )
        for s in imd.run_ladder()
    ]


def write_topology(model: str, directory: Path) -> Path:
    """A .top holding nothing but water: no solute exists in this system."""
    directory.mkdir(parents=True, exist_ok=True)
    water_itp = "spce.itp" if model == "spce" else "spc.itp"
    path = directory / f"{model}.top"
    path.write_text(
        f'#include "{FORCEFIELD}/forcefield.itp"\n'
        f'#include "{FORCEFIELD}/{water_itp}"\n'
        "\n"
        "[ system ]\n"
        f"{protocol.N_WATERS} {model.upper()} water, benchmark box\n"
        "\n"
        "[ molecules ]\n"
        f"SOL                 {protocol.N_WATERS}\n"
    )
    return path


def render_mdp(stage: MdpStage, model: str) -> str:
    """One .mdp, matched block for block to the corresponding .imd."""
    if stage.minimisation:
        integrator = f"""integrator               = steep
nsteps                   = {stage.steps}
emtol                    = 100.0
emstep                   = 0.01
"""
    else:
        integrator = f"""integrator               = md
dt                       = {protocol.TIMESTEP}
nsteps                   = {stage.steps}
"""

    if stage.minimisation:
        coupling = "tcoupl                   = no\npcoupl                   = no\n"
    else:
        pressure = (
            f"""pcoupl                   = berendsen
pcoupltype               = isotropic
tau_p                    = {protocol.TAU_P}
compressibility          = {protocol.COMPRESSIBILITY_BAR}
ref_p                    = {protocol.PRESSURE_BAR}
"""
            if stage.pressure_coupling
            else "pcoupl                   = no\n"
        )
        coupling = (
            f"""tcoupl                   = berendsen
tc-grps                  = System
tau_t                    = {protocol.TAU_T}
ref_t                    = {stage.temperature}
"""
            + pressure
        )

    velocities = (
        f"""gen_vel                  = yes
gen_temp                 = {stage.temperature}
gen_seed                 = {protocol.seed(model)}
continuation             = no
"""
        if stage.generate_velocities
        else "gen_vel                  = no\ncontinuation             = yes\n"
    )

    return f"""; {protocol.N_WATERS} {model.upper()} water -- {stage.name}
; mirror of the GROMOS protocol: dt {protocol.TIMESTEP} ps, cutoff {protocol.CUTOFF} nm,
; reaction field with epsilon_rf {protocol.EPSILON_RF:.0f}
{integrator}
; --- neighbour searching and interactions ---
cutoff-scheme            = Verlet
pbc                      = xyz
verlet-buffer-tolerance  = -1
nstlist                  = {protocol.PAIRLIST_UPDATE}
rlist                    = {protocol.CUTOFF}
coulombtype              = reaction-field
rcoulomb                 = {protocol.CUTOFF}
epsilon_r                = 1
epsilon_rf               = {protocol.EPSILON_RF:.0f}
vdwtype                  = cut-off
rvdw                     = {protocol.CUTOFF}
; Potential-shift (the default) would change the reported LJ energy and so the
; heat of vaporisation; GROMOS applies no shift.
vdw-modifier             = none
DispCorr                 = no

; --- constraints ---
; Water is rigid by SETTLE, which is what the GROMOS SOLVENTCONSTR SHAKE does.
constraints              = none
constraint-algorithm     = lincs

; --- coupling ---
{coupling}{velocities}
; --- centre of mass motion ---
comm-mode                = linear
comm-grps                = System
nstcomm                  = {protocol.NSCM}

; --- output ---
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = {stage.write_every}
compressed-x-grps        = System
nstenergy                = {stage.write_every}
nstlog                   = {stage.write_every}
"""


def write_all(directory: Path, model: str) -> list[Path]:
    """Write the topology and every .mdp for one model."""
    directory.mkdir(parents=True, exist_ok=True)
    write_topology(model, directory)
    written = []
    for stage in stages():
        path = directory / f"{stage.name}.mdp"
        path.write_text(render_mdp(stage, model))
        written.append(path)
    return written


def slurm_script(model: str, remote_dir: str) -> str:
    """A Setonix job that walks the whole ladder for one model on one GCD."""
    ladder = stages()
    lines = [
        "#!/bin/bash",
        f"#SBATCH --account={SETONIX['account']}",
        f"#SBATCH --partition={SETONIX['partition']}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={SETONIX['cpus_per_task']}",
        f"#SBATCH --gpus-per-node={SETONIX['gpus']}",
        "#SBATCH --time=24:00:00",
        f"#SBATCH --job-name=w_{model}",
        "#SBATCH --output=%x-%j.out",
        "",
        *(f"module load {m}" for m in SETONIX["modules"]),
        "",
        "# mdrun's internal thread count must agree with what SLURM allocated.",
        f"export OMP_NUM_THREADS={SETONIX['cpus_per_task']}",
        f"cd {remote_dir}",
        "set -e",
        "",
    ]

    previous = f"water_{protocol.N_WATERS}.gro"
    for stage in ladder:
        # -nb/-bonded only: -pme gpu is fatal under a reaction field.
        offload = "" if stage.minimisation else " -nb gpu -bonded gpu"
        lines += [
            f"gmx grompp -f {stage.name}.mdp -c {previous} -p {model}.top "
            f"-o {stage.name}.tpr -maxwarn 1",
            f"gmx mdrun -deffnm {stage.name} -ntmpi 1 -ntomp {SETONIX['cpus_per_task']}{offload}",
            "",
        ]
        previous = f"{stage.name}.gro"
    return "\n".join(lines)
