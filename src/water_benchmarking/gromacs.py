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
#: Values taken from gromacs_pipeline's 'setonix-gpu' host, which is the configuration
#: the peptide campaign actually runs -- note the GPU allocation is a separate
#: account (m72-gpu), not the m72 the CPU work bills to.
SETONIX = {
    "ssh_alias": "setonix",
    "remote_root": "/scratch/m72/mstroet/water_benchmarking",
    "account": "m72-gpu",
    "partition": "gpu",
    "gpus": 1,
    "cpus_per_task": 8,
    "modules": ("gromacs-amd-gfx90a/2026.1-mixed",),
    "binary": "gmx_mpi",
    # The whole ladder is ~10.1M steps; at the measured 0.0008 s/step for a system
    # three times this size that is well under three hours, so six is ample.
    "walltime": "06:00:00",
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
; Let GROMACS size the Verlet buffer (the default tolerance) rather than pinning
; rlist to the cutoff.  A zero buffer is not "closer to GROMOS": pairs that drift
; inside the cutoff between list updates are simply missed, which blew up SETTLE
; at step 38 of eq1.  The buffer only widens the neighbour LIST -- interactions
; are still cut at rvdw/rcoulomb = {protocol.CUTOFF} nm, so the energy is unchanged.
verlet-buffer-tolerance  = 0.005
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
    binary = SETONIX["binary"]
    threads = SETONIX["cpus_per_task"]
    # srun with one rank bound to the nearest GCD; this is the launcher shape the
    # peptide GROMACS runs use on this machine.
    launcher = (
        f"srun -l -u -N 1 -n 1 -c {threads} "
        f"--gpus-per-node={SETONIX['gpus']} --gpus-per-task=1 --gpu-bind=closest"
    )

    lines = [
        "#!/bin/bash",
        f"#SBATCH --account={SETONIX['account']}",
        f"#SBATCH --partition={SETONIX['partition']}",
        # Deliberately no --cpus-per-task: Setonix's GPU partition allocates 8 cores
        # per GPU implicitly and its submit filter rejects any explicit CPU request
        # ("cannot explicitly request CPU resources for GPU allocation").
        "#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node={SETONIX['gpus']}",
        f"#SBATCH --gpus-per-node={SETONIX['gpus']}",
        f"#SBATCH --time={SETONIX['walltime']}",
        f"#SBATCH --job-name=w_{model}",
        "#SBATCH --output=%x-%j.out",
        "",
        *(f"module load {m}" for m in SETONIX["modules"]),
        "",
        "# mdrun's -ntomp must agree with what SLURM allocated, or it oversubscribes.",
        f"export OMP_NUM_THREADS={threads}",
        f"cd {remote_dir}",
        "set -euo pipefail",
        "",
    ]

    previous = f"water_{protocol.N_WATERS}.gro"
    for stage in ladder:
        # -nb only.  Two flags that the peptide runs use are wrong here:
        #   -pme gpu is a fatal error under a reaction field (there is no PME);
        #   -bonded gpu is a fatal error in rigid water, which has no bonded
        #     interactions whatsoever -- SETTLE is a constraint, not a bonded term
        #     ("Bonded interactions can not be computed on a GPU").
        offload = "" if stage.minimisation else " -nb gpu"
        # -maxwarn 3 covers exactly three warnings, every one of them a deliberate
        # choice made to match GROMOS rather than an oversight:
        #   1. Berendsen thermostat is deprecated -- but it is what GROMOS weak
        #      coupling is, and swapping in V-rescale would break the comparison.
        #   2. Berendsen barostat likewise (NPT stages only).
        #   3. The GROMOS force fields were parametrised with a twin-range scheme;
        #      the peptide protocol deliberately runs single-range, as does this.
        lines += [
            f"{binary} grompp -f {stage.name}.mdp -c {previous} -p {model}.top "
            f"-o {stage.name}.tpr -maxwarn 3",
            f"{launcher} {binary} mdrun -deffnm {stage.name} "
            f"-ntomp {threads}{offload}",
            "",
        ]
        previous = f"{stage.name}.gro"
    lines.append('echo "LADDER COMPLETE"')
    return "\n".join(lines)


def submit(model: str, local_dir: Path, dry_run: bool = False) -> str:
    """Stage one model's inputs to Setonix and sbatch the ladder.

    Returns the SLURM job id.  The whole ladder is a single job: at ~13 minutes
    per nanosecond on one GCD the entire run fits inside one allocation, so there
    is nothing to gain from chaining fourteen dependent jobs.
    """
    import subprocess

    remote_dir = f"{SETONIX['remote_root']}/{model}"
    host = SETONIX["ssh_alias"]
    ssh_options = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=30"]

    script = local_dir / "run.slurm"
    script.write_text(slurm_script(model, remote_dir))

    inputs = sorted(local_dir.glob("*.mdp")) + [
        local_dir / f"{model}.top",
        local_dir / f"water_{protocol.N_WATERS}.gro",
        script,
    ]
    if dry_run:
        print(f"scp -> {host}:{remote_dir}: {[f.name for f in inputs]}")
        print(f"ssh {host} 'cd {remote_dir} && sbatch run.slurm'")
        return ""

    subprocess.run(["ssh", *ssh_options, host, f"mkdir -p {remote_dir}"], check=True)
    subprocess.run(
        ["scp", "-q", *ssh_options, *[str(f) for f in inputs], f"{host}:{remote_dir}/"],
        check=True,
    )
    result = subprocess.run(
        ["ssh", *ssh_options, host, f"cd {remote_dir} && sbatch run.slurm"],
        capture_output=True, text=True, check=True,
    )
    # "Submitted batch job 12345"
    return result.stdout.strip().split()[-1]
