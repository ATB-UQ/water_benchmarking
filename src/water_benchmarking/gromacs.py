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

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

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
    #: Bulk transfer goes to the data movers (setonix-dm*), not the login node.
    #: Pawsey provides them for exactly this and asks that large transfers use them,
    #: so a multi-GB pull does not sit on a shared, throttled login node.  Measured
    #: from this host on 2026-08-29: the data mover is *not* faster per connection
    #: (~19 MB/s either way) -- the ceiling is per-stream, not the link -- which is
    #: why collect() also runs its pulls concurrently: 4 streams reach ~72 MB/s.
    "transfer_alias": "setonix-dm",
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
    """The same ladder as the GROMOS run, stage for stage -- with one difference.

    The two engines run the ten production nanoseconds differently, and the shared
    ladder is written for the GROMOS shape.  GROMOS runs them as ten *independent
    replicates* from the same equilibrated box, so each must draw fresh velocities
    or they would be ten copies of one trajectory; that is what
    imd.run_ladder sets generate_velocities for.  This driver instead *chains* its
    segments -- slurm_script grompps each one from the previous segment's .gro --
    so a chained segment must continue the velocities it inherits, not throw them
    away.  Regenerating them would restart the dynamics ten times from a fresh
    Maxwell distribution, at the same seed each time, which is neither a chain nor
    a set of replicates, and would leave a re-equilibration transient at the head of
    every segment that the GROMACS path (unlike the GROMOS one) does not discard.
    """
    ladder = []
    for s in imd.run_ladder():
        production = s.name.startswith("md_")
        ladder.append(
            MdpStage(
                s.name, s.steps, s.temperature, s.pressure_coupling,
                s.generate_velocities and not production,
                s.write_every, s.minimisation,
            )
        )
    return ladder


#: Water models this package ships its own .itp for, rather than taking one from
#: the stock force field.  The file is copied into the run directory and included
#: by name, so it travels to Setonix with the rest of the inputs.
PACKAGED_ITP = {"opc3.itp"}


def water_itp_source(name: str) -> Path:
    """The packaged .itp of that name, as installed beside this module."""
    return Path(__file__).parent / "data" / name


def write_topology(model: str, directory: Path) -> Path:
    """A .top holding nothing but water: no solute exists in this system.

    The model's [ moleculetype ] comes either from the stock gromos54a7.ff or, for
    a model that force field does not carry, from a copy of this package's own .itp
    placed beside the .top.  There is deliberately no default: an unknown model
    used to fall through to spc.itp, which produced a topology that grompp accepts,
    mdrun runs and every analysis reports -- as SPC, under the other model's name.
    """
    directory.mkdir(parents=True, exist_ok=True)
    water_itp = protocol.model(model).gromacs_itp
    if water_itp in PACKAGED_ITP:
        include = water_itp
        shutil.copyfile(water_itp_source(water_itp), directory / water_itp)
    else:
        include = f"{FORCEFIELD}/{water_itp}"
    path = directory / f"{model}.top"
    path.write_text(
        f'#include "{FORCEFIELD}/forcefield.itp"\n'
        f'#include "{include}"\n'
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
    # A packaged .itp is included by bare name, so it has to travel with the .top;
    # Setonix's GROMACS has no copy of it.
    packaged = protocol.model(model).gromacs_itp
    if packaged in PACKAGED_ITP:
        inputs.append(local_dir / packaged)
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


#: Written by trjconv into the TITLE; harmless, but a reminder the file is derived.
TRJCONV_GROUP = "System"


LOCAL_GMX = GROMACS_PREFIX / "bin" / "gmx"
LOCAL_GMXLIB = GROMACS_PREFIX / "share" / "gromacs" / "top"


#: Concurrent scp streams used to pull results back.  A single ssh connection
#: tops out around 19 MB/s to Pawsey from this host regardless of endpoint or
#: cipher -- the limit is per-stream, not the link -- and throughput scales very
#: nearly linearly: 2 streams measured 31 MB/s, 4 measured 72 MB/s.  Four keeps a
#: 2.2 GB model pull under a minute without opening a connection per file.
TRANSFER_STREAMS = 4


def collect(model: str, local_dir: Path) -> list[Path]:
    """Pull one model's results back from Setonix, over the data movers.

    The compressed .xtc comes back rather than a converted .g96: the same 1 ns
    segment is ~220 MB as .xtc and ~2.8 GB as text, and the conversion is done
    locally anyway.  Ten segments per model is then a 2.2 GB transfer instead of 28.

    Two things about how it is pulled, both measured rather than assumed:

    * **Not through the login node.**  `transfer_alias` is a data mover; Pawsey
      provides them for bulk transfer and asks that it not go through the shared
      login nodes.  It is not faster per connection -- that was measured -- it is
      simply the right endpoint, and it leaves the login node free.
    * **In parallel.**  One connection reaches ~19 MB/s whatever the endpoint, so
      the pull is split across TRANSFER_STREAMS connections, which is ~3.8x.
    """
    import subprocess
    from concurrent.futures import ThreadPoolExecutor

    local_dir.mkdir(parents=True, exist_ok=True)
    remote = f"{SETONIX['remote_root']}/{model}"
    host = SETONIX.get("transfer_alias", SETONIX["ssh_alias"])

    def scp(sources: Sequence[str]) -> subprocess.CompletedProcess:
        # One unquoted source argument per pattern: scp hands each to the remote
        # shell, which expands the glob.  Joining them into a single quoted string
        # makes scp look for one absurdly long literal filename instead.
        return subprocess.run(
            ["scp", "-q", "-o", "BatchMode=yes",
             *[f"{host}:{remote}/{p}" for p in sources], str(local_dir)],
            capture_output=True, text=True,
        )

    def pull(patterns: Sequence[str], required: bool) -> None:
        """Split the patterns over concurrent connections and run them together."""
        groups = [list(patterns[i::TRANSFER_STREAMS]) for i in range(TRANSFER_STREAMS)]
        with ThreadPoolExecutor(max_workers=TRANSFER_STREAMS) as pool:
            results = list(pool.map(scp, [g for g in groups if g]))
        failed = [r for r in results if r.returncode]
        if failed and required:
            raise RuntimeError(
                f"could not pull {list(patterns)} from {host}: "
                + "; ".join(r.stderr.strip() for r in failed)
            )

    # Per segment rather than one glob per file type, so the streams carry roughly
    # equal shares: the .xtc is ~220 MB and the .tpr and .edr are a rounding error,
    # so splitting by type would put the whole transfer on one connection.
    segments = [f"md_{n:02d}" for n in range(1, protocol.PRODUCTION_SEGMENTS + 1)]
    pull([f"{seg}.{ext}" for seg in segments for ext in ("xtc", "tpr", "edr")],
         required=True)
    # Nice to have, and absent often enough not to fail the transfer over.
    pull(["eq3.edr", "*.out"], required=False)

    return sorted(local_dir.glob("md_*.xtc"))


def to_g96(xtc: Path, tpr: Path, output: Path, remote_host: str | None = None) -> Path:
    """Convert a GROMACS trajectory to the g96 format trc.py reads.

    `-pbc mol` is not optional.  Without it GROMACS writes atoms wrapped
    individually, which splits roughly 80 of the 2048 molecules across the
    periodic boundary in every frame.  Those molecules then have a 5 nm "O-H bond",
    which silently ruins the box dipole (so the dielectric constant) and the
    molecular vectors (so the rotational correlation times).  assert_whole_molecules
    below is the check that this was actually done.
    """
    import subprocess

    # gmx_mpi is what the Setonix module provides; the local install is a serial
    # gmx at a fixed path and is not on PATH.
    binary = "gmx_mpi" if remote_host else str(LOCAL_GMX)
    command = (
        f"printf '{TRJCONV_GROUP}\n' | {binary} trjconv "
        f"-f {xtc} -s {tpr} -o {output} -pbc mol"
    )
    if remote_host:
        modules = " && ".join(f"module load {m}" for m in SETONIX["modules"])
        subprocess.run(
            ["ssh", "-o", "BatchMode=yes", remote_host, f"{modules} && {command}"],
            check=True, capture_output=True, text=True,
        )
    else:
        environment = dict(os.environ, GMXLIB=str(LOCAL_GMXLIB))
        subprocess.run(command, shell=True, check=True, capture_output=True,
                       text=True, env=environment)
    return output


def _legend_columns(xvg_text: str) -> dict:
    """Map each energy term to its column in the .xvg, by reading the legends.

    gmx energy writes the selected terms in its own internal order, not the order
    they were asked for, so pairing request order with column order silently swaps
    values -- Density and Potential came back as each other's numbers. The legends
    are the only reliable statement of what is in which column.
    """
    columns = {}
    for line in xvg_text.splitlines():
        if not line.startswith("@ s"):
            continue
        marker, _, remainder = line.partition(" legend ")
        if not remainder:
            continue
        index = int(marker.split()[1][1:])
        columns[remainder.strip().strip('"')] = index + 1   # column 0 is time
    return columns


def energy_series(edr_files: Sequence[Path], properties: Sequence[str]) -> dict:
    """Per-frame series for the named energy terms, via `gmx energy`.

    The GROMACS counterpart of ene_ana; density and potential energy come from
    here for both engines rather than from the trajectory.
    """
    import subprocess
    import tempfile

    series: dict[str, list] = {name: [] for name in properties}
    environment = dict(os.environ, GMXLIB=str(LOCAL_GMXLIB))
    for edr in edr_files:
        with tempfile.TemporaryDirectory() as work:
            output = Path(work) / "e.xvg"
            selection = "\n".join(properties) + "\n\n"
            subprocess.run(
                [str(LOCAL_GMX), "energy", "-f", str(edr), "-o", str(output)],
                input=selection, capture_output=True, text=True, check=True,
                env=environment,
            )
            text = output.read_text()
            data = np.loadtxt(output.open(), comments=("#", "@"))
            data = np.atleast_2d(data)
            for name, column in _legend_columns(text).items():
                if name in series:
                    series[name].append(data[:, column])

    missing = [name for name, values in series.items() if not values]
    if missing:
        raise KeyError(f"gmx energy returned no column for {missing}")
    return {name: np.concatenate(values) for name, values in series.items()}


def assert_whole_molecules(frame, tolerance: float = 0.02, r_oh: float = 0.1) -> None:
    """Fail loudly if a converted frame has molecules split across the boundary.

    Cheap, and the failure it catches is otherwise invisible: nothing downstream
    errors on a broken molecule, the numbers just come out wrong.

    ``r_oh`` is the model's constrained O-H distance -- pass it rather than relying
    on the SPC default.  OPC3's 0.0979 nm happens to sit inside the tolerance
    around 0.1, so the default would still pass, but only by 2% of the 20% margin
    this check is meant to have.
    """
    import numpy as np

    bond = np.linalg.norm(frame.positions[:, 1] - frame.positions[:, 0], axis=1)
    broken = int((np.abs(bond - r_oh) > tolerance).sum())
    if broken:
        raise AssertionError(
            f"{broken} molecules are split across the periodic boundary at "
            f"t = {frame.time} ps -- the trajectory was converted without `-pbc mol`"
        )
