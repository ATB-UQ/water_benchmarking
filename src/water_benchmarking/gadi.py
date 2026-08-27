"""Run the GROMOS ladder on Gadi through the gromos_job_wrapper ssh/PBS shim.

gadi_md.sh is a drop-in replacement for the md binary: it stages the inputs to
/scratch, renders a PBS script whose walltime it derives from NSTLIM, submits with
`qsub -W block=true`, and pulls the results back with size verification.  Nothing
here reimplements any of that -- this module only decides what to run, in what
order, and with how many ranks.

Segments are chained: each one starts from the previous one's final configuration,
so the ladder is strictly sequential within a model.  The two models are
independent and are meant to be run as two concurrent chains.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from . import box, forcefield, imd, protocol

#: Ranks to ask for.  24 is the latency choice: efficiency at 24 is only ~46%, but
#: a production segment is on the critical path and equilibration is not.
PRODUCTION_CORES = 24
EQUILIBRATION_CORES = 8

#: Seconds per step used to size the PBS walltime.  gadi_md.sh defaults to 0.08,
#: which is the 23k-atom peptide figure; this system is a quarter of that size, so
#: the default would ask for ~22 h and queue badly.  Refined from the eq3 log by
#: seconds_per_step_from_log() before production is launched.
INITIAL_SECONDS_PER_STEP = 0.02

COMPLETION_MARKER = "MD++ finished successfully"

#: Measured gadi wall time per ns against rank count, from
#: gromos_job_wrapper/deployment/BENCHMARKS.md.  Scaling is nowhere near linear --
#: 24 ranks is 1.67x faster than 8, not 3x (58% efficiency against 32%) -- so
#: converting an 8-rank measurement into a 24-rank walltime by the ratio of core
#: counts under-requests by about a third, and a segment killed just short of its
#: restart is the exact failure the peptide runs hit when GJW_GADI_SPS was 0.05.
GADI_HOURS_PER_NS = {1: 17.3, 4: 5.26, 8: 3.71, 16: 2.76, 24: 2.22}


def rank_scaling(from_cores: int, to_cores: int) -> float:
    """How much the per-step cost changes moving between two rank counts."""
    return GADI_HOURS_PER_NS[to_cores] / GADI_HOURS_PER_NS[from_cores]


@dataclass
class Segment:
    name: str
    model: str
    imd_file: Path
    conf_in: Path
    conf_out: Path
    trajectory: Path
    energies: Path
    log: Path
    cores: int


def prepare(model: str, run_dir: Path) -> Path:
    """Write topology, box and every .imd for one model into its run directory."""
    run_dir.mkdir(parents=True, exist_ok=True)
    forcefield.build_topology(model, run_dir)
    box.build(run_dir)
    imd.write_all(run_dir, model)
    return run_dir


def segments(model: str, run_dir: Path) -> list[Segment]:
    """The ordered ladder, with each stage reading the previous stage's output."""
    stages = imd.run_ladder()
    result = []
    previous = run_dir / f"water_{protocol.N_WATERS}.cnf"
    equilibrated = run_dir / f"{protocol.EQ_STAGES[-1][0]}.cnf"
    for stage in stages:
        production = stage.name.startswith("md_")
        segment = Segment(
            name=stage.name,
            model=model,
            imd_file=run_dir / f"{stage.name}.imd",
            # Equilibration is a ladder, production is a fan: every replicate
            # starts from the equilibrated box and they are all independent.
            conf_in=equilibrated if production else previous,
            conf_out=run_dir / f"{stage.name}.cnf",
            trajectory=run_dir / f"{stage.name}.trc",
            energies=run_dir / f"{stage.name}.tre",
            log=run_dir / f"{stage.name}.log",
            cores=PRODUCTION_CORES if production else EQUILIBRATION_CORES,
        )
        result.append(segment)
        previous = segment.conf_out
    return result


def seconds_per_step_from_log(log: Path, steps: int, scaling: float = 1.0) -> float | None:
    """Read md++'s own timing out of a finished log, to size later walltimes.

    Returns None if the log has no timing line, in which case the caller should
    keep the conservative default rather than guess.
    """
    if not log.exists():
        return None
    for line in log.read_text().splitlines():
        if "Wall time simulation" in line:
            seconds = float(line.split(":")[-1])
            return seconds / steps * scaling
    return None


def run_segment(
    segment: Segment,
    topology: Path,
    seconds_per_step: float = INITIAL_SECONDS_PER_STEP,
    dry_run: bool = False,
) -> None:
    """Submit one segment and block until it comes back."""
    if segment.conf_out.exists() and _finished(segment.log):
        return  # already done; the ladder is restartable

    environment = dict(os.environ)
    environment.update(
        {
            "GJW_GADI_NCORES": str(segment.cores),
            "GJW_GADI_SPS": f"{seconds_per_step:g}",
            # PBS truncates -N at 15 characters and gadi_md.sh keeps the tail, so
            # the model must be in the name or both models' jobs are indistinguishable
            # in the queue -- "w_spce_md_10" is 12 characters and survives intact.
            "GJW_JOB_NAME": f"w_{segment.model}_{segment.name}",
        }
    )

    command = [
        str(protocol.GADI_MD_SHIM),
        "@topo", str(topology.resolve()),
        "@conf", str(segment.conf_in.resolve()),
        "@input", str(segment.imd_file.resolve()),
        "@fin", str(segment.conf_out.resolve()),
        "@trc", str(segment.trajectory.resolve()),
        "@tre", str(segment.energies.resolve()),
    ]
    if dry_run:
        print(" ".join(f"{k}={v}" for k, v in sorted(environment.items()) if k.startswith("GJW_")))
        print(" ".join(command))
        return

    with open(segment.log, "w") as log_handle:
        subprocess.run(command, env=environment, stdout=log_handle,
                       stderr=subprocess.STDOUT, check=True)

    if not _finished(segment.log):
        raise RuntimeError(
            f"{segment.name}: md++ did not report success; see {segment.log}"
        )
    for output in (segment.trajectory, segment.energies):
        _compress(output)


def _finished(log: Path) -> bool:
    return log.exists() and COMPLETION_MARKER in log.read_text()


def _compress(path: Path) -> None:
    """Gzip a trajectory in place; 0.1 ps sampling makes these large."""
    if path.exists() and shutil.which("gzip"):
        subprocess.run(["gzip", "-f", str(path)], check=True)


def run_model(model: str, run_dir: Path, dry_run: bool = False) -> None:
    """Run the whole ladder for one model."""
    prepare(model, run_dir)
    topology = run_dir / f"{model}.top"
    ladder = segments(model, run_dir)

    seconds_per_step = INITIAL_SECONDS_PER_STEP
    for segment in ladder:
        run_segment(segment, topology, seconds_per_step, dry_run=dry_run)
        if segment.name == "eq3" and not dry_run:
            # eq3 is the first stage run at production settings; its measured rate
            # sizes every production walltime, at the production rank count.
            measured = seconds_per_step_from_log(
                segment.log, protocol.EQ_STAGES[-1][1],
                scaling=rank_scaling(EQUILIBRATION_CORES, PRODUCTION_CORES),
            )
            if measured:
                seconds_per_step = measured * 1.3   # headroom for queue-time variation


def submit_replicate(segment: Segment, topology: Path, seconds_per_step: float) -> Path:
    """Submit one production replicate without waiting for it.

    GJW_SUBMIT_ONLY leaves a job record beside the log; collect_replicate turns
    that back into local files once the job has run.  Blocking would serialise the
    ten replicates, which is the entire thing being avoided.
    """
    environment = dict(os.environ)
    environment.update({
        "GJW_GADI_NCORES": str(segment.cores),
        "GJW_GADI_SPS": f"{seconds_per_step:g}",
        "GJW_JOB_NAME": f"w_{segment.model}_{segment.name}",
        "GJW_SUBMIT_ONLY": "1",
    })
    command = [
        str(protocol.GADI_MD_SHIM),
        "@topo", str(topology.resolve()),
        "@conf", str(segment.conf_in.resolve()),
        "@input", str(segment.imd_file.resolve()),
        "@fin", str(segment.conf_out.resolve()),
        "@trc", str(segment.trajectory.resolve()),
        "@tre", str(segment.energies.resolve()),
    ]
    with open(segment.log, "w") as log_handle:
        subprocess.run(command, env=environment, stdout=log_handle,
                       stderr=subprocess.STDOUT, check=True)
    return Path(f"{segment.log}.job.json")


def collect_replicate(segment: Segment) -> bool:
    """Pull one finished replicate back; False if it is not ready yet."""
    record = Path(f"{segment.log}.job.json")
    if not record.exists():
        return False
    if segment.conf_out.exists() and _finished(segment.log):
        return True
    result = subprocess.run(
        [str(protocol.GADI_MD_SHIM), "--collect", str(record)],
        capture_output=True, text=True,
    )
    if result.returncode:
        return False
    for output in (segment.trajectory, segment.energies):
        _compress(output)
    return True


def submit_production(model: str, run_dir: Path, seconds_per_step: float | None = None) -> list:
    """Fan every production replicate out at once.

    Ten 1 ns replicates queued together finish in roughly the wall time of one,
    where chaining them takes ten times as long. Water equilibrates in
    picoseconds, so independent replicates from a common equilibrated box are as
    good as one continuous run for equilibrium averages -- and better for error
    bars, since they are genuinely independent.
    """
    topology = run_dir / f"{model}.top"
    if seconds_per_step is None:
        measured = seconds_per_step_from_log(
            run_dir / "eq3.log", protocol.EQ_STAGES[-1][1],
            scaling=rank_scaling(EQUILIBRATION_CORES, PRODUCTION_CORES),
        )
        seconds_per_step = (measured or INITIAL_SECONDS_PER_STEP) * 1.3

    records = []
    for segment in segments(model, run_dir):
        if not segment.name.startswith("md_"):
            continue
        records.append(submit_replicate(segment, topology, seconds_per_step))
    return records
