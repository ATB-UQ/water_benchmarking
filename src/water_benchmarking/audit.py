"""Assemble a self-contained audit record for one run, so the published numbers
can be checked after the bulk trajectory data is gone.

A production trajectory is ~1 GB per nanosecond and this benchmark has 35 GB of
them; they are the one thing here that is both enormous and, once the analysis has
run, never read again.  What a reader actually needs in order to check a result is
much smaller: what was asked of the engine, what the engine said it did, the
energies the thermodynamic numbers come from, and a checksum of everything that was
removed so a restored copy can be proved identical.  That is what this writes.

The record follows `audit/spc_gromos/`, which was assembled by hand first:

    <out>/inputs/       .imd / .mdp, topology, starting configuration (gzipped)
    <out>/logs/         engine logs, trimmed; scheduler output
    <out>/provenance/   SHA256SUMS.txt over audited *and* excluded files, runs.csv
    <out>/results/      the analysed numbers for this run

**Energies are kept, trajectories are not.**  Density and heat of vaporisation are
recomputable from the .edr / .tre.gz at any time, because those are a few MB.  The
diffusion coefficient, the rotational correlation times and the dielectric constant
are not: they need the coordinates, and the record keeps their computed values and
the checksum of the trajectory they came from instead.  That asymmetry is deliberate
and is stated in the README each record carries.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

#: Files worth keeping verbatim: everything that says what was run, plus the
#: energies.  Ordered so the glob for a stage's inputs stays readable.
INPUT_PATTERNS = ("*.imd", "*.mdp", "*.top", "*.itp", "*.slurm")
#: Kept because they are small and are what density and dH_vap are computed from.
ENERGY_PATTERNS = ("*.edr", "*.tre.gz")
#: The bulk.  Excluded from the record, checksummed before removal.
TRAJECTORY_PATTERNS = ("*.xtc", "*.trc.gz", "*.trc")

#: A log is mostly periodic energy blocks.  Both engines mark them: md++ opens each
#: with a TIMESTEP block, GROMACS with a "Step Time" line.  The header before the
#: first one is what records the build, the host and every parameter the engine
#: actually parsed, and is kept whole.
BLOCK_MARKERS = (re.compile(r"^TIMESTEP\s*$"), re.compile(r"^ *Step +Time\s*$"))


def _blocks(lines: list[str]) -> list[int]:
    return [i for i, line in enumerate(lines)
            if any(m.match(line) for m in BLOCK_MARKERS)]


def trim_log(text: str, keep_tail: int = 40) -> str:
    """Header + first energy block + a marker + last block + the ending.

    A 1 ns md++ log is 46 MB, of which 99.9% is ten thousand energy blocks that say
    nothing the .tre does not.  Trimming rather than filtering keeps the log
    readable as a log: what the engine was asked to do survives in full, and the
    marker says exactly how much was removed rather than leaving a silent gap.
    """
    lines = text.splitlines(keepends=True)
    starts = _blocks(lines)
    if len(starts) < 3:
        return text

    block_len = starts[1] - starts[0]
    head = lines[: starts[0] + block_len]
    last = lines[starts[-1]: starts[-1] + block_len]
    tail = lines[max(starts[-1] + block_len, len(lines) - keep_tail):]
    removed = len(starts) - 2
    marker = (
        f"\n[audit] {removed} intermediate energy blocks removed "
        f"({block_len} lines each); the first and last are kept.  The full series "
        f"is in the .tre/.edr beside this record.\n\n"
    )
    return "".join(head) + marker + "".join(last) + "".join(tail)


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass
class Record:
    run_dir: Path
    out_dir: Path
    audited: list      # (relative path, sha256, bytes)
    excluded: list     # the same, for files not carried into the record


def _matching(run_dir: Path, patterns) -> list[Path]:
    found: list[Path] = []
    for pattern in patterns:
        found.extend(sorted(run_dir.glob(pattern)))
    return found


def build(run_dir: Path, out_dir: Path, results: dict | None = None) -> Record:
    """Write the record for one run directory and return what it covers."""
    run_dir, out_dir = Path(run_dir), Path(out_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)
    for sub in ("inputs", "logs", "provenance", "results"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    audited, excluded = [], []

    def record(path: Path, into: list) -> None:
        into.append((path.name, sha256(path), path.stat().st_size))

    # --- inputs: verbatim, plus the starting configuration compressed -----------
    for path in _matching(run_dir, INPUT_PATTERNS):
        shutil.copy2(path, out_dir / "inputs" / path.name)
        record(path, audited)
    for path in _matching(run_dir, ("water_*.cnf", "water_*.gro", "start.gro")):
        target = out_dir / "inputs" / (path.name + ".gz")
        with open(path, "rb") as src, gzip.open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)
        record(path, audited)

    # --- logs: trimmed, and the scheduler's own output --------------------------
    for path in sorted(run_dir.glob("*.log")):
        if path.name.endswith(".job.json"):
            continue
        (out_dir / "logs" / path.name).write_text(
            trim_log(path.read_text(errors="replace"))
        )
        record(path, audited)
    for path in sorted(run_dir.glob("*.out")):
        shutil.copy2(path, out_dir / "logs" / path.name)
        record(path, audited)

    # --- provenance: the scheduler's per-job records -----------------------------
    for path in sorted(run_dir.glob("*.job.json")):
        shutil.copy2(path, out_dir / "provenance" / path.name)
        record(path, audited)

    # --- energies stay in the run directory, but are named and checksummed -------
    for path in _matching(run_dir, ENERGY_PATTERNS):
        record(path, audited)
    # --- and the trajectories, which the record deliberately does not carry -------
    for path in _matching(run_dir, TRAJECTORY_PATTERNS):
        record(path, excluded)

    sums = out_dir / "provenance" / "SHA256SUMS.txt"
    sums.write_text(
        f"# Audit record for {run_dir}\n"
        "# 'audited' files are carried in this record or kept beside the run;\n"
        "# 'excluded' are the trajectories, removed after this was written.  A\n"
        "# restored copy can be proved identical with `sha256sum -c`.\n\n"
        + "".join(f"{h}  {n}\n" for n, h, _ in audited)
        + "\n# excluded (trajectory data, not retained)\n"
        + "".join(f"{h}  {n}\n" for n, h, _ in excluded)
    )

    with open(out_dir / "provenance" / "files.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "sha256", "bytes", "retained"])
        for name, digest, size in audited:
            writer.writerow([name, digest, size, "yes"])
        for name, digest, size in excluded:
            writer.writerow([name, digest, size, "no"])

    if results is not None:
        (out_dir / "results" / "aggregate.json").write_text(json.dumps(results, indent=1))

    return Record(run_dir=run_dir, out_dir=out_dir, audited=audited, excluded=excluded)
