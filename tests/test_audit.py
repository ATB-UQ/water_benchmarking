"""The audit record is what survives the trajectories, so it has to be right.

Once the .xtc and .trc.gz are deleted these files are the only evidence that the
published numbers came from the runs they claim to.  A record that silently omitted
a file, or checksummed it wrongly, would not be discovered until someone tried to
check a result and could not.
"""
import gzip

from water_benchmarking import audit


def _md_log(blocks: int) -> str:
    """An md++-shaped log: header, then `blocks` periodic energy blocks."""
    head = "MD++\n========\nversion    :     1.6.0\nPARAMETERS\n...\nEND\n"
    body = "".join(
        f"TIMESTEP\n{i} {i * 0.5}\nEND\nENERGIES\nE_Total : {-76786.0 + i}\nEND\n"
        for i in range(blocks)
    )
    return head + body + "MD++ finished successfully\n"


def test_trim_keeps_the_header_and_both_end_blocks():
    text = _md_log(500)
    trimmed = audit.trim_log(text)

    assert "version    :     1.6.0" in trimmed        # the header survives whole
    assert "PARAMETERS" in trimmed
    assert "E_Total : -76786.0" in trimmed            # first block
    assert f"E_Total : {-76786.0 + 499}" in trimmed   # last block
    assert "MD++ finished successfully" in trimmed    # and how it ended
    assert "498 intermediate energy blocks removed" in trimmed
    assert len(trimmed) < len(text) / 10


def test_trim_leaves_a_short_log_alone():
    """Nothing to gain, and a truncated marker in a 2-block log would mislead."""
    text = _md_log(2)
    assert audit.trim_log(text) == text


def test_record_carries_inputs_and_checksums_what_it_drops(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    (run / "md_01.mdp").write_text("integrator = md\n")
    (run / "spce.top").write_text("#include \"x.itp\"\n")
    (run / "water_2048.gro").write_text("box\n")
    (run / "md_01.log").write_text(_md_log(50))
    (run / "md_01.edr").write_bytes(b"energies")
    (run / "md_01.xtc").write_bytes(b"a trajectory, notionally enormous")

    record = audit.build(run, tmp_path / "record")
    out = tmp_path / "record"

    assert (out / "inputs" / "md_01.mdp").exists()
    assert (out / "inputs" / "spce.top").exists()
    # The starting box is carried compressed rather than dropped.
    with gzip.open(out / "inputs" / "water_2048.gro.gz", "rb") as handle:
        assert handle.read() == b"box\n"
    # The log is carried, trimmed.
    assert "blocks removed" in (out / "logs" / "md_01.log").read_text()

    names = {n for n, _, _ in record.audited}
    assert {"md_01.mdp", "spce.top", "water_2048.gro", "md_01.log", "md_01.edr"} <= names
    # The energy file is audited (it stays beside the run); the trajectory is not.
    assert [n for n, _, _ in record.excluded] == ["md_01.xtc"]

    sums = (out / "provenance" / "SHA256SUMS.txt").read_text()
    assert "md_01.xtc" in sums, "a dropped file must still be checksummed"
    assert "excluded (trajectory data, not retained)" in sums
    assert audit.sha256(run / "md_01.xtc") in sums


def test_every_file_appears_exactly_once_in_the_manifest(tmp_path):
    """A file counted twice, or not at all, makes the manifest useless as a check."""
    run = tmp_path / "run"
    run.mkdir()
    for name in ("a.mdp", "b.top", "c.edr", "d.xtc", "e.slurm"):
        (run / name).write_bytes(b"x")

    record = audit.build(run, tmp_path / "record")
    listed = [n for n, _, _ in record.audited] + [n for n, _, _ in record.excluded]
    assert sorted(listed) == ["a.mdp", "b.top", "c.edr", "d.xtc", "e.slurm"]
    assert len(listed) == len(set(listed))
