"""The trajectory reader is shared by both engines, so it has to be exact."""
import numpy as np
import pytest

from water_benchmarking import box, trc


def _write_trajectory(path, frames, edge=4.0):
    text = ["TITLE", "test", "END"]
    for index, positions in enumerate(frames):
        text += ["TIMESTEP", f"{index:15d}{index * 0.1:20.9f}", "END", "POSITIONRED"]
        for molecule in positions:
            for atom in molecule:
                text.append("".join(f"{c:15.9f}" for c in atom))
        text += ["END", "GENBOX", "    1",
                 "".join(f"{edge:15.9f}" for _ in range(3)),
                 "".join(f"{90.0:15.9f}" for _ in range(3)),
                 "".join(f"{0.0:15.9f}" for _ in range(3)),
                 "".join(f"{0.0:15.9f}" for _ in range(3)), "END"]
    path.write_text("\n".join(text) + "\n")


def test_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    frames = rng.normal(size=(3, 5, 3, 3))
    path = tmp_path / "t.trc"
    _write_trajectory(path, frames)

    read = list(trc.read_frames(path, n_molecules=5))
    assert len(read) == 3
    assert np.allclose(read[1].positions, frames[1], atol=1e-8)
    assert read[2].time == pytest.approx(0.2)
    assert read[0].edge == pytest.approx(4.0)


def test_read_all_drops_repeated_boundary_frames(tmp_path):
    """A segment restarts from the previous final frame; the duplicate must go.

    Keeping it would add a zero-displacement step at every segment boundary and
    bias the diffusion coefficient low.
    """
    rng = np.random.default_rng(1)
    first = rng.normal(size=(3, 2, 3, 3))
    second = rng.normal(size=(3, 2, 3, 3))
    a, b = tmp_path / "a.trc", tmp_path / "b.trc"
    _write_trajectory(a, first)
    _write_trajectory(b, second)          # restarts at t = 0.0, overlapping a

    frames = list(trc.read_all([a, b], n_molecules=2))
    times = [f.time for f in frames]
    assert times == sorted(times)
    assert len(times) == len(set(times))
