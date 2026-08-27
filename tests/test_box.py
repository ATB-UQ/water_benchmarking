"""The box cutter must produce a periodic box md++ will accept."""
from pathlib import Path

import numpy as np
import pytest

from water_benchmarking import box, protocol


@pytest.fixture(scope="module")
def water_box():
    return box.cut(protocol.SOURCE_WATER_BOX, protocol.N_WATERS)


def test_exact_molecule_count(water_box):
    assert water_box.n_molecules == protocol.N_WATERS


def test_box_clears_the_cutoff(water_box):
    # A half-box below the cutoff would break the minimum image convention: a
    # molecule would interact with two images of the same partner.
    assert water_box.edge / 2 > protocol.CUTOFF


def test_molecules_are_rigid(water_box):
    """Wrapping must move whole molecules; SHAKE fails at step 0 otherwise."""
    oh1 = np.linalg.norm(water_box.positions[:, 1] - water_box.positions[:, 0], axis=1)
    oh2 = np.linalg.norm(water_box.positions[:, 2] - water_box.positions[:, 0], axis=1)
    hh = np.linalg.norm(water_box.positions[:, 2] - water_box.positions[:, 1], axis=1)
    assert np.allclose(oh1, 0.1, atol=1e-4)
    assert np.allclose(oh2, 0.1, atol=1e-4)
    assert np.allclose(hh, 0.163299, atol=1e-4)


def test_no_clashes(water_box):
    assert box.min_image_oxygen_distance(water_box) > 0.22


def test_density_is_liquid_water(water_box):
    assert 950.0 < water_box.density < 1010.0


def test_written_cnf_has_the_blocks_md_needs(tmp_path, water_box):
    box.write_cnf(water_box, tmp_path / "b.cnf", "test")
    text = (tmp_path / "b.cnf").read_text()
    for block in ("POSITION", "LATTICESHIFTS", "GENBOX"):
        assert block in text, f"{block} missing; md++ refuses the configuration"
