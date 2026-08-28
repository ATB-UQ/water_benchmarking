"""grompp is the only thing that can say the topology is really well formed.

Everything in test_models.py checks the files this package writes against what it
meant to write.  That cannot catch the class of mistake that matters most here:
the packaged OPC3 .itp restates amber19sb's sigma/epsilon parameters as the C6/C12
gromos54a7.ff's [ defaults ] require, and a combination rule that disagrees with
the numbers is not a syntax error -- it is a topology that builds, runs, and is
wrong.  So this runs the real preprocessor and reads the parameters back out of
the binary .tpr it produces.

Skipped when the local GROMACS is not installed; it is not on Setonix that this
needs to pass, it is here, before ten nanoseconds are spent.
"""
import os
import shutil
import subprocess

import pytest

from water_benchmarking import box, gromacs, protocol

pytestmark = pytest.mark.skipif(
    not gromacs.LOCAL_GMX.exists(), reason=f"no local GROMACS at {gromacs.LOCAL_GMX}"
)

#: The three warnings the run script also passes -maxwarn for, every one of them a
#: deliberate choice: the two deprecated Berendsen couplings (which are what GROMOS
#: weak coupling is) and the GROMOS twin-range parametrisation notice, which this
#: protocol answers by running single-range on purpose.
MAXWARN = "3"


@pytest.fixture(scope="module")
def water_gro(tmp_path_factory):
    directory = tmp_path_factory.mktemp("box")
    box.build(directory)
    return directory / f"water_{protocol.N_WATERS}.gro"


def _grompp(model, stage, directory, water_gro):
    gromacs.write_all(directory, model)
    shutil.copy(water_gro, directory)
    tpr = directory / f"{stage}.tpr"
    result = subprocess.run(
        [str(gromacs.LOCAL_GMX), "grompp",
         "-f", f"{stage}.mdp", "-p", f"{model}.top", "-c", water_gro.name,
         "-o", tpr.name, "-maxwarn", MAXWARN],
        cwd=directory, capture_output=True, text=True,
        env={**os.environ, "GMXLIB": str(gromacs.LOCAL_GMXLIB)},
    )
    assert result.returncode == 0, result.stderr[-3000:]
    return tpr


def _dump(tpr):
    result = subprocess.run(
        [str(gromacs.LOCAL_GMX), "dump", "-s", tpr.name],
        cwd=tpr.parent, capture_output=True, text=True,
        env={**os.environ, "GMXLIB": str(gromacs.LOCAL_GMXLIB)},
    )
    assert result.returncode == 0, result.stderr[-3000:]
    return result.stdout


@pytest.mark.parametrize("model", sorted(protocol.MODELS))
@pytest.mark.parametrize("stage", ["emin", "md_01"])
def test_every_model_preprocesses(model, stage, tmp_path, water_gro):
    """A model whose .top grompp rejects cannot be run; find that out here."""
    _grompp(model, stage, tmp_path / f"{model}_{stage}", water_gro)


def test_opc3_tpr_carries_the_published_parameters(tmp_path, water_gro):
    """Read the model back out of the .tpr: charges, geometry and the LJ pair.

    The C6/C12 assertion is the point of this file.  Had the .itp been written in
    amber's sigma/epsilon convention and included under gromos54a7.ff's comb-rule 1
    -- which is what happens if you copy the stock file -- grompp would have taken
    0.3174 as C6 and 0.6837 as C12 and reported no error at all.
    """
    dump = _dump(_grompp("opc3", "md_01", tmp_path / "opc3", water_gro))
    model = protocol.MODELS["opc3"]

    oxygen, hydrogen = None, None
    for line in dump.splitlines():
        if "atom[     0]" in line:
            oxygen = line
        elif "atom[     1]" in line:
            hydrogen = line
    assert f"q={model.charges[0]: .5e}".replace(" ", "") in oxygen.replace(" ", "")
    assert f"q={model.charges[1]: .5e}".replace(" ", "") in hydrogen.replace(" ", "")

    sigma, eps = 0.317427035094, 0.683690704000
    lj = [line for line in dump.splitlines() if "LJ_SR" in line]
    c6 = float(lj[0].split("c6=")[1].split(",")[0])
    c12 = float(lj[0].split("c12=")[1].strip())
    assert c6 == pytest.approx(4 * eps * sigma ** 6, rel=1e-6)
    assert c12 == pytest.approx(4 * eps * sigma ** 12, rel=1e-6)
    # Hydrogen has no Lennard-Jones term, so every other pair must be zero.
    for line in lj[1:]:
        assert "c6= 0.00000000e+00" in line and "c12= 0.00000000e+00" in line

    settle = [line for line in dump.splitlines() if "SETTLE" in line and "doh=" in line]
    doh = float(settle[0].split("doh=")[1].split(",")[0])
    dhh = float(settle[0].split("dhh=")[1].strip())
    assert doh == pytest.approx(model.r_oh, rel=1e-6)
    assert dhh == pytest.approx(model.r_hh, rel=1e-6)
