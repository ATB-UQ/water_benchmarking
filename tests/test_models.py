"""The model registry, and the topology each engine builds from it.

None of this is arithmetic anyone can eyeball either: a water model that is wrong
in the topology still runs, still converges and still produces a full set of
plausible numbers under the wrong name.  These tests exist because that failure is
silent -- write_topology used to fall through to spc.itp for any model it did not
recognise, so `--model opc3` would have benchmarked SPC and labelled it OPC3.
"""
from pathlib import Path

import numpy as np
import pytest

from water_benchmarking import forcefield, gromacs, protocol


# --- the registry --------------------------------------------------------------

def test_every_model_is_neutral():
    for name, model in protocol.MODELS.items():
        assert len(model.charges) == 3, name
        assert abs(sum(model.charges)) < 1e-9, name
        # A water model's two hydrogens are equivalent by construction.
        assert model.charges[1] == model.charges[2], name


def test_geometry_is_consistent_with_the_hoh_angle():
    """r_HH must be the chord r_OH subtends at the model's angle."""
    for name, model in protocol.MODELS.items():
        angle = 2 * np.arcsin(model.r_hh / (2 * model.r_oh))
        assert 100.0 < np.degrees(angle) < 115.0, f"{name}: {np.degrees(angle):.2f} deg"


def test_charges_are_the_published_ones():
    assert protocol.MODELS["spc"].charges == (-0.82, 0.41, 0.41)
    assert protocol.MODELS["spce"].charges == (-0.8476, 0.4238, 0.4238)
    # Izadi & Onufriev 2016, J. Chem. Phys. 145:074501.
    assert protocol.MODELS["opc3"].charges == (-0.895170, 0.447585, 0.447585)


def test_expected_charges_covers_every_model():
    """The GROMACS dielectric path reads this dict, not the registry directly."""
    assert set(forcefield.EXPECTED_CHARGES) == set(protocol.MODELS)
    for name, charges in forcefield.EXPECTED_CHARGES.items():
        assert charges == protocol.MODELS[name].charges


def test_adding_a_model_does_not_move_the_published_seeds():
    """spc and spce were run at these seeds; a third model must not renumber them.

    seed() indexed sorted(MODELS) once, so adding "opc3" -- which sorts first --
    would have shifted both of the published models by one.
    """
    assert protocol.seed("spc") == 770_000
    assert protocol.seed("spce") == 770_100
    assert protocol.seed("opc3") == 770_200
    assert protocol.seed("spce", replicate=7) == 770_107


def test_seeds_are_unique_across_models_and_replicates():
    seeds = [protocol.seed(m, r)
             for m in protocol.MODELS
             for r in range(protocol.PRODUCTION_SEGMENTS)]
    assert len(set(seeds)) == len(seeds)


def test_engine_membership():
    assert protocol.engines_for("spc") == ("gromos", "gromacs")
    assert protocol.engines_for("opc3") == ("gromacs",)
    assert protocol.models_for_engine("gromos") == ["spc", "spce"]
    assert protocol.models_for_engine("gromacs") == ["opc3", "spc", "spce"]


def test_unknown_model_names_the_alternatives():
    with pytest.raises(ValueError, match="unknown model"):
        protocol.model("tip3p")


# --- GROMOS is refused for a model 54A7 does not carry --------------------------

def test_gromos_topology_is_refused_for_a_gromacs_only_model(tmp_path):
    """Before make_top is reached: there is no @solv block to ask for."""
    with pytest.raises(ValueError, match="no 54A7 solvent building block"):
        forcefield.build_topology("opc3", tmp_path)


# --- the GROMACS topology ------------------------------------------------------

def test_stock_models_include_the_force_field_itp(tmp_path):
    for model, itp in (("spc", "spc.itp"), ("spce", "spce.itp")):
        text = gromacs.write_topology(model, tmp_path).read_text()
        assert f'#include "{gromacs.FORCEFIELD}/{itp}"' in text
        assert not (tmp_path / itp).exists()   # taken from the force field, not copied


def test_opc3_gets_its_own_itp_and_not_the_spc_fallback(tmp_path):
    text = gromacs.write_topology("opc3", tmp_path).read_text()
    assert '#include "opc3.itp"' in text
    assert "spc.itp" not in text
    assert (tmp_path / "opc3.itp").exists(), "the .itp must sit beside the .top for scp"


def test_an_unknown_model_raises_instead_of_silently_becoming_spc(tmp_path):
    with pytest.raises(ValueError, match="unknown model"):
        gromacs.write_topology("tip4p", tmp_path)


def _itp_fields(text: str, section: str) -> list[list[str]]:
    rows, inside = [], False
    for line in text.splitlines():
        line = line.split(";")[0].strip()
        # "#ifndef FLEXIBLE" and friends are directives, not rows of the section
        # they sit in the middle of.
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            inside = line.replace("[", "").replace("]", "").strip() == section
            continue
        if inside:
            rows.append(line.split())
    return rows


def test_packaged_opc3_itp_matches_the_published_parameters():
    """C6 and C12 are re-derived from sigma and eps, which is where they came from.

    The stock file is in amber19sb.ff's comb-rule 2 (sigma/eps) form and cannot be
    included under gromos54a7.ff's comb-rule 1; restating it by hand is exactly the
    step where a transcription error would go unnoticed.
    """
    sigma, eps = 0.317427035094, 0.683690704000     # amber19sb.ff/ffnonbonded.itp
    text = gromacs.water_itp_source("opc3.itp").read_text()

    types = {row[0]: row for row in _itp_fields(text, "atomtypes")}
    assert float(types["OW_opc3"][-2]) == pytest.approx(4 * eps * sigma ** 6, rel=1e-9)
    assert float(types["OW_opc3"][-1]) == pytest.approx(4 * eps * sigma ** 12, rel=1e-9)
    # Hydrogen carries no Lennard-Jones term in any of these models.
    assert float(types["HW_opc3"][-2]) == 0.0
    assert float(types["HW_opc3"][-1]) == 0.0

    model = protocol.MODELS["opc3"]
    atoms = _itp_fields(text, "atoms")
    assert [float(row[6]) for row in atoms] == list(model.charges)
    # The topology's masses must be the ones box.py weights the centre of mass with.
    from water_benchmarking.box import MASSES
    assert [float(row[7]) for row in atoms] == pytest.approx(list(MASSES))

    settle = _itp_fields(text, "settles")[0]
    assert float(settle[2]) == pytest.approx(model.r_oh)
    assert float(settle[3]) == pytest.approx(model.r_hh)


def test_topology_declares_the_whole_box(tmp_path):
    """A count that disagrees with the .gro is a grompp error, not a silent one --
    but it is still worth pinning, since the box size is what sets the cutoff."""
    text = gromacs.write_topology("opc3", tmp_path).read_text()
    assert f"SOL                 {protocol.N_WATERS}" in text


# --- the two engines run production differently ---------------------------------

def test_gromacs_production_chains_rather_than_regenerating_velocities():
    """The shared ladder is written for GROMOS; GROMACS must not inherit this bit.

    GROMOS runs the ten production nanoseconds as independent replicates from one
    equilibrated box, so each draws fresh velocities.  GROMACS chains its segments
    from the previous .gro, so regenerating velocities would discard the state it
    is meant to be continuing -- ten times, at the same seed each time.  The .mdp
    files of the published SPC and SPC/E runs chain; the generator had drifted.
    """
    from water_benchmarking import imd

    production = [s for s in gromacs.stages() if s.name.startswith("md_")]
    assert len(production) == protocol.PRODUCTION_SEGMENTS
    assert not any(s.generate_velocities for s in production)

    # ... while GROMOS still runs them as replicates, which is the whole reason
    # the two differ here.
    gromos_production = [s for s in imd.run_ladder() if s.name.startswith("md_")]
    assert all(s.generate_velocities for s in gromos_production)


def test_equilibration_still_generates_velocities_once():
    """eq1 is where the velocities come from, for both engines."""
    eq1 = [s for s in gromacs.stages() if s.name == "eq1"]
    assert len(eq1) == 1 and eq1[0].generate_velocities


def test_production_mdp_continues_the_previous_segment(tmp_path):
    text = (tmp_path / "md_01.mdp")
    gromacs.write_all(tmp_path, "opc3")
    body = text.read_text()
    assert "gen_vel                  = no" in body
    assert "continuation             = yes" in body
    assert "gen_seed" not in body


def test_a_model_cannot_claim_gromos_without_a_building_block():
    """The registry must not be able to contradict itself."""
    with pytest.raises(ValueError, match="cannot be one of its engines"):
        protocol.Model(
            name="bogus", charges=(-0.8, 0.4, 0.4), r_oh=0.1, r_hh=0.1633,
            seed_index=9, gromos_block=None, engines=("gromos", "gromacs"),
        )


def test_analyse_prints_non_numeric_values_without_raising(capsys, monkeypatch):
    """`analyse` used to crash formatting dielectric_relation, a str, with ":.6g".

    It raised only after the entire trajectory analysis had run, and nothing is
    written to disk on the way, so a whole model's analysis was lost to a print.
    """
    import argparse
    from water_benchmarking import cli, report

    results = report.Results(
        model="opc3", engine="gromacs",
        values={"density": 994.4, "dielectric_relation": "conducting (eps = 1 + y)"},
    )
    # The analysis itself is tested elsewhere; this is about the printing.
    monkeypatch.setattr(cli, "analyse_run", lambda *a, **k: results)
    cli.cmd_analyse(argparse.Namespace(
        model="opc3", engine="gromacs", root=Path("/nonexistent"), collect=False,
    ))
    out = capsys.readouterr().out
    assert "994.4" in out and "conducting" in out


# --- the self-polarisation correction -------------------------------------------

def test_dipole_moments_match_the_published_values():
    """Derived from charges and geometry, so it cannot drift from what is simulated.

    Table I of Izadi & Onufriev 2016 gives SPC/E 2.35 D and OPC3 2.43 D.
    """
    from water_benchmarking.analysis import density_hov

    assert density_hov.dipole_moment("spce") == pytest.approx(2.35, abs=0.005)
    assert density_hov.dipole_moment("opc3") == pytest.approx(2.43, abs=0.005)


def test_polarisation_formula_reproduces_the_published_spce_value():
    """(mu - mu_gas)^2 / 2*alpha must give Berendsen 1987's 5.22 kJ/mol for SPC/E.

    That is the calibration: the same formula then gives OPC3's 7.03, which is what
    reconciles this protocol's raw 51.66 with the paper's 44.89.
    """
    from water_benchmarking.analysis import density_hov

    assert density_hov.self_polarisation_energy("spce") == pytest.approx(5.22, abs=0.05)
    assert density_hov.self_polarisation_energy("opc3") == pytest.approx(7.03, abs=0.05)
    assert density_hov.POLARISATION_CORRECTION["opc3"] == pytest.approx(
        density_hov.self_polarisation_energy("opc3"), abs=0.05
    )


def test_spc_is_deliberately_left_uncorrected():
    """Its dipole is enhanced too, but the published SPC values are uncorrected."""
    from water_benchmarking.analysis import density_hov

    assert "spc" not in density_hov.POLARISATION_CORRECTION
    assert density_hov.self_polarisation_energy("spc") > 3.0   # it would be 3.76
