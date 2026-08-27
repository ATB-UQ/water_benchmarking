"""Each analysis is checked against a case whose answer is known analytically.

These are the tests that matter: the trajectory analyses produce numbers nobody
can eyeball, so the only defence against a wrong factor or a missed unit
conversion is a synthetic system with a known result.
"""
import numpy as np
import pytest

from water_benchmarking.analysis import dielectric, diffusion, errors, rotation


# --- diffusion -----------------------------------------------------------------

def test_msd_recovers_a_known_diffusion_coefficient():
    rng = np.random.default_rng(0)
    n_frames, n_molecules, dt = 4000, 200, 0.1
    d_true = 2.5e-9                       # m^2 s^-1
    sigma = np.sqrt(2 * (d_true * 1e18 / 1e12) * dt)   # nm per dimension per step

    steps = rng.normal(0, sigma, (n_frames - 1, n_molecules, 3))
    path = np.concatenate([np.zeros((1, n_molecules, 3)), np.cumsum(steps, axis=0)])

    msd = diffusion._msd_fft(path)
    lags = np.arange(n_frames) * dt
    result = diffusion.diffusion_from_msd(lags, msd, edge=4.0)
    assert result.d_pbc == pytest.approx(d_true, rel=0.05)


def test_yeh_hummer_correction_is_positive_and_shrinks_with_box_size():
    """A bigger box needs a smaller correction; both must raise D, never lower it."""
    lags = np.arange(2000) * 0.1
    msd = 6 * (2.5e-9 * 1e18 / 1e12) * lags        # perfectly linear
    small = diffusion.diffusion_from_msd(lags, msd, edge=3.0)
    large = diffusion.diffusion_from_msd(lags, msd, edge=6.0)
    assert small.d_corrected > small.d_pbc
    assert (small.d_corrected - small.d_pbc) > (large.d_corrected - large.d_pbc)


def test_unwrap_follows_a_molecule_across_the_boundary():
    edge = 4.0
    # A molecule stepping steadily right, wrapped back into the box each time.
    true_path = np.arange(0, 12, 0.5).reshape(-1, 1, 1) * np.array([[[1.0, 0.0, 0.0]]])
    wrapped = np.mod(true_path, edge)
    edges = np.full(len(wrapped), edge)
    assert np.allclose(diffusion.unwrap(wrapped, edges), true_path, atol=1e-9)


# --- rotation ------------------------------------------------------------------

def test_frozen_vectors_stay_perfectly_correlated():
    rng = np.random.default_rng(1)
    vector = rng.normal(size=(1, 300, 3))
    vector /= np.linalg.norm(vector, axis=2, keepdims=True)
    frozen = np.repeat(vector, 200, axis=0)

    c1, c2 = rotation.correlation_functions(frozen, 50)
    assert np.allclose(c1, 1.0, atol=1e-6)
    assert np.allclose(c2, 1.0, atol=1e-6)


def test_rotational_diffusion_gives_the_textbook_tau1_over_tau2_ratio():
    """For isotropic rotational diffusion tau_1 / tau_2 = 3 exactly."""
    rng = np.random.default_rng(2)
    n_frames, n_molecules, dt = 4000, 300, 0.05
    tau2_target = 2.0
    rotational_d = 1.0 / (6.0 * tau2_target)
    sigma = np.sqrt(2 * rotational_d * dt)

    u = rng.normal(size=(n_molecules, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    trajectory = np.empty((n_frames, n_molecules, 3), dtype=np.float32)
    trajectory[0] = u
    for frame in range(1, n_frames):
        u = u + rng.normal(0, sigma, (n_molecules, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        trajectory[frame] = u

    c1, c2 = rotation.correlation_functions(trajectory, int(15 / dt))
    lags = np.arange(len(c2)) * dt

    def decay_time(correlation):
        window = (lags >= 1.0) & (lags <= 8.0) & (correlation > 0)
        return -1.0 / np.polyfit(lags[window], np.log(correlation[window]), 1)[0]

    assert decay_time(c2) == pytest.approx(tau2_target, rel=0.10)
    assert decay_time(c1) / decay_time(c2) == pytest.approx(3.0, rel=0.10)


def test_molecular_vectors_are_unit_length_and_dipole_bisects():
    positions = np.array([[[0.0, 0.0, 0.0],
                           [0.1, 0.0, 0.0],
                           [-0.0333, 0.0943, 0.0]]])
    vectors = rotation.molecular_vectors(positions)
    for name, vector in vectors.items():
        assert np.linalg.norm(vector[0]) == pytest.approx(1.0), name
    # The bisector is equidistant in angle from both O-H bonds.
    oh1 = vectors["OH"][0]
    oh2 = (positions[0, 2] - positions[0, 0]) / np.linalg.norm(positions[0, 2] - positions[0, 0])
    dipole = vectors["dipole"][0]
    assert np.dot(dipole, oh1) == pytest.approx(np.dot(dipole, oh2), abs=1e-3)


# --- dielectric ----------------------------------------------------------------

def test_reaction_field_relation_inverts_exactly():
    for epsilon in (1.0, 20.0, 65.0, 78.4):
        y = (epsilon - 1) * (2 * 61 + 1) / (2 * 61 + epsilon)
        assert dielectric._solve_epsilon(y, 61.0) == pytest.approx(epsilon)


def test_a_frozen_box_has_no_dipole_fluctuation_and_epsilon_one():
    """Zero fluctuation must give eps = 1, not a division by zero."""
    assert dielectric._solve_epsilon(0.0, 61.0) == pytest.approx(1.0)


def test_total_dipole_is_charge_weighted_and_neutral_under_translation():
    rng = np.random.default_rng(3)
    positions = rng.normal(size=(50, 3, 3))
    charges = np.array([-0.82, 0.41, 0.41])

    dipole = dielectric.total_dipole(positions, charges)
    # A neutral system's dipole must not depend on where the origin is.
    shifted = dielectric.total_dipole(positions + np.array([1.0, -2.0, 0.5]), charges)
    assert np.allclose(dipole, shifted, atol=1e-9)


# --- errors --------------------------------------------------------------------

def test_block_averaging_reports_a_larger_error_for_correlated_data():
    rng = np.random.default_rng(4)
    independent = rng.normal(0, 1, 8192)
    # An AR(1) series with the same variance but a long correlation time.
    correlated = np.empty(8192)
    correlated[0] = rng.normal()
    for i in range(1, 8192):
        correlated[i] = 0.98 * correlated[i - 1] + rng.normal(0, 0.2)

    assert errors.block_average(correlated).error > errors.block_average(independent).error


# --- trajectory conversion -----------------------------------------------------

def test_broken_molecules_are_caught_not_silently_analysed():
    """A trajectory converted without `-pbc mol` must fail, not produce numbers.

    GROMACS wraps atoms individually, splitting ~80 of 2048 molecules across the
    boundary in every frame. Nothing downstream raises on that -- the dipole and
    the molecular vectors just come out wrong -- so the guard is the only defence.
    """
    from water_benchmarking import gromacs
    from water_benchmarking.trc import Frame

    whole = np.tile(np.array([[0.10, 0.10, 0.10],
                              [0.20, 0.10, 0.10],
                              [0.15, 0.19, 0.10]]), (8, 1, 1))
    gromacs.assert_whole_molecules(Frame(0.0, 0, whole, 4.0))

    split = whole.copy()
    split[3, 1, 0] += 4.0          # one hydrogen wrapped to the far face
    with pytest.raises(AssertionError, match="pbc mol"):
        gromacs.assert_whole_molecules(Frame(0.0, 0, split, 4.0))


def test_dielectric_reports_the_stable_relation_and_flags_the_unstable_one():
    """y = 67 is literature SPC under eps = 1 + y and eps = 150 under Neumann(61).

    The two must both be reported, and the Neumann sensitivity must say the
    second number is not to be trusted: near its pole a 10% change in y moves
    eps by fifty.
    """
    rng = np.random.default_rng(5)
    n = 4000
    volume = 62.78
    # Draw box dipoles with the variance that gives y = 67 at this V and T.
    target_y = 67.0
    variance_si = target_y * 3 * dielectric.VACUUM_PERMITTIVITY * volume * 1e-27 \
        * dielectric.BOLTZMANN * 298.15
    sigma_e_nm = np.sqrt(variance_si / 3) / (dielectric.ELEMENTARY_CHARGE * 1e-9)
    dipoles = rng.normal(0, sigma_e_nm, (n, 3))

    result = dielectric.from_dipoles(dipoles, np.full(n, volume), np.arange(n) * 0.1)
    assert result.y == pytest.approx(target_y, rel=0.08)
    assert result.epsilon == pytest.approx(1 + target_y, rel=0.08)
    assert result.epsilon_neumann > 120
    assert result.neumann_sensitivity > 2.5
