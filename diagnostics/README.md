# Dielectric-constant diagnostics

Eight 1 ns GROMACS runs of 2048 SPC waters (16384 for the 8x box) on Setonix, each changing
exactly one thing relative to the benchmark protocol, made to find out why the protocol
gave eps ~ 140 against a published ~65. All start from the protocol's own eq3 structure.
The `.mdp` files here are the ones that ran; `run*.slurm` are the jobs.

`manifest.json` is what `water-bench diagnostics` reads to analyse them with exactly the
code used for the main runs, so every row of the diagnostics table carries all five
properties rather than only the one that prompted the runs.

The last two run the reaction field on the CPU nonbonded kernel instead of the GPU one, to
find out whether the missing eps_rf dependence of the box-dipole fluctuation in the GPU runs
(GROMOS shows the dependence; GROMACS-GPU does not) is a property of the GPU kernel or of
GROMACS's reaction field as such.
