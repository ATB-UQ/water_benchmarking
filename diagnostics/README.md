# Dielectric-constant diagnostics

Six 1 ns GROMACS runs of 2048 SPC waters (16384 for the 8x box) on Setonix, each changing
exactly one thing relative to the benchmark protocol, made to find out why the protocol
gave eps ~ 140 against a published ~65. All start from the protocol's own eq3 structure.
The `.mdp` files here are the ones that ran; `run*.slurm` are the jobs.

`manifest.json` is what `water-bench diagnostics` reads to analyse them with exactly the
code used for the main runs, so every row of the diagnostics table carries all five
properties rather than only the one that prompted the runs.
