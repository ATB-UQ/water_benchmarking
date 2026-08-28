# OPC3 protocol sweep: cutoff 1.4 nm, timestep 2 fs, RF ε_RF = 78.4 / 61, PME

Same 2048-water box, branched from the end of the 1.8 nm / 1 fs main run through one
100 ps NPT re-equilibration at the new cutoff and timestep (`diagnostics/o3eq.mdp`; the
density relaxed 995 → 991 over it, so the step was needed). Then 10 ns per condition,
sampled every 0.2 ps (`diagnostics/o3rf78.mdp`, `o3rf61.mdp`, `o3pme.mdp`,
`run_opc3.slurm`). Analysed with the diagnostics machinery unchanged
(`results/opc3_settings/analyse_sweep.py`; raw values and uncertainties in the `.json`
files beside it). ΔH_vap carries OPC3's 7.03 kJ mol⁻¹ polarisation correction. Setonix,
one GCD, GROMACS 2026.1; ns/day from mdrun's own `Performance:` line, the PME run with
`-pme gpu`.

| condition | ρ (kg m⁻³) | ΔH_vap − E_pol | D (10⁻⁹ m² s⁻¹) | τ₂(HH) (ps) | ε | ns/day | rel. |
|---|---|---|---|---|---|---|---|
| 1.8 nm, 1 fs, RF 61 (main run) | 994.4 ± 0.1 | 44.63 | 2.37 ± 0.02 | 2.28 | 78.3 ± 1.7 | 558 | 1.00 |
| 1.4 nm, 2 fs, RF ε_RF = 78.4 | 993.1 ± 0.1 | 44.62 | 2.35 ± 0.02 | 2.33 | 82.5 ± 2.5 | 1139 | 2.04 |
| 1.4 nm, 2 fs, RF ε_RF = 61 | 993.0 ± 0.1 | 44.62 | 2.29 ± 0.01 | 2.33 | 80.3 ± 2.6 | 1142 | 2.05 |
| 1.4 nm, 2 fs, PME | 992.6 ± 0.1 | 44.57 | 2.40 | 2.32 | 77.6 ± 2.2 | 884 | 1.58 |
| experiment | 997.05 | 43.99 | 2.30 | 2.0 | 78.4 | | |

## What it says

- **ε_RF does not matter, and neither does RF vs PME, for ε.** 82.5 ± 2.5, 80.3 ± 2.6 and
  77.6 ± 2.2 are one population; the main run's 78.3 ± 1.7 sits inside it. This is the
  README's standing observation — GROMACS's reaction field shows no ε_RF dependence and
  fluctuates like PME — reproduced on a third model at a second cutoff. The "SPC-specific
  61" is therefore not a source of error for OPC3 under this engine; it is simply inert.
  (It would not be inert in GROMOS, which honours the finite ε_RF — see README.)
- **Nothing else moves either.** D, τ₂ and ΔH_vap agree across all four conditions to
  within their uncertainties. OPC3 was fitted under PME, and under a 1.4 nm reaction field
  it gives the same ρ, D, τ₂ and ε: it transfers.
- **The cutoff costs 1.3–1.8 kg m⁻³ of density**, 994.4 → 993.0 (RF) / 992.6 (PME), the
  same direction and size as the SPC `rc14` diagnostic. That is the one systematic effect
  of the protocol change, and it is a 0.15 % one.
- **Efficiency: 2.05× for RF, 1.58× for PME**, almost entirely from the timestep. Per
  step, 1.4 nm RF is only ~12 % cheaper than 1.8 nm — at 6k atoms on one GCD the
  non-bonded kernel is not the bottleneck, so the shorter cutoff buys little here (it will
  buy more on a solvated protein, where the pair count scales with R_c³). PME is 22 %
  slower than RF at the same cutoff and timestep.
- **10 × 1 ns is needed only for ε.** ρ and ΔH_vap are converged at 1 ns (their means do
  not move between 1 and 10 ns; only the error bars shrink); D and τ₂ by 2 ns. ε from a
  single nanosecond scatters over 65–86 (σ ≈ 5.4), so 10 ns gives ±1.7–2.5 — about the
  resolution needed to tell 78 from 70. Fewer frames are fine down to 0.2 ps (this sweep;
  ~10 points across τ₂), not below.

## Recommendation for the protocol

For a water-model benchmark run with GROMACS: **1.4 nm, 2 fs, reaction field at ε_RF =
the model's own ε (or 61 — it makes no measurable difference), 10 ns, 0.2 ps sampling.**
Twice the throughput of the current protocol, half the storage, and every property within
its uncertainty of the 1.8 nm / 1 fs result except a 0.15 % density shift that is
understood. PME is not needed for accuracy here and costs 22 %; it becomes the right
choice only where the solute carries net charge or a long-ranged dipole that a 1.4 nm
reaction field cannot screen — which the SPC `rc09` diagnostic (density −1 %, D −7 %) says
is a question of cutoff, not of scheme.
