# water_benchmarking

Benchmarks **SPC** and **SPC/E** water against experiment under the simulation protocol the
ATB peptide/protein validation uses, with both **GROMOS** (md++ 1.6.0, Gadi) and **GROMACS**
(2026.1, Setonix GPU) so that a disagreement can be attributed to the engine. Five properties
at 298.15 K and 1 atm: density, heat of vaporisation, self-diffusion, rotational correlation
time, static dielectric constant.

## Protocol

Defined once in [`protocol.py`](src/water_benchmarking/protocol.py) and rendered into both the
GROMOS `.imd` and the GROMACS `.mdp`.

| | |
|---|---|
| System | 2048 waters, cubic, a ≈ 3.98 nm (1024 would put the half-box below the cutoff) |
| Thermostat / barostat | Berendsen, 298.15 K / 1 atm, τ_T 0.1 ps, τ_P 0.5 ps |
| Timestep | 1 fs; rigid water (SHAKE / SETTLE) |
| Cutoff | single-range 1.8 nm |
| Electrostatics | reaction field, ε_RF = 61 (the peptide-protocol value, for both models) |
| Ladder | emin → 10 ps NVT 50 K → 20 ps NVT → 100 ps NPT → 10 × 1 ns NPT |
| Sampling | 0.1 ps |

Production is ten 1 ns runs per model. GROMACS chains them; GROMOS runs them as concurrent
replicates from the same equilibrated box with distinct velocity seeds, which finishes in the
wall time of one.

## Results

10 ns per model and engine (GROMOS: ten 1 ns replicates; GROMACS: ten chained 1 ns segments).
Full tables, per-run details and figures in [`results/`](results/).

| Property | SPC GROMOS | SPC GROMACS | SPC/E GROMOS | SPC/E GROMACS | Experiment |
|---|---|---|---|---|---|
| Density (kg m⁻³) | 974.8 | 976.0 | 996.4 | 997.4 | 997.05 |
| ΔH_vap (kJ mol⁻¹) | 44.16 | 44.19 | 49.24 | 49.27 | 43.99 |
| ΔH_vap, SPC/E polarisation-corrected | – | – | 44.02 | 44.05 | 43.99 |
| D (10⁻⁹ m² s⁻¹, Yeh–Hummer corrected) | 4.45 | 4.31 | 2.78 | 2.75 | 2.30 |
| τ₂(HH) (ps) | 1.16 | 1.17 | 1.96 | 1.98 | 2.0 |
| ε | 68 ± 2 | 66 ± 2 | 70 ± 2 | 69 ± 2 | 78.4 |

The two engines agree to within 0.1 % on density and ΔH_vap, 1–3 % on D and τ₂, and 2–3 % on ε.
Both models behave as published: SPC/E reproduces density and (corrected) ΔH_vap to within 0.1 %,
diffuses 20 % too fast and reorients at the experimental rate; SPC is 2 % light, diffuses nearly
twice too fast and reorients in half the time. Both underestimate ε by 11–15 %. Statistical errors
are below the last digit shown except where given.

**The dielectric constant is read through the boundary condition each engine realises.** What is
measured is the dimensionless box-dipole fluctuation y = (⟨M²⟩−⟨M⟩²)/(3ε₀VkT); the relation that
turns it into ε depends on the electrostatic boundary. GROMOS honours the finite ε_RF = 61: its
fluctuation is 35 % below the conducting-boundary value, exactly as Neumann's relation
(ε−1)(2ε_RF+1)/(2ε_RF+ε) = y predicts for ε ≈ 70, and that relation gives 68 / 70. GROMACS's
reaction field shows no ε_RF dependence — ε_RF = 61 and ε_RF = ∞ fluctuate alike, on the GPU and
CPU kernels both, as does PME — so the conducting-boundary relation ε = 1 + y describes it, giving
66 / 69. Why GROMACS's reaction field does not realise the finite-ε_RF boundary that GROMOS's does
is not resolved here; the natural next test is GROMOS at ε_RF → ∞. Each read through the
other's relation is wrong by a factor of ~2 (GROMOS 44, GROMACS 140), and the Neumann relation is
in any case unstable near y ≈ 65 (a pole at y = 123, dε/dy ≈ 2.5–5). `results/summary.md` records y
and both readings beside every number. Controls varying thermostat, cutoff (0.9–1.8 nm) and box
(2048–16384 waters) leave y unchanged within scatter — and leave density, ΔH_vap, D and τ₂ within
1 % — except a 0.9 nm cutoff, which lowers the density 1 %, D 7 % and slows reorientation 8 %;
see [`results/diagnostics.md`](results/diagnostics.md).

## What a pure-solvent system needs

- `make_top` with no `@seq` (an empty one is rejected); `NPM 0` with `NTC 1` — md++ refuses solute
  SHAKE with no solute, and the rigid geometry comes from `SOLVENTCONSTR` regardless.
- The `ene_ana` library shipped with the md++ build; the `gromos_job_wrapper` copy shares its
  `ENEVERSION` stamp and cannot parse the `.tre`.
- gromos++ `epsilon`, `diffus` and `check_box` gather first, and every gather method needs a
  solute; none run on a solvent-only system. The transport and dielectric analyses are in this
  package, and GROMACS is the cross-check (single-point energy agrees with GROMOS to 0.013 %).
- GROMACS: `vdw-modifier = none` (the default shift changes the potential energy and so ΔH_vap);
  keep the Verlet buffer (no buffer loses pairs and SETTLE fails); `-nb gpu` only — `-pme gpu` is
  fatal under a reaction field and `-bonded gpu` fatal with no bonded terms; convert with
  `trjconv -pbc mol` or ~80 molecules per frame straddle the boundary.
- Analyses are streamed one segment at a time and validated against cases with known answers
  (`tests/`).

## Usage

```bash
water-bench build                  # topology, box and inputs for both engines
water-bench run --model spc        # GROMOS ladder on Gadi (equilibration, then replicate fan)
water-bench submit-gromacs --model spc
water-bench analyse --model spc --engine gromacs
water-bench diagnostics            # the protocol-variation runs
water-bench report                 # results/summary.{md,csv} and figures
```

Runs live under `/ssd1_nas_md/water_benchmarking/<model>/<engine>/`.
