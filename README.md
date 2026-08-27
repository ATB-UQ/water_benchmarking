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

GROMACS, 10 ns per model. GROMOS: pending. Full tables and figures in [`results/`](results/).

| Property | SPC | SPC/E | Experiment |
|---|---|---|---|
| Density (kg m⁻³) | 976.0 | 997.4 | 997.05 |
| ΔH_vap (kJ mol⁻¹) | 44.19 | 49.27 (44.05 with the SPC/E polarisation correction) | 43.99 |
| D (10⁻⁹ m² s⁻¹, Yeh–Hummer corrected) | 4.31 | 2.75 | 2.30 |
| τ₂(HH) (ps) | 1.17 | 1.99 | 2.0 |
| ε | 66.3 ± 1.7 | 69.4 ± 1.5 | 78.4 |

Both models behave as published: SPC/E reproduces density and (corrected) ΔH_vap to within
0.1 %, diffuses 20 % too fast and reorients at the experimental rate; SPC is 2 % light, diffuses
nearly twice too fast and reorients in half the time. Both underestimate ε by 12–15 %.
Statistical errors are below the last digit shown except where given.

**The dielectric constant is reported through ε = 1 + y**, where y = (⟨M²⟩−⟨M⟩²)/(3ε₀VkT).
The textbook finite-ε_RF relation (Neumann) gives ~140 for the same trajectories: at ε_RF = 61 it
has a pole at y = 123, fifty units from where water sits, and it presumes an ε_RF dependence of
the fluctuation that the force field cannot carry (k_rf at ε_RF = 61 is within 2.4 % of the
conducting-boundary value). Controls varying thermostat, cutoff (0.9–1.8 nm), box (2048–16384
waters) and boundary condition (ε_RF = 61, ∞, PME) all give y = 59–76; see
[`results/diagnostics.md`](results/diagnostics.md). Anyone computing ε from an ATB reaction-field
trajectory with the finite-ε_RF formula will get roughly double the true value.

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
