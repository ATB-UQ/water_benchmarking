# water_benchmarking

Benchmarks the two classical water models ATB relies on — **SPC** and **SPC/E** — against
experiment, under the exact simulation protocol the peptide/protein validation campaign uses.
Every property is computed with both **GROMOS** (md++ 1.6.0, MPI, on Gadi) and **GROMACS**
(GPU, on Setonix), so a disagreement can be attributed to the engine rather than to the setup.

Five properties, all at 298.15 K and 1 atm: density, heat of vaporisation, self-diffusion,
rotational correlation times, and the static dielectric constant.

## Why it exists

The peptide runs put ATB force fields in SPC water at a 1.8 nm single-range cutoff with a
reaction field. Nothing had ever measured what that water model does *on its own* under that
protocol — the pure-liquid machinery in `gromos_job_wrapper` covers density and heat of
vaporisation for ATB solutes only, and nothing in the platform computed diffusion, rotational
relaxation or a dielectric constant at all.

## The protocol

Defined once, in [`protocol.py`](src/water_benchmarking/protocol.py), and rendered into both a
GROMOS `.imd` and a GROMACS `.mdp` so the two engines cannot drift apart.

| | |
|---|---|
| System | 2048 water molecules, cubic, edge ≈ 3.98 nm |
| Temperature / pressure | 298.15 K, 1 atm (Berendsen, τ_T 0.1 ps, τ_P 0.5 ps) |
| Timestep | 1 fs |
| Cutoff | single-range 1.8 nm (RCUTP = RCUTL = RCRF = ASHAPE) |
| Electrostatics | reaction field, ε_RF = 61 |
| Constraints | rigid water (SHAKE in GROMOS, SETTLE in GROMACS) |
| Ladder | emin → 10 ps NVT at 50 K → 20 ps NVT at 298 K → 100 ps NPT → 10 × 1 ns NPT |
| Sampling | 0.1 ps (τ₂ of water is ≈ 2 ps, so the 50 ps of the peptide runs is far too coarse) |

**2048 is not a free choice.** 1024 waters at bulk density give a ≈ 3.1 nm, and a half-box of
1.55 nm is *below* the 1.8 nm cutoff, which breaks the minimum image convention. 2048 is the
smallest power of two that clears it (half-box 1.99 nm).

**ε_RF is 61 for both models**, because that is what the peptide protocol uses. SPC/E is more
often published with ε_RF ≈ 71; the point here is to benchmark the protocol ATB actually runs.

## What Step 0 established

Four things had to be checked before any of this could work, and each one changed the design:

1. **`make_top` needs no `@seq` at all.** An empty sequence is rejected outright ("Cannot find
   building block for"), but omitting it gives exactly the wanted topology — `SOLUTEATOM` with
   NRP 0 plus the requested solvent block.
2. **md++ accepts `NPM 0`, but not with solute SHAKE.** `NTC 2` fails with *"solvent only
   simulation does not work with SHAKE for solute"*. The protocol uses `NTC 1`; nothing is lost,
   because the rigid water geometry comes from `SOLVENTCONSTR` via `NTCS` either way.
3. **The `ene_ana` library must be the one shipped with this md++ build**
   (`/opt/gromos/1.6.0/share/md++/ene_ana.md++.lib`). The copy in `gromos_job_wrapper` carries the
   same `ENEVERSION` stamp (2023-04-15) and still fails to parse this build's `.tre`
   ("Tried to read an integer for NUM_EDS_STATES"). The version stamp does not identify the layout.
4. **gromos++ `epsilon`, `diffus` and `check_box` cannot be used here at all.** All three gather
   the trajectory first, and every gather method requires at least one solute molecule:
   *"the cog gather method requires at least one solute molecule"*. In a solvent-only system they
   refuse to run. This is why the transport and dielectric analyses are implemented in this package
   rather than shelled out — and why **GROMACS is the cross-check**, not gromos++.

## Cross-engine validation

Single-point potential energy on the identical 2048-molecule SPC configuration:

| engine | U_pot (kJ mol⁻¹) |
|---|---|
| GROMOS md++ 1.6.0 | −70875.4 |
| GROMACS 2026.1 | −70884.8 |

0.013 %, or 0.005 kJ mol⁻¹ per molecule — the two engines see the same Hamiltonian. The residual
is the Verlet cluster list against GROMOS's exact cutoff.

`vdw-modifier = none` is load-bearing: the default `Potential-shift` changes the reported LJ
energy, and the heat of vaporisation is computed straight from it.

**The Verlet buffer must be left on.** Pinning `rlist` to the cutoff
(`verlet-buffer-tolerance = -1`) looks like the closer match to GROMOS's exact cutoff, but it is
not "more faithful" -- it means pairs drifting inside the cutoff between list updates are simply
missed, and eq1 died of an unsettleable water at step 38. With the default tolerance GROMACS picks
`rlist = 1.802` nm, a 2 pm buffer that only widens the neighbour *list*; interactions are still cut
at 1.8 nm and the single-point energy above is unchanged to the last digit.

## The analyses

All four are in [`analysis/`](src/water_benchmarking/analysis/) and each is validated in
`tests/test_analysis.py` against a case with a known answer.

- **Density / ΔH_vap** — `ene_ana`. For a rigid non-polarisable model U_gas = 0, so
  ΔH_vap = −⟨U_pot⟩/N + RT and no vacuum leg is needed. SPC/E additionally reports the value with
  Berendsen's self-polarisation correction (−5.22 kJ mol⁻¹), which is how SPC/E is normally quoted.
- **Diffusion** — MSD by the Fast Correlation Algorithm (FFT; the direct form is O(n²) and hopeless
  at 10⁵ frames), Einstein fit over 10–100 ps. Reported both as simulated and with the Yeh–Hummer
  finite-size correction, which is ≈ +7 % for this box and must never be hidden.
- **Rotational relaxation** — C₁ and C₂ for the O–H, H–H and dipole vectors. C₂ is obtained from
  the autocorrelation of the second-rank tensor u⊗u, since ⟨(u(0)·u(t))²⟩ is not an
  autocorrelation of u itself.
- **Dielectric** — ⟨M²⟩−⟨M⟩² with the Neumann reaction-field relation, *not* the vacuum Kirkwood
  formula (which would understate ε by ~15 % at ε_RF = 61). A running estimate is reported because
  ε converges slowly: one box gives one sample per frame regardless of how many molecules it holds.

## Usage

```bash
water-bench build                          # topology, box and inputs for both engines
water-bench run --model spc                # GROMOS ladder on Gadi (blocking, one chain)
water-bench run --model spce               # run the second model concurrently
water-bench analyse --model spc            # properties for one finished run
water-bench report                         # the comparison table across models and engines
```

Runs land under `/ssd1_nas_md/water_benchmarking/<model>/<engine>/`; trajectories stay there and
only the summary tables and plots are committed.

GROMOS goes to Gadi through `gromos_job_wrapper/deployment/gadi_md.sh`, which stages inputs,
sizes the PBS walltime from `NSTLIM`, and pulls results back with size verification. Walltime is
sized from `INITIAL_SECONDS_PER_STEP = 0.02` and then refined from the measured eq3 rate — the
shim's own default of 0.08 s/step is the 23k-atom peptide figure and would request ~22 h for a
system a quarter that size.

GROMACS runs on one Setonix GCD per model (a whole GPU node is *slower* at 6k atoms), with
`-nb gpu -bonded gpu` and no PME offload — `-pme gpu` is a fatal error under a reaction field.

## The dielectric constant: an open finding (2026-08-27)

Four of the five properties come out where the literature says they should.
The fifth does not, and the discrepancy is in the *simulation*, not the analysis:

| SPC, 2048 waters, 298 K | ε |
|---|---|
| this protocol (RF ε_RF = 61, R_c = 1.8 nm), GROMACS, 10 ns | **140** |
| same, Berendsen → v-rescale / C-rescale, 1 ns | 141 |
| same, reaction field → PME, 1 ns | 76 |
| published SPC (RF or Ewald, R_c ≈ 0.9–1.4 nm) | 65 ± 5 |

SPC/E gives 155 under the protocol against a published ~71 — the same ×2.2.

What is ruled out: the analysis (`gmx dipoles -epsilonRF 61` agrees with the
in-house code to 5 % on the same frames; the molecular dipole is 2.274 D
exactly); the thermostat (v-rescale gives the same answer); the state point
(298.1 K, 2 bar, 976 kg m⁻³); and sampling length — the box dipole decorrelates
in 6–8 ps, so each 1 ns segment carries ~100 independent samples, the
per-segment values scatter 137–148, and the cumulative estimate is flat from
5 ns on. Truncation would bias ε *low* in any case.

What is implicated: the reaction field at this geometry. The RF box's ⟨ΔM²⟩ is
only 13 % below the PME box's, where a true ε ≈ 70 under Neumann's relation
requires ~36 % below. R_c/L here is 1.8/3.98 = 0.45, hard against the L/2
limit; the SPC/RF literature sits near 0.3. Diagnostics in flight: the same box
at R_c = 1.4 and 0.9 nm, and the 1.8 nm cutoff in an 8× box (R_c/L = 0.23).
If ε falls toward 70 as R_c/L falls, the number is an artefact of running a
1.8 nm reaction field in a box only just big enough to hold it — which the
peptide boxes (minwall 2.6 nm, L ≳ 6 nm) do not do, but any small solvated
system under this protocol would.
