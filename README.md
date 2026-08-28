# water_benchmarking

Benchmarks **SPC**, **SPC/E** and **OPC3** water against experiment under the simulation protocol
the ATB peptide/protein validation uses. SPC and SPC/E are run with both **GROMOS** (md++ 1.6.0,
Gadi) and **GROMACS** (2026.1, Setonix GPU) so that a disagreement can be attributed to the
engine; **OPC3 is run with GROMACS only** — it has no 54A7 solvent building block, and the
two-engine comparison has already established that the engines agree. Five properties at
298.15 K and 1 atm: density, heat of vaporisation, self-diffusion, rotational correlation time,
static dielectric constant.

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

## Models

| Model | q_O / q_H (e) | r_OH / r_HH (nm) | Engines | Source |
|---|---|---|---|---|
| SPC | −0.82 / +0.41 | 0.1 / 0.163299 | GROMOS, GROMACS | 54A7 `H2O` |
| SPC/E | −0.8476 / +0.4238 | 0.1 / 0.163299 | GROMOS, GROMACS | 54A7 `H2OE` |
| OPC3 | −0.895170 / +0.447585 | 0.0978882 / 0.1598507 | GROMACS | Izadi & Onufriev 2016, *J. Chem. Phys.* **145**:074501 |

Defined once in [`protocol.py`](src/water_benchmarking/protocol.py) as a `Model` registry, which is
what says which engines a model is run with, what charges weight its box dipole and which seeds it
draws. The seed index is a fixed field rather than a position in the sorted registry, so adding a
model cannot renumber the seeds of runs already published.

**OPC3 is the one model whose parameters this package ships itself**, as
[`data/opc3.itp`](src/water_benchmarking/data/opc3.itp). GROMACS does carry it, in `amber19sb.ff`,
but that force field declares `comb-rule 2` (σ/ε) with `fudgeQQ 0.8333` while every topology here is
built on `gromos54a7.ff`, which declares `comb-rule 1` (C6/C12) — grompp refuses the mixture. For a
box of nothing but water there is a single interacting pair, so the packaged file simply restates
the same parameters in C6/C12 (C6 = 4εσ⁶, C12 = 4εσ¹², re-derived and checked in
`tests/test_models.py`) and the .top is otherwise unchanged. The file is copied beside the .top and
staged to Setonix with it, since Setonix's GROMACS has no copy of it.

OPC3 was fitted under Ewald/PME while this protocol runs a 1.8 nm reaction field at ε_RF = 61, so
its published ε and D are not exactly the values to expect here — both are boundary- and
cutoff-sensitive, see [`results/diagnostics.md`](results/diagnostics.md) — and `experiment.py`
carries correspondingly wide literature ranges. Table III of the paper gives, at 298.16 K and
1 bar: ρ = 0.996 ± 0.001 g cm⁻³, ε = 78.4 ± 1, D = 2.30 ± 0.02 × 10⁻⁹ m² s⁻¹, ΔH_vap = 10.73 ±
0.004 kcal mol⁻¹ (44.89 kJ mol⁻¹). It reports no rotational correlation time, so τ₂ is the one
property with no published value to check against; the deviation table marks it `[?]` rather than
leaving it looking confirmed.

## Results

10 ns per model and engine (GROMOS: ten 1 ns replicates; GROMACS: ten chained 1 ns segments).
Full tables, per-run details and figures in [`results/`](results/).

| Property | SPC GROMOS | SPC GROMACS | SPC/E GROMOS | SPC/E GROMACS | OPC3 GROMACS | Experiment |
|---|---|---|---|---|---|---|
| Density (kg m⁻³) | 974.8 | 976.0 | 996.4 | 997.4 | 994.4 | 997.05 |
| ΔH_vap (kJ mol⁻¹) | 44.16 | 44.19 | 49.24 | 49.27 | 51.66 | 43.99 |
| ΔH_vap, polarisation-corrected | – | – | 44.02 | 44.05 | 44.63 | 43.99 |
| D (10⁻⁹ m² s⁻¹, Yeh–Hummer corrected) | 4.45 | 4.31 | 2.78 | 2.75 | 2.37 | 2.30 |
| τ₂(HH) (ps) | 1.16 | 1.17 | 1.96 | 1.98 | 2.28 | 2.0 |
| ε | 68 ± 2 | 66 ± 2 | 70 ± 2 | 69 ± 2 | 78.3 | 78.4 |

**OPC3 closes the dielectric gap the other two share.** ε = 78.3 against experiment's 78.4, where
SPC/E gives 69 and SPC 66 — an 11–15 % shortfall in the property that screens charged residues.
D is 3 % high where SPC/E is 20 % high and SPC nearly double; density and corrected ΔH_vap are
within 1.4 %. It reproduces its own published ε (78.4 ± 1) to 0.1 % despite having been fitted
under PME at an 8 Å cutoff with a dispersion correction, and being run here under a 1.8 nm
reaction field at ε_RF = 61 — a boundary condition it never saw. τ₂ is the one property the paper
does not report, so it is unchecked against the model rather than confirmed (`[?]` in the
deviation table); at +14 % it reorients slightly too slowly, against SPC's −42 %.

The polarisation correction applies to SPC/E **and OPC3**, and is derived rather than tabulated:
E_pol = (μ − μ_gas)²/2α with μ computed from each model's own charges and geometry, giving 5.24
for SPC/E (against Berendsen's published 5.22) and 7.03 for OPC3. Izadi & Onufriev apply it too —
their SPC/E entry of 10.43 kcal mol⁻¹ (43.64 kJ mol⁻¹) is a corrected value, not the ~49 a raw
run gives. SPC is deliberately left uncorrected, since the published SPC numbers it is checked
against are.

### The polarisation correction: for and against

The full first-principles derivation — why a rigid model's gas is the unphysical state, where the
linear-response factor of ½ comes from, and what the argument does not establish — is in
[`docs/polarisation_correction.md`](docs/polarisation_correction.md). The short version:

A fixed-charge model carries a liquid-phase dipole (SPC/E 2.35 D, OPC3 2.43 D) that is larger
than the gas-phase molecule's 1.85 D. In the real liquid that enhancement is induced, and the
energy to induce it — the self-polarisation energy E_pol = (μ − μ_gas)²/2α — is paid on
condensation and recovered on vaporisation. A rigid model never pays it: its molecules leave the
liquid still carrying the enhanced dipole, so its raw ΔH_vap = −U_liq + RT contains a binding
energy the real liquid does not release. Berendsen, Grigera & Straatsma (*J. Phys. Chem.* 1987,
**91**, 6269) introduced the term with SPC/E, fitting the model so that ΔH_vap *after* subtracting
E_pol matches experiment; that is what let SPC/E take the larger dipole (and the better density,
diffusion and dielectric constant) that SPC could not.

*The case for.* It is the thermodynamically consistent comparison: the model's liquid-phase
energy is the quantity that should match ΔH_vap − E_pol, not ΔH_vap. Fitting a rigid model to the
raw ΔH_vap forces a dipole too small for the liquid — this is exactly why SPC and TIP3P underbind,
diffuse too fast and (for SPC) read 15 % low on ε. Vega & Abascal (*Phys. Chem. Chem. Phys.* 2011,
**13**, 19663) make the same argument for TIP4P/2005, whose raw ΔH_vap is ~1.5 kcal mol⁻¹ over
experiment by design; and Leontyev & Stuchebrukhov (*J. Chem. Phys.* 2009, **130**, 085102;
*Phys. Chem. Chem. Phys.* 2011, **13**, 2613) generalise it: a non-polarisable force field is
implicitly a molecule in an electronic continuum, and its energies must be read through that
continuum.

*The case against.* The correction is not a property of the model but of two external numbers,
μ_gas and α, and of the assumption that the whole dipole enhancement is inductive and linear.
Different choices of α (1.44 versus 1.47 Å³) or of the liquid-phase reference dipole move it by
several tenths of a kJ mol⁻¹, and it says nothing about the higher multipoles that also change
on condensation. A raw ΔH_vap is at least an unambiguous property of the model as simulated.
More practically, most force fields were fitted to the **uncorrected** ΔH_vap — SPC, TIP3P,
GROMOS 53A6/54A7 (Oostenbrink *et al.*, *J. Comput. Chem.* 2004, **25**, 1656), OPLS and the
organic-liquid benchmarks built on them (Caleman *et al.*, *J. Chem. Theory Comput.* 2012, **8**,
61) — so applying the correction post hoc to those breaks the comparison the parameterisation
was validated against. That is why this benchmark applies it per model, following each model's
own convention, rather than everywhere.

*Is water special?* No — but its numbers are. E_pol scales as (Δμ)²/α, and water combines an
unusually large enhancement (Δμ ≈ 0.5–0.6 D) with a small polarisability (1.44 Å³), giving
5–7 kJ mol⁻¹, i.e. 12–16 % of ΔH_vap. For most organic liquids Δμ is smaller and α several times
larger, so the term drops to a fraction of a kJ mol⁻¹ and is routinely ignored; only strongly
hydrogen-bonded, small, polar molecules (methanol, amides) reach a few kJ mol⁻¹, and even there
the force fields in common use were fitted without it. For a protein force field the question
does not arise at all: the solute is never vaporised, and what matters is that the *water* model's
own convention is used consistently with the solute–water parameters that were tuned to it.

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
water-bench build                  # topology, box and inputs for every engine each model uses
water-bench build --models opc3    # ... or one model; OPC3 skips GROMOS and says so
water-bench run --model spc        # GROMOS ladder on Gadi (equilibration, then replicate fan)
water-bench submit-gromacs --model opc3
water-bench analyse --model spc --engine gromacs
water-bench diagnostics            # the protocol-variation runs
water-bench report                 # results/summary.{md,csv} and figures
```

Runs live under `/ssd1_nas_md/water_benchmarking/<model>/<engine>/`.
