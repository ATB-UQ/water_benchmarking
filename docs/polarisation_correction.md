# Why a rigid water model's ΔH_vap needs a polarisation correction

The empirical case first, then the first-principles derivation behind it, then its limits.
Notation: μ_g is the gas-phase dipole of the real molecule (1.855 D), μ is the fixed dipole of
the model, α the molecular polarisability (1.44 Å³), E_pol ≡ (μ − μ_g)²/2α.

## 1. The case, in one table: fitting to the *uncorrected* ΔH_vap is the inconsistent choice

A rigid model's raw ΔH_vap must exceed experiment by the energy of distorting a molecule
into its liquid-phase state — equation (6) of the derivation in §2 below, which puts the
electronic part of that energy at E_pol = (μ − μ_g)²/2α with no free parameter. The
consequence is visible in the models themselves. Suppose a model is instead fitted so
that its raw ΔH_vap equals experiment. Then U_inter,model = U_inter +
E_dist, i.e. the model's liquid is *under-bound* by E_dist ≈ 5–7 kJ mol⁻¹ (12–16 % of the
cohesive energy). The only way to under-bind with a Coulomb-plus-Lennard-Jones form is a
smaller dipole. That is the trajectory of SPC (μ = 2.27 D) and TIP3P (2.35 D with a
TIP3P geometry), and their well-known consequences follow: too little cohesion, too fast
diffusion (SPC nearly 2× experiment in this benchmark, TIP3P 2.4×), too weak orientational
correlation, ε 15 % low. Berendsen, Grigera & Straatsma (*J. Phys. Chem.* 1987, **91**,
6269) introduced E_pol precisely to escape this: SPC/E keeps the SPC geometry and raises
the dipole to 2.35 D, fits ΔH_vap − E_pol to experiment, and gains density, D and ε at a
stroke. The modern models that score best on the Vega–Abascal benchmark — TIP4P/2005,
fitted with the Berendsen term included (Abascal & Vega, *J. Chem. Phys.* 2005, **123**,
234505; Vega & Abascal, *Phys. Chem. Chem. Phys.* 2011, **13**, 19663) — and the
optimised 3- and 4-point models OPC3/OPC and TIP3P-FB/TIP4P-FB all sit at μ ≈ 2.4–2.5 D,
i.e. they implicitly or explicitly accept (6). The rest of this note derives (6) and (9),
and then says what the derivation does *not* establish.

This benchmark shows the same thing directly. Under one protocol:

| model | μ (D) | E_pol | raw ΔH_vap | corrected | experiment |
|---|---|---|---|---|---|
| SPC | 2.27 | (3.76) | 44.2 | — | 43.99 |
| SPC/E | 2.35 | 5.24 | 49.27 | 44.03 | 43.99 |
| OPC3 | 2.43 | 7.03 | 51.66 | 44.63 | 43.99 |

The two enhanced-dipole models are 12 % and 17 % over experiment raw and within 0.1 % and
1.4 % corrected, *while also* reproducing density, D and (for OPC3) ε. SPC matches the raw
ΔH_vap and gets the rest wrong. A single physically motivated constant per model, with no
free parameter, reconciles the three — which is the behaviour one expects of a genuine
term rather than a fudge.

## 2. The supporting derivation

Three steps: what ΔH_vap is (§2.1), what a rigid model computes instead (§2.2), and why
the difference has the closed form (μ − μ_g)²/2α (§2.3).

### 2.1 What ΔH_vap actually measures

Per mole, ΔH_vap = H_gas − H_liq. The vapour at 1 atm is close to ideal, so
H_gas = U_gas + RT, and the pV term of the liquid is negligible (0.002 kJ mol⁻¹), so
H_liq ≈ U_liq. Take the energy of an isolated, relaxed molecule as zero:

    ΔH_vap = RT − U_liq                                                    (1)

U_liq is the total potential energy of the liquid per molecule, and for real water it has
two parts. Molecules in the liquid are not the molecules of the gas: their electron
clouds are polarised by the field of their neighbours (μ rises from 1.85 D to roughly
2.5–3 D), and their geometry is slightly stretched. Call the energy stored in that
distortion, per molecule, E_dist > 0, and the interaction energy between the *already
distorted* molecules U_inter < 0:

    U_liq,real = U_inter + E_dist                                          (2)

so that

    ΔH_vap,exp = RT − U_inter − E_dist                                     (3)

On vaporisation the real molecule relaxes; E_dist is *returned*. That is why it appears
with a minus sign in (3): the liquid is less bound than U_inter alone says, by the cost
of distorting its molecules.

### 2.2 What a rigid non-polarisable model computes

A rigid fixed-charge model has no intramolecular degrees of freedom. Its molecules carry
the same enhanced dipole μ in the liquid and in the gas, and its total energy is purely
intermolecular:

    U_liq,model = U_inter,model                                            (4)
    ΔH_vap,model = RT − U_inter,model                                      (5)

The model is built so that its molecules interact like the polarised molecules of the
real liquid — that is the whole point of the enhanced dipole. If it succeeds,
U_inter,model ≈ U_inter, and comparing (5) with (3):

    ΔH_vap,model − ΔH_vap,exp = E_dist                                     (6)

The raw model value must exceed experiment by the distortion energy, **even for a model
whose liquid is perfect**. The model's vapour is a gas of pre-polarised molecules that
cannot relax, which is a different and higher-energy state than real vapour. The
correction does not repair a defect of the liquid; it accounts for the fact that the
model's *gas* is unphysical — which is why, in the table of §1, the models that fit the raw
value are the ones whose liquids are wrong. Subtracting E_dist compares like with like: the
liquid-phase energy of polarised molecules, model against real.

### 2.3 Why E_dist = (μ − μ_g)²/2α: the linear-response half

This is not a fitted fudge. It follows from linear response, and the factor of two is the
familiar "half" of induction energetics.

Consider a molecule of polarisability α in a local field E from its neighbours. Its
induced dipole is μ_ind = αE, and the work done against the molecule's own restoring
force to build that dipole is

    W_dist = ∫₀^{μ_ind} (m/α) dm = μ_ind²/2α = ½αE²                       (7)

The interaction of that induced dipole with the field is −μ_ind·E = −αE². The net
induction contribution to the energy of a *polarisable* molecule at self-consistency is
therefore

    −αE² + ½αE² = −½αE²                                                     (8)

Now the fixed-charge model. It represents the polarised molecule by a permanent dipole
μ = μ_g + Δμ with Δμ playing the part of ⟨μ_ind⟩. Its pair energy contains the full
interaction −Δμ·E, but no term (7): nothing was spent to make Δμ, because it was there
from the start. Per molecule the model therefore over-binds, relative to a polarisable
description at the same mean field, by exactly

    E_pol = W_dist = Δμ²/2α = (μ − μ_g)²/2α                                (9)

This identifies the electronic part of E_dist in (6) with the model's own Δμ. It is the
same result as the statement that a linear dielectric in a field has free energy −½αE²
while a fixed dipole of the same size has −αE²; the "missing half" is the self-energy.

Two features of (9) are worth noting:

- **It uses the model's own μ, not the real liquid dipole.** The model's Δμ is the
  induction *it* has pre-paid, so (9) is its own self-energy. Inserting the real liquid
  dipole (~2.9 D) would give ~20 kJ mol⁻¹, which is a statement about real water, not
  about the model, and is not the quantity in (6). (That the model's μ ≈ 2.35–2.43 D is
  below the real ~2.9 D is itself understood: a non-polarisable model's charges are
  *effective* charges, screened by the electronic dielectric ε_el ≈ 1.78 of the medium
  they implicitly sit in — Leontyev & Stuchebrukhov, *J. Chem. Phys.* 2009, **130**,
  085102; *Phys. Chem. Chem. Phys.* 2011, **13**, 2613 — and 2.9/√1.78 ≈ 2.2 D. The
  correction is consistent within that effective picture.)
- **It leaves everything except the energy bookkeeping untouched.** Forces at the mean
  field are identical whether the dipole is induced or fixed; (9) is a constant per
  molecule. Structure, density, diffusion, rotational relaxation, the dielectric constant,
  and every free-energy difference in which the number of water molecules is conserved
  (hydration free energies, binding, conformational equilibria) are unaffected. Only
  quantities that compare the liquid with the model's gas — ΔH_vap, vapour pressure, the
  coexistence curve — carry the term.

## 3. The limits — what the argument does not establish

- **Linear response, mean field.** (9) treats Δμ as a fixed mean induced dipole. Real
  induction fluctuates with the instantaneous field, and the variance contributes to the
  energy; the correction captures the mean only. It is also isotropic scalar α, which is a
  good approximation for water (α is nearly isotropic, 1.42–1.47 Å³) but not general.
- **The inputs are external.** μ_g is precise, but the appropriate α is the molecular
  polarisability *in the liquid*, which differs from the gas-phase value by a few per cent,
  and reasonable choices move E_pol by ±0.2–0.3 kJ mol⁻¹. The corrected number inherits
  that uncertainty; the raw number is an exact property of the model.
- **Electronic only.** E_dist in (2) also includes the intramolecular strain of the
  stretched liquid-phase geometry (r_OH longer by ~0.01 Å, angle opened by ~1°). It is
  small — of order 1 kJ mol⁻¹ or less — but it is not in (9) and rigid models cannot
  represent it.
- **Nuclear quantum effects are absorbed elsewhere.** Real ΔH_vap includes the change in
  zero-point and quantised vibrational energy between liquid and gas, which a classical
  model cannot compute and fits implicitly. Path-integral studies (e.g. Habershon,
  Markland & Manolopoulos, *J. Chem. Phys.* 2009, **131**, 024501) put the classical–
  quantum difference at a few kJ mol⁻¹ — comparable to E_pol. A classical model fitted to
  experiment is therefore always an effective model, and E_pol is not the only implicit
  term; it is the only one with a closed form.
- **Convention, not law, decides whether to apply it.** Because a model's other
  parameters were tuned against *some* ΔH_vap, the comparison that is meaningful for a
  given model is the one its authors made. SPC, TIP3P, GROMOS 53A6/54A7 (Oostenbrink
  *et al.*, *J. Comput. Chem.* 2004, **25**, 1656), OPLS and the organic-liquid benchmarks
  built on them (Caleman *et al.*, *J. Chem. Theory Comput.* 2012, **8**, 61) used the raw
  value; SPC/E, TIP4P/2005 and OPC3 (Izadi & Onufriev, *J. Chem. Phys.* 2016, **145**,
  074501, whose SPC/E entry of 10.43 kcal mol⁻¹ is the corrected one) used the corrected
  value. Applying either to a model built on the other breaks its own validation. This
  benchmark follows each model's convention, and reports both numbers.

## 4. Is the argument specific to water?

No. Equations (1)–(9) hold for any molecule whose condensed-phase charge distribution is
enhanced over its gas-phase one, and every fixed-charge force field pre-polarises to some
degree. What is specific to water is the size: E_pol ∝ Δμ²/α, and water pairs one of the
largest relative enhancements of any small molecule (Δμ ≈ 0.5–0.6 D on a 1.85 D base)
with one of the smallest polarisabilities. Typical organic liquids have smaller Δμ and
α several times larger, so E_pol falls below a kJ mol⁻¹ and is customarily neglected; small
strongly hydrogen-bonded molecules (methanol, formamide) are the intermediate cases where
it reaches a few kJ mol⁻¹ and the neglect is a choice. For a solute in a protein force
field the issue is moot for ΔH_vap, which is never measured, but the same physics
reappears as the effective-charge (electronic-continuum) scaling of ions and polar groups,
which is the same "missing ε_el" seen from the solute's side.
