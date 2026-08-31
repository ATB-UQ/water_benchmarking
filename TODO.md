# TODO — SPC / SPC/E water benchmark (GROMOS on Gadi, GROMACS-GPU on Setonix)

> **Status 2026-08-27.** Everything up to submission is built, tested and verified locally:
> package scaffold, sibling registration, both topologies, the 2048-water box, the full .imd and
> .mdp ladders, the Gadi driver, all four analyses (18 tests pass), and the report. Step 0 is
> complete and changed three things — see README "What Step 0 established". A local `emin` on the
> real box finished successfully and relieved the seam contacts (0.2345 -> 0.2534 nm).
>
> **Update, 2026-08-28.** Done. Both engines complete and analysed for both models (GROMACS on
> Setonix GPU, GROMOS on Gadi as ten concurrent 1 ns replicates); `results/summary.md`,
> `results/diagnostics.md` and the figures are final and pushed. The dielectric constant is read
> through the boundary each engine realises (GROMOS: Neumann at ε_RF = 61; GROMACS: conducting,
> confirmed on GPU and CPU kernels alike) — see README. Open question, not pursued here: why
> GROMACS's reaction field does not realise the finite-ε_RF boundary GROMOS's does (next test:
> GROMOS at ε_RF → ∞). Still to do in the superproject: commit the `siblings.txt`/`.gitignore` rows.

> **Update, 2026-08-28 (OPC3).** A third model is registered: **OPC3** (Izadi & Onufriev 2016,
> *J. Chem. Phys.* 145:074501), **GROMACS only** — 54A7 has no building block for it, and the
> two-engine agreement is already established. Everything up to submission is built and verified
> locally: the model registry, the packaged `data/opc3.itp` (amber19sb's parameters restated in
> gromos54a7's C6/C12 convention), per-engine build gating, and 20 new tests including a real
> `grompp` + `gmx dump` round trip that reads charges, C6/C12 and the SETTLE geometry back out of
> the .tpr. SPC and SPC/E inputs regenerate byte-identically. Remaining:
>
> **Bug found while adding OPC3, fixed:** `gromacs.stages()` inherited
> `generate_velocities=True` on all ten production segments from the shared
> `imd.run_ladder()`, which is written for the GROMOS shape (ten independent
> replicates). GROMACS *chains* its segments, so this regenerated velocities at the
> head of every one, at the same seed each time — neither a chain nor a set of
> replicates. The `.mdp` files actually run on Setonix chain correctly, so the
> generator had drifted from them; with the fix, SPC regenerates byte-identical to
> all 14 published files and SPC/E to 13 of 14. The one remaining difference is
> `spce/eq1.mdp`'s `gen_seed` (published 770001, code 770100): the published run
> predates the 100-spacing that keeps replicate seeds from colliding across models.
> Statistically equivalent, not worth re-running, but it means the SPC/E
> equilibration cannot be reproduced bit-for-bit from the current code.
>
> - [x] built, submitted (Setonix 47753597, 26m50s), collected, analysed — every property with a
>       published value inside its literature range; results in the README table
> - [x] protocol sweep at 1.4 nm / 2 fs / RF 78.4, RF 61, PME (Setonix 47758274) —
>       `results/opc3_settings.md`; recommendation: 1.4 nm, 2 fs, RF, 10 ns, 0.2 ps sampling
> - [ ] `water-bench report` — regenerate `results/summary.md` with the `OPC3/gromacs` column
>       (needs the 1.8 nm analysis re-run through `report`, ~1 h of streaming)
> - [x] `LITERATURE["opc3"]` filled in from Table III of the paper (ρ 0.996 ± 0.001, ε 78.4 ± 1,
>       D 2.30 ± 0.02, ΔH_vap 10.73 ± 0.004 kcal/mol = 44.89 kJ/mol). The paper reports no
>       rotational correlation time, so `tau2_HH` stays absent and prints `[?]` — unchecked, which
>       is the honest reading, rather than silently looking confirmed.
> - [x] README results table filled in for OPC3

> **Update, 2026-08-31 (data retention).** The raw trajectories are gone: 35.7 GB
> locally and 13 GB on Setonix `/scratch`, all 69 files SHA-256'd first into
> `audit/`. What remains is 1.2 GB of inputs, energies, logs and final structures.
> Density and dH_vap stay recomputable from the kept `.edr`/`.tre.gz` (verified:
> OPC3 reproduces 994.4 / 51.66 exactly); D, tau2 and eps do not, and rest on the
> audit records plus the checksums. `audit/README.md` states the asymmetry, and
> `water_benchmarking.audit.build()` regenerates a record for any future run.

Ordered checklist. Full rationale in `~/.claude/plans/benchmark-two-classical-md-whimsical-dahl.md`.
Fixed decisions: **N = 2048** (a ≈ 3.95 nm; 1024 violates minimum image at 1.8 nm), **10 × 1 ns
production per model/engine**, water as GROMOS **solvent**, 298.15 K / 1 atm, ε_RF = 61 for both
models (the peptide-protocol value; SPC/E literature would use ~71 — document, don't change).

Reference `.imd` to copy settings from (peptide protocol):
`/ssd1_nas_md/protein_validation/runs/peptides_v3/ATB_protein/gb1/md_rep_1_run_1_to_10_rep_1/md_rep_1_run_1_to_10_rep_1_run_1.imd`
Cross-engine settings audit: `/ssd1_nas_md/protein_validation/settings.md`.

---

## 0. Repo scaffold

- [x] `git init`; `pyproject.toml` (package `water_benchmarking`, src layout, console script
      `water-bench`), `README.md`, `.gitignore` (`runs/`, `results/*.trc*`, `*.tre*`).
- [x] Register as a sibling: row in `/home/atb/ATB/siblings.txt`
      (`water_benchmarking https://github.com/ATB-UQ/water_benchmarking main editable`) + entry in the
      root `.gitignore` sibling block; `bash /home/atb/ATB/scripts/sync_siblings.sh --check`.
- [x] `uv pip install --python /home/atb/ATB/.venv/bin/python -e /home/atb/ATB/water_benchmarking --no-deps`
- [ ] Create the GitHub remote under ATB-UQ (manual; local repo + commit exist).
- [x] Run root: `mkdir -p /ssd1_nas_md/water_benchmarking/{spc,spce}/{gromos,gromacs}`.

## 1. Step-0 feasibility checks (local, `/opt/gromos/1.6.0/bin`) — do these first

- [x] Pure-solvent topology:
      `make_top @build 54A7.mtb @param 54A7.ifp @seq "" @solv H2O > spc.top` (mtb/ifp in
      `/home/atb/ATB/gromos_job_wrapper/src/gromos_job_wrapper/lib/`). Repeat with `@solv H2OE` → `spce.top`.
      If an empty `@seq` is refused: write the topology from Python (solute blocks zero-count,
      `SOLVENTATOM`/`SOLVENTCONSTR` copied from the mtb solvent block).
- [x] md++ with `SYSTEM NPM 0 NSM 2048`: 10 serial steps on the cut box (`md @topo @conf @input`).
      **Fallback** if NPM=0 is rejected: 1 solute `H2O`/`H2OE` block + 2047 solvent — record as a deviation.
- [x] `check_box` on the cut box; `ene_ana` parses the 10-step `.tre` with
      `gromos_job_wrapper/.../lib/ene_ana.md++.lib` (tre-version gate).
- [x] Confirm gromos++ `epsilon` sums solvent dipoles (run it on the 10-step `.trc`; a non-zero
      `<M^2>` proves it). If solute-only, the Python implementation (§5) is the only ε path.

## 2. Box (`src/water_benchmarking/box.py`)

- [x] Source: `gromos_job_wrapper/src/gromos_job_wrapper/lib/H2O_box.g96` (5384 SPC, cubic 5.4937 nm).
- [x] Per-molecule COM; rank by max(|x|,|y|,|z|) about box centre; keep the innermost 2048;
      edge `a = (2048 / ρ_N,source)^(1/3)` ≈ 3.98 nm; wrap; write `water_2048.cnf` with cubic `GENBOX 1`.
- [x] Asserts: exactly 2048 molecules, `a/2 > 1.8` nm, min O–O across all images > 0.23 nm.
- [x] Same coordinates for SPC and SPC/E. Also write `water_2048.gro` (six decimals, GROMACS names
      `OW HW1 HW2`, resname `SOL`) for §6.
- [x] Test: `tests/test_box.py`.

## 3. GROMOS inputs (`imd.py`, `protocol.py`)

- [x] Template = peptide `.imd` minus solute blocks. Per block:
      `STEP DT 0.001`; `BOUNDCOND NTB 1 NDFMIN 6`; `SYSTEM NPM 0 NSM 2048`; `FORCE` all 1, `NEGR 1 NRE 6144`;
      `CONSTRAINT NTC 2 / SHAKE 1e-5 / solvent SHAKE 1e-5`; `PAIRLIST 1 5 1.8 1.8 0.4 0`;
      `NONBONDED NLRELE 1 APPAK 0 RCRF 1.8 EPSRF 61 NSLFEXCL 1 / NSHAPE -1 ASHAPE 1.8 …` (rest verbatim);
      `MULTIBATH` weak coupling, **1 bath** 298.15/0.1, DOFSET 1 → last 6144, bath 1;
      `PRESSURESCALE 2 1 4.575e-4 0.5 2`, PRES0 0.06102 diag; `COMTRANSROT 1000`; `PRINTOUT 500 0`.
      No `POSITIONRES`, no `ROTTRANS`.
- [x] Stages: `emin` (`ENERGYMIN NTEM 1 NCYC 1 DELE 0.01 DX0 0.01 DXM 0.05 NMIN 100`, 1000 steps,
      `PAIRLIST algorithm 0`) → `eq1` NVT 50 K, 10 k steps, `INITIALISE NTIVEL 1 TEMPI 50 IG 770000+{0|1}`
      → `eq2` NVT 298.15 K, 20 k → `eq3` NPT 298.15 K, 100 k (100 ps) → `md_01..md_10` NPT, 1 000 000 steps,
      `WRITETRAJ NTWX 100 NTWE 100` (eq: 1000/1000).
- [x] Test: rendered blocks parse (round-trip a keyword grep), NSTLIM/NTWX per stage.

## 4. GROMOS on Gadi (`gadi.py`, `cli.py run`)

- [x] One segment = one call of `/home/atb/ATB/gromos_job_wrapper/deployment/gadi_md.sh`
      with `@topo @conf @input @fin @trc @tre` (absolute paths under `/ssd1_nas_md/water_benchmarking/<model>/gromos/`).
      Shim is a drop-in `md`: stages inputs, renders PBS (`-P m72 -q normal -l storage=scratch/m72`,
      `mem = 4 GB × ncpus`, walltime = `NSTLIM × GJW_GADI_SPS + 600 s`), `qsub -W block=true`, pulls back
      with size verification. Never edit the shim in place (rename-replace only).
- [x] Env: `GJW_GADI_NCORES=8` for emin/eq, `24` for production; `GJW_JOB_NAME=w_{spc|spce}_{stage}`
      (≤15 chars); `GJW_GADI_SPS`: start at `0.02`, then **set from eq3's `Wall time simulation`**
      (expect ~0.005–0.01 s/step at 24 ranks; the default 0.08 would request 22 h).
- [x] Chain: `emin → eq1 → eq2 → eq3 → md_01 … md_10`, each `@conf` = previous `@fin`; gzip `.trc/.tre`
      after each pull; abort the chain unless the log ends `finished successfully`
      (shape: `runs/peptides_v3/.../md_rep_1_run_1_to_10_rep_1/*_local.sh`, `runs/peptides_v3/gromos_rep.sh`).
- [ ] Launch the two models as two independent `nohup` chains: `water-bench run --engine gromos --model spc`
      and `… --model spce`. Optional non-blocking mode: `GJW_SUBMIT_ONLY=1` + `gadi_md.sh --collect <log>.job.json`.
- [ ] After eq3: density within ~3 % of 997 kg m⁻³ before launching production.

## 5. Analysis (`trc.py`, `analysis/`)

- [x] `trc.py`: streaming reader for (gzipped) GROMOS `.trc` — `TIMESTEP`, `POSITIONRED`, `GENBOX`;
      yields `(t, box, xyz[6144,3])`. Test on a hand-written 2-frame file.
- [x] `density_hov.py`: `ene_ana @prop densit totpot pressu boxvol` over all 10 `.tre.gz`;
      ΔH_vap = −⟨U_pot⟩/N + RT (rigid model ⇒ no gas-phase leg); for SPC/E also report the
      self-polarisation-corrected value (−5.22 kJ mol⁻¹). Errors: block averaging + KS equilibration
      (`ks_convergence_analysis`, `block_averaging` siblings; see
      `gromos_job_wrapper/.../helpers/data_processing.py:361 ks_error_analysis`).
- [x] `diffusion.py`: MSD of OW, unwrapped via nearest-image displacement between consecutive frames;
      Einstein slope over 10–100 ps lags; report D_PBC and Yeh–Hummer D∞ = D_PBC + 2.837 k_B T/(6π η L),
      η = 0.89 mPa s, correction shown separately. Cross-check: `diffus @topo @pbc r @time 0 0.1 @dim x y z @atoms s:OW @traj …`.
- [x] `rotation.py`: per-molecule unit vectors OH (both), HH, dipole (bisector); C_l(t) = ⟨P_l(u(0)·u(t))⟩,
      l = 1, 2, out to 20 ps, multiple time origins; τ_l by integral and by single-exponential fit of 1–10 ps.
- [x] `dielectric.py`: M = Σ q r per frame (molecules whole under SHAKE; no gathering);
      (ε−1)(2ε_RF+1)/(2ε_RF+ε) = (⟨M²⟩−⟨M⟩²)/(3 ε₀ V k_B T), ε_RF = 61; running estimate vs time.
      Cross-check: `epsilon @topo @pbc r @temp 298.15 @e_rf 61 @traj …` (only if §1 showed it includes solvent).
- [x] `experiment.py` (298.15 K, 1 atm, cite in code): ρ 997.05 kg m⁻³; ΔH_vap 43.99 kJ mol⁻¹;
      D 2.30×10⁻⁹ m² s⁻¹; ε 78.4; τ₂(HH) ≈ 2.0 ps, τ₂(OH) ≈ 1.95 ps (NMR); τ₁(dipole) ≈ 8.3 ps (Debye).
      Literature model values for sanity: SPC ρ 972–985, ε ≈ 65, D 3.9–4.3; SPC/E ρ ≈ 998, ε ≈ 71, D 2.5–2.8.
- [x] Synthetic tests: fixed vectors → C₂ = 1; random walk → known D; fixed M → known ε.

## 6. GROMACS-GPU mirror on Setonix (`gromacs.py`)

- [ ] **First**: `/home/atb/ATB/gromacs_pipeline` working tree is broken (`runners.py`: `BatchRunner`
      deleted but still referenced at line ~234; `LocalRunner.submit` calls removed `_ssh/_script`).
      Either finish that refactor or `git stash` and work from HEAD before importing anything.
- [x] Topology: minimal `.top` including `gromos54a7.ff/forcefield.itp` + `spc.itp` | `spce.itp`
      (`/home/atb/opt/gromacs/2026.1/share/gromacs/top/gromos54a7.ff/`), `[ molecules ] SOL 2048`;
      SETTLE. Coordinates: `water_2048.gro` from §2. `gmx grompp -maxwarn 0`.
- [x] `.mdp` (mirror of §3; pipeline template `gromacs_pipeline/src/gromacs_pipeline/mdp.py`):
      `dt 0.001`, `cutoff-scheme Verlet`, `rlist = rcoulomb = rvdw = 1.8`, `coulombtype reaction-field`,
      `epsilon_r 1`, `epsilon_rf 61`, `vdwtype cut-off`, **`vdw-modifier none`** (default Potential-shift
      shifts the reported LJ energy → wrong ΔH_vap), `DispCorr no`, `tcoupl berendsen tau_t 0.1 ref_t 298.15`,
      `pcoupl berendsen isotropic tau_p 0.5 compressibility 2.755e-5 ref_p 1.01325`,
      `comm-mode linear nstcomm 1000`, `nstxout-compressed = nstenergy = 100` (eq: 1000),
      `gen_vel yes gen_temp 50 gen_seed 770000+{0|1}` in eq1 only, `continuation yes` after.
      Stages: `em` (steep, 1000) → `eq1` NVT 50 K 10 k → `eq2` NVT 298.15 K 20 k → `eq3` NPT 100 k →
      `md_01..md_10` 1 000 000 steps each, chained with `-cpi`.
- [x] Host `setonix-gpu` from `gromacs_pipeline/src/gromacs_pipeline/config.py` (account m72; one GCD per
      run — a whole node is slower at this size; run spc and spce concurrently on separate GCDs).
      `mdrun -nb gpu -bonded gpu` **only** — `-pme gpu` is fatal under reaction field;
      `OMP_NUM_THREADS` must equal `-ntomp`.
- [ ] Pull back `.xtc`/`.edr`/`.log` to `/ssd1_nas_md/water_benchmarking/<model>/gromacs/`.
- [ ] Analysis on the same code path: `gmx trjconv -f seg.xtc -s seg.tpr -o seg.g96` → `trc.py` reads
      it (POSITIONRED + GENBOX). Density/ΔH_vap from `gmx energy` (Density, Potential; pipeline
      `analysis.py energy_summary`). Cross-checks: `gmx msd`, `gmx rotacf -P 2`, `gmx dipoles -epsRF 61`.

## 7. Report (`report.py`, `water-bench report`)

- [x] `results/summary.csv` + `summary.md`: rows = property; columns = SPC/GROMOS, SPC/GROMACS,
      SPC-E/GROMOS, SPC-E/GROMACS, experiment, % deviation; every value with a block-average error.
- [ ] Plots: ρ(t), MSD(t), C₂(t) per vector, ε running average, per model × engine.
- [ ] Flag engine–engine differences beyond the error bars (expected signatures: Verlet buffer vs
      exact cutoff, SHAKE vs SETTLE, RF self-term handling).
- [ ] README: protocol table, deviations from the peptide protocol (1 bath, no solute, ε_RF 61 for SPC/E),
      how to rerun. Commit tables + plots; trajectories stay on `/ssd1_nas_md`.

## Done when

- [ ] §1 checks all pass and are recorded in README.
- [ ] `pytest tests/` green.
- [ ] Both GROMOS chains and both GROMACS chains finished 10/10 segments with `finished successfully` /
      clean `.log`.
- [ ] ε running average flat over the last 3 ns for all four runs; Python D and ε agree with
      `diffus`/`epsilon`/`gmx` cross-checks on one shared segment.
- [ ] Summary table filled, within literature ranges for both models.