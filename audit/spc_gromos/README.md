# Data audit — SPC production, GROMOS

A self-contained record of the ten SPC production simulations reported in the benchmark, and
of the analysis applied to them, assembled so the published numbers can be checked without the
bulk trajectory data. Everything here is derived from the run directory
`/ssd1_nas_md/water_benchmarking/spc/gromos`, which is where the excluded files remain.

Scope: **SPC only, GROMOS engine, production stage only.** The SPC/E runs, the GROMACS runs and
the protocol diagnostics are reported in the top-level `results/`, not here.

## Provenance

| | |
|---|---|
| Engine | GROMOS md++ 1.6.0 (MPI), build Mon Jul 13 16:53:48 AEST 2026 |
| Binary | `/home/564/mzs564/gromos_mpi/bin/md_mpi` (built per `gromos_job_wrapper/deployment/build_gromosxx.sh`) |
| Machine | NCI Gadi, `normal` queue, project m72, 24 MPI ranks per replicate |
| Submitted | 2026-08-27T17:43:48 (all ten within ~10 minutes, run concurrently) |
| Driver | `water_benchmarking.gadi.submit_production` via `gromos_job_wrapper/deployment/gadi_md.sh` |
| Force field | GROMOS 54A7 (`54A7.mtb` / `54A7.ifp`), solvent building block `H2O` |

| replicate | PBS job            | node              | wall (s) | s/step   | ns/day | completed |
|-----------|--------------------|-------------------|----------|----------|--------|-----------|
| md_01     | 177616435.gadi-pbs | gadi-cpu-clx-1303 | 16319.1  | 0.016319 | 5.29   | True      |
| md_02     | 177616448.gadi-pbs | gadi-cpu-clx-2199 | 16007.3  | 0.016007 | 5.4    | True      |
| md_03     | 177616466.gadi-pbs | gadi-cpu-clx-0249 | 15933.9  | 0.015934 | 5.42   | True      |
| md_04     | 177616517.gadi-pbs | gadi-cpu-clx-0597 | 15764.7  | 0.015765 | 5.48   | True      |
| md_05     | 177616521.gadi-pbs | gadi-cpu-clx-0597 | 16144.2  | 0.016144 | 5.35   | True      |
| md_06     | 177616527.gadi-pbs | gadi-cpu-clx-2041 | 16559.3  | 0.016559 | 5.22   | True      |
| md_07     | 177616534.gadi-pbs | gadi-cpu-clx-1046 | 15342.9  | 0.015343 | 5.63   | True      |
| md_08     | 177616541.gadi-pbs | gadi-cpu-clx-2857 | 14904.8  | 0.014905 | 5.8    | True      |
| md_09     | 177616550.gadi-pbs | gadi-cpu-clx-0004 | 15006.3  | 0.015006 | 5.76   | True      |
| md_10     | 177616554.gadi-pbs | gadi-cpu-clx-2949 | 15598.4  | 0.015598 | 5.54   | True      |

Every replicate ended with md++'s own `MD++ finished successfully`; the aggregate throughput was
54.9 ns/day across the ten concurrent jobs.

## System and protocol

2048 SPC water molecules, cubic box, no solute (`NPM 0`, `NSM 2048`). The starting configuration
was cut from the equilibrated 5384-molecule library box and taken through
emin → 10 ps NVT at 50 K → 20 ps NVT at 298.15 K → 100 ps NPT before production.

| Setting | Value | Where |
|---|---|---|
| Timestep | 1 fs | `STEP` |
| Length per replicate | 1 000 000 steps = 1 ns | `STEP` |
| Temperature | 298.15 K, Berendsen weak coupling, τ_T 0.1 ps, one bath | `MULTIBATH` |
| Pressure | 1 atm (0.06102 kJ mol⁻¹ nm⁻³), Berendsen, τ_P 0.5 ps, isotropic, molecular virial | `PRESSURESCALE` |
| Cutoff | single-range 1.8 nm, grid pairlist, updated every 5 steps | `PAIRLIST` |
| Electrostatics | reaction field, ε_RF = 61, RCRF = ASHAPE = 1.8 nm | `NONBONDED` |
| Constraints | `NTC 1` (solvent only); rigid geometry from `SOLVENTCONSTR`, SHAKE tol 10⁻⁵ | `CONSTRAINT` |
| Output | coordinates and energies every 100 steps (0.1 ps) | `WRITETRAJ` |
| COM motion | removed every 1000 steps | `COMTRANSROT` |

The ten replicates are **identical but for the velocity seed** (`INITIALISE` IG = 770001…770010),
each drawing fresh Maxwell velocities at 298.15 K from the same equilibrated box. They are
therefore independent, and the spread across them is a real error estimate rather than a
correlated one. `diff inputs/md_01.imd inputs/md_07.imd` shows only the title and the seed.

## Results, per replicate

Each replicate analysed on its own, with the same code that produced the aggregate. The leading
10 % of each is discarded as re-thermalisation after the fresh velocities.

| replicate  | seed   | rho           | U/N             | ΔH_vap         | P         | D_pbc         | D_corr        | τ₂(HH)        | τ₂(OH)        | τ₁(µ)         | y            | ε          |
|------------|--------|---------------|-----------------|----------------|-----------|---------------|---------------|---------------|---------------|---------------|--------------|------------|
|            |        | kg m⁻³        | kJ mol⁻¹        | kJ mol⁻¹       | atm       | 10⁻⁹ m² s⁻¹   | 10⁻⁹ m² s⁻¹   | ps            | ps            | ps            |              |            |
| md_01      | 770001 | 974.99        | -41.682         | 44.161         | 0.3       | 4.3412        | 4.5163        | 1.1589        | 1.0585        | 2.8065        | 42.59        | 66.15      |
| md_02      | 770002 | 974.57        | -41.6775        | 44.156         | 2.1       | 4.2523        | 4.4274        | 1.1585        | 1.0577        | 2.8055        | 43.625       | 68.6       |
| md_03      | 770003 | 974.9         | -41.6815        | 44.16          | 0.0       | 4.3618        | 4.5369        | 1.165         | 1.0611        | 2.8057        | 41.651       | 63.98      |
| md_04      | 770004 | 974.69        | -41.6797        | 44.159         | 0.1       | 4.2059        | 4.381         | 1.1622        | 1.0599        | 2.8204        | 46.342       | 75.36      |
| md_05      | 770005 | 974.77        | -41.681         | 44.16          | 2.0       | 4.2975        | 4.4726        | 1.1619        | 1.0598        | 2.7977        | 39.715       | 59.65      |
| md_06      | 770006 | 975.0         | -41.6809        | 44.16          | 1.1       | 4.2722        | 4.4473        | 1.1634        | 1.0594        | 2.822         | 43.918       | 69.31      |
| md_07      | 770007 | 974.49        | -41.6805        | 44.159         | 3.3       | 4.2528        | 4.4278        | 1.1643        | 1.0593        | 2.8042        | 42.2         | 65.24      |
| md_08      | 770008 | 974.51        | -41.6753        | 44.154         | 0.2       | 4.1742        | 4.3493        | 1.1595        | 1.0618        | 2.7791        | 44.908       | 71.73      |
| md_09      | 770009 | 975.07        | -41.6772        | 44.156         | 1.5       | 4.1609        | 4.336         | 1.1566        | 1.0573        | 2.799         | 40.192       | 60.7       |
| md_10      | 770010 | 974.72        | -41.6799        | 44.159         | 1.0       | 4.4388        | 4.6139        | 1.1636        | 1.0598        | 2.8298        | 44.823       | 71.52      |
| mean ± SEM |        | 974.77 ± 0.07 | -41.680 ± 0.001 | 44.158 ± 0.001 | 1.2 ± 0.3 | 4.276 ± 0.028 | 4.451 ± 0.028 | 1.161 ± 0.001 | 1.059 ± 0.000 | 2.807 ± 0.005 | 43.00 ± 0.67 | 67.2 ± 1.6 |

`D_pbc` is the Einstein slope as simulated, fitted over 100–500 ps; `D_corr` adds the Yeh–Hummer
finite-size correction (+0.175 × 10⁻⁹ m² s⁻¹ for this box). `y` is the dimensionless box-dipole
fluctuation (⟨M²⟩−⟨M⟩²)/(3ε₀VkT); `ε` is `y` read through Neumann's relation at ε_RF = 61, which is
the relation appropriate to GROMOS's reaction field (see the top-level README).

Spread across the ten replicates: density 974.49–975.07 kg m⁻³, ΔH_vap
44.15–44.16 kJ mol⁻¹, D_corr 4.34–4.61 × 10⁻⁹ m² s⁻¹,
τ₂(HH) 1.16–1.17 ps, ε 59.6–75.4.

## What is included, and what is not

| | size | |
|---|---|---|
| `inputs/` | 285 KB | topology, starting configuration (gzipped), all 14 `.imd` files |
| `logs/` | 376 KB | md++ logs, trimmed |
| `provenance/` | 16 KB | PBS job records, run table, SHA-256 of every audited file |
| `results/` | 3 KB | per-replicate table, aggregate statistics |
| **excluded** | **11 GB** | trajectories (`md_*.trc.gz`, 1.07 GB each), energy trajectories (`md_*.tre.gz`, 3.6 MB each), final configurations (`md_*.cnf`, 1 MB each) |

The excluded files stay in the run directory and are listed with their SHA-256 sums in
`provenance/SHA256SUMS.txt` (37 files), so this audit is tied to specific data rather than
to files of the same name.

**The logs are trimmed, not filtered.** A 1 ns log is 45 MB, of which 99.9 % is 10 000 energy
blocks printed every 100 steps. Each log here keeps md++'s complete header — the build, the host,
and every parameter md++ actually parsed (topology: 2048 solvents and 0 solute; force field:
Coulomb reaction field, grid pairlist, rectangular boundary, molecular virial; SHAKE; pressure
scaling; integration) — then the **first** energy block, a marker line stating exactly how many
blocks were removed, the **last** energy block, the timing, and the completion marker. Total
460 MB → 376 KB.

## Verifying

```bash
# the excluded bulk data is unchanged since the audit was taken
cd /ssd1_nas_md/water_benchmarking/spc/gromos && sha256sum -c <path-to>/provenance/SHA256SUMS.txt

# the replicates differ only by seed
diff audit/spc_gromos/inputs/md_01.imd audit/spc_gromos/inputs/md_07.imd

# every replicate completed
grep -c 'finished successfully' audit/spc_gromos/logs/md_*.log

# re-derive the aggregate from the per-replicate table
python -c "import csv,statistics as s; r=list(csv.DictReader(open('audit/spc_gromos/results/per_replicate.csv'))); \
print(s.mean(float(x['density_kg_m3']) for x in r))"

# regenerate everything from the run directory (needs the excluded files)
water-bench analyse --model spc --engine gromos
```
