# Water model benchmark

2048 water molecules, 298.15 K, 1 atm, 10 x 1 ns production.
Single-range cutoff 1.8 nm, reaction field eps_rf = 61, timestep 1 fs.

## Results

| Property                 | Unit          | OPC3/gromacs     | SPC/gromos       | SPC/gromacs     | SPCE/gromos      | SPCE/gromacs     | Experiment | Source                                                    |
|--------------------------|---------------|------------------|------------------|-----------------|------------------|------------------|------------|-----------------------------------------------------------|
| Density                  | kg m^-3       | 994.4 +/- 0.072  | 974.8 +/- 0.071  | 976.0 +/- 0.078 | 996.4 +/- 0.086  | 997.4 +/- 0.069  | 997.0      | Kell 1975, J. Chem. Eng. Data 20:97                       |
| dH_vap                   | kJ mol^-1     | 51.66 +/- 0.0016 | 44.16 +/- 0.0013 | 44.19 +/- 0.001 | 49.24 +/- 0.0016 | 49.27 +/- 0.0014 | 43.99      | Wagner & Pruss 2002, IAPWS-95                             |
| dH_vap - E_pol           | kJ mol^-1     | 44.63            | -                | -               | 44.02            | 44.05            | 43.99      | Wagner & Pruss 2002, IAPWS-95                             |
| Self-diffusion D         | 1e-9 m^2 s^-1 | 2.37 +/- 0.004   | 4.45 +/- 0.0032  | 4.31 +/- 0.014  | 2.78 +/- 0.021   | 2.75 +/- 0.004   | 2.30       | Holz et al. 2000, PCCP 2:4740                             |
| Rot. corr. time tau2(HH) | ps            | 2.28             | 1.16             | 1.17            | 1.96             | 1.98             | 2.00       | NMR relaxation; Ludwig 2001, Angew. Chem. 40:1808         |
| Rot. corr. time tau1(mu) | ps            | 5.36             | 2.81             | 2.83            | 4.60             | 4.63             | 8.30       | Debye tau_D (collective); Ronne et al. 1997, JCP 107:5319 |
| Dielectric constant      | -             | 78.3 +/- 2.4     | 68.1 +/- 1.9     | 66.3 +/- 1.7    | 69.8 +/- 2.2     | 69.4 +/- 1.5     | 78.4       | Fernandez et al. 1997, J. Phys. Chem. Ref. Data 26:1125   |

## Deviation from experiment

`[!]` marks a value outside the published range for that model, which points
at the setup rather than at the model -- unless a note below says otherwise.
`[?]` marks a property with no published range on record for that model, so
the run is unchecked against the model rather than confirmed to reproduce it.
tau1(mu) is omitted: the experimental Debye time is a collective quantity and
not directly comparable to the single-molecule correlation time simulated.

| Property                 | OPC3/gromacs | SPC/gromos  | SPC/gromacs | SPCE/gromos | SPCE/gromacs |
|--------------------------|--------------|-------------|-------------|-------------|--------------|
| Density                  | -0.3%        | -2.2%       | -2.1%       | -0.1%       | +0.0%        |
| dH_vap                   | +17.4%       | +0.4%  [!]  | +0.5%  [!]  | +11.9%      | +12.0%       |
| dH_vap - E_pol           | +1.4%        | -           | -           | +0.1%       | +0.1%        |
| Self-diffusion D         | +3.0%        | +93.5%  [!] | +87.2%      | +20.8%      | +19.4%       |
| Rot. corr. time tau2(HH) | +13.9%  [?]  | -41.9%      | -41.5%      | -1.9%       | -0.8%        |
| Dielectric constant      | -0.1%        | -13.2%      | -15.5%      | -10.9%      | -11.5%       |

- SPC/gromos, hov: published SPC dH_vap (41-44) is mostly at ~0.9 nm cutoffs; the 1.8 nm cutoff here recovers more attractive energy, so a value at the top of the range is the protocol, not an error
- SPC/gromacs, hov: published SPC dH_vap (41-44) is mostly at ~0.9 nm cutoffs; the 1.8 nm cutoff here recovers more attractive energy, so a value at the top of the range is the protocol, not an error

## Dielectric constant: relation used per run

y = (<M^2> - <M>^2) / (3 eps0 V kT) is what is measured; the relation that
turns it into eps depends on the boundary the engine realises (see README).

- OPC3/gromacs: conducting boundary, eps = 1 + y  (y = 77.3; conducting 78.3, Neumann 209)
- SPC/gromos: Neumann, eps_rf = 61  (y = 43.4; conducting 44.4, Neumann 68)
- SPC/gromacs: conducting boundary, eps = 1 + y  (y = 65.3; conducting 66.3, Neumann 140)
- SPCE/gromos: Neumann, eps_rf = 61  (y = 44.1; conducting 45.1, Neumann 70)
- SPCE/gromacs: conducting boundary, eps = 1 + y  (y = 68.4; conducting 69.4, Neumann 155)

## Figures

![dielectric](dielectric.png)
![msd](msd.png)
![c2_HH](c2_HH.png)
![density](density.png)
