# Protocol diagnostics (SPC, 298.15 K, 1 atm)

Each run changes one thing relative to the benchmark protocol; all are analysed with
the same code as the main runs. Units: rho kg m^-3, dH_vap kJ mol^-1, D 1e-9 m^2 s^-1
(Yeh-Hummer corrected), tau2 ps. y is the dimensionless box-dipole fluctuation;
eps = 1 + y is the conducting-boundary relation, eps Neumann(61) the finite-eps_rf one
(see README: the latter sits near its pole here and is not the number to quote).

| run                                             | changed from protocol                                        | R_c (nm) | electrostatics  | N     | rho   | dH_vap | D    | tau2(HH) | y    | eps=1+y | eps Neumann(61) |
|-------------------------------------------------|--------------------------------------------------------------|----------|-----------------|-------|-------|--------|------|----------|------|---------|-----------------|
| thermostat v-rescale, barostat C-rescale (1 ns) | Berendsen T/P coupling -> v-rescale / C-rescale              | 1.8      | RF eps_rf = 61  | 2048  | 976.1 | 44.19  | 4.34 | 1.17     | 65.9 | 66.9    | 143             |
| cutoff 1.4 nm (1 ns)                            | R_c 1.8 -> 1.4 nm (R_c/L 0.45 -> 0.35)                       | 1.4      | RF eps_rf = 61  | 2048  | 974.6 | 44.16  | 4.25 | 1.17     | 59.3 | 60.3    | 116             |
| cutoff 0.9 nm (1 ns)                            | R_c 1.8 -> 0.9 nm (R_c/L 0.45 -> 0.23)                       | 0.9      | RF eps_rf = 61  | 2048  | 964.9 | 44.17  | 4.02 | 1.26     | 68.1 | 69.1    | 153             |
| 16384 waters (8x box) (1 ns)                    | 2048 -> 16384 waters, L 3.98 -> 7.96 nm (R_c/L 0.45 -> 0.23) | 1.8      | RF eps_rf = 61  | 16384 | 975.7 | 44.19  | 4.36 | 1.17     | 75.6 | 76.6    | 197             |
| RF eps_rf = infinity (1 ns)                     | eps_rf 61 -> infinity (conducting boundary; k_rf +2.4%)      | 1.8      | RF eps_rf = inf | 2048  | 976.0 | 44.19  | 4.32 | 1.17     | 61.9 | 62.9    | 126             |
| PME (1 ns)                                      | reaction field -> particle-mesh Ewald                        | 1.8      | PME             | 2048  | 975.7 | 44.18  | 4.37 | 1.17     | 75.6 | 76.6    | 197             |
| experiment                                      | -                                                            | -        | -               | -     | 997.0 | 43.99  | 2.30 | 2.00     | -    | 78.4    | -               |
