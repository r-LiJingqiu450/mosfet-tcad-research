### H1: reference quality (A-R, body B1)

| State | ψ̂ (P1 L2) | G_DD ref (S/cm) | G_DD P4/L2−1 | G_DD L2/L1−1 | G_SD ref (S/cm) | G_SD P4/L2−1 | G_SD L2/L1−1 | max \|forms\| L0/L1/L2 | max KCL | B0 bracket L0/L1/L2 | H1 |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---|---|
| AR-1 | 11.997 | 3.527579e-12 | -1.42e-03 | -1.42e-03 | 3.504148e-12 | -1.42e-03 | -1.42e-03 | 1.4e-10/4.6e-10/9.6e-10 | 9.6e-10 | ok/ok/ok | **pass** |
| AR-2 | 15.997 | 1.616219e-10 | -1.83e-03 | -1.81e-03 | 1.615985e-10 | -1.83e-03 | -1.82e-03 | 9.3e-11/1.4e-10/1.1e-10 | 1.4e-10 | ok/ok/ok | **pass** |
| AR-3 | 21.017 | 2.109804e-08 | -2.33e-03 | -2.30e-03 | 2.109802e-08 | -2.33e-03 | -2.30e-03 | 6.0e-11/1.1e-10/1.9e-11 | 1.1e-10 | ok/ok/ok | **pass** |
| AR-4 | 12.000 | 2.550194e-09 | -1.30e-04 | -1.17e-04 | 2.545113e-09 | -1.30e-04 | -1.17e-04 | 3.9e-11/3.6e-11/3.0e-11 | 3.9e-11 | ok/ok/ok | **pass** |
| AR-5 | 16.000 | 1.203374e-07 | -1.25e-04 | -1.19e-04 | 1.203323e-07 | -1.25e-04 | -1.19e-04 | 2.9e-11/1.7e-11/1.0e-10 | 1.0e-10 | ok/ok/ok | **pass** |
| AR-6 | 21.020 | 1.597378e-05 | -9.99e-05 | -1.18e-04 | 1.597377e-05 | -9.99e-05 | -1.18e-04 | 1.8e-11/2.0e-12/7.6e-11 | 7.6e-11 | ok/ok/ok | **pass** |

H1 summary: 6 of 6 pass (AR-1, AR-2, AR-3, AR-4, AR-5, AR-6); kill (fewer than 4) = False.

### Reference values and tolerances (section 7)

| State | G_DD tol (rel) | G_SD tol (rel) | u_disc G_DD / G_SD | u_round G_DD / G_SD | u_solver | observed order G_DD / G_SD | P4 refined vs base G_DD | ΔL_SD (nm) | ΔL tol (nm) | G_DB (S/cm) | G_DB/G_DD |
|---|---:|---:|---|---|---:|---|---:|---:|---:|---:|---:|
| AR-1 | 2.85e-03 | 2.85e-03 | 1.4e-03 / 1.4e-03 | 9.6e-10 / 1.4e-14 | 6.6e-14 | 2.22/2.23 | -2.78e-04 | 34.209 | 1.326 | 2.343e-14 | 6.64e-03 |
| AR-2 | 3.65e-03 | 3.65e-03 | 1.8e-03 / 1.8e-03 | 1.1e-10 / 7.1e-15 | 7.7e-14 | 2.23/2.23 | -3.28e-04 | 38.515 | 1.685 | 2.343e-14 | 1.45e-04 |
| AR-3 | 4.66e-03 | 4.66e-03 | 2.3e-03 / 2.3e-03 | 1.9e-11 / 4.2e-14 | 4.4e-14 | 2.22/2.22 | -3.89e-04 | 40.844 | 2.138 | 2.344e-14 | 1.11e-06 |
| AR-4 | 2.60e-04 | 2.61e-04 | 1.3e-04 / 1.3e-04 | 3.0e-11 / 5.2e-14 | 1.5e-13 | 1.98/1.98 | +5.59e-05 | 364.229 | 0.166 | 5.081e-12 | 1.99e-03 |
| AR-5 | 2.50e-04 | 2.50e-04 | 1.2e-04 / 1.2e-04 | 1.0e-10 / 5.1e-15 | 1.2e-13 | 1.96/1.96 | +5.17e-05 | 381.533 | 0.154 | 5.102e-12 | 4.24e-05 |
| AR-6 | 2.36e-04 | 2.36e-04 | 1.2e-04 / 1.2e-04 | 7.6e-11 / 5.2e-14 | 1.5e-13 | 1.95/1.95 | -6.13e-06 | 395.180 | 0.143 | 5.129e-12 | 3.21e-07 |

### H-B: body effect at level 2 (B1 vs B0)

| State | Δ_DD | band | in | Δ_SD | band | in | G_DB/G₀ | P4 Δ_DD | P4 Δ_SD |
|---|---:|---|---|---:|---|---|---:|---:|---:|
| AR-1 | +4.19e-03 | [1.2e-03, 1.2e-02] | yes | -2.48e-03 | [-7.2e-03, -7.8e-04] | yes | 6.67e-03 | +4.19e-03 | -2.48e-03 |
| AR-2 | +9.11e-05 | [2.6e-05, 2.6e-04] | yes | -5.39e-05 | [-1.6e-04, -1.7e-05] | yes | 1.45e-04 | +9.12e-05 | -5.40e-05 |
| AR-3 | +6.98e-07 | [2.0e-07, 2.0e-06] | yes | -4.13e-07 | [-1.2e-06, -1.3e-07] | yes | 1.11e-06 | +6.95e-07 | -4.14e-07 |
| AR-4 | +1.60e-03 | [2.4e-04, 6.1e-03] | yes | -3.95e-04 | [-1.1e-03, -1.0e-04] | yes | 2.00e-03 | +1.60e-03 | -3.95e-04 |
| AR-5 | +3.40e-05 | [5.0e-06, 1.3e-04] | yes | -8.42e-06 | [-2.3e-05, -2.2e-06] | yes | 4.24e-05 | +3.40e-05 | -8.42e-06 |
| AR-6 | +2.57e-07 | [3.8e-08, 9.7e-07] | yes | -6.39e-08 | [-1.7e-07, -1.7e-08] | yes | 3.21e-07 | +2.57e-07 | -6.39e-08 |

H-B summary: Δ_DD in band 6/6, Δ_SD in band 6/6; supported (≥ 5 of 6 each) = True.

### H2: held-out robustness (level 1)

| State | ψ̂ | G_DD (S/cm) | G_SD (S/cm) | G_DB (S/cm) | (a) G_DD forms | (b) KCL | (c) B0 bracket | H2 |
|---|---:|---:|---:|---:|---:|---:|---|---|
| HO-1 | 14.993 | 6.177613e-11 | 6.175267e-11 | 2.346e-14 | 6.0e-11 | 6.0e-11 | ok | **pass** |
| HO-2 | 17.993 | 1.121119e-09 | 1.121096e-09 | 2.347e-14 | 2.7e-10 | 2.7e-10 | ok | **pass** |
| HO-3 | 15.000 | 4.569247e-08 | 4.568738e-08 | 5.097e-12 | 2.5e-11 | 2.5e-11 | ok | **pass** |
| HO-4 | 18.000 | 8.401174e-07 | 8.401123e-07 | 5.114e-12 | 5.8e-11 | 5.8e-11 | ok | **pass** |
| HO-5 | 11.994 | 1.717058e-12 | 1.683965e-12 | 3.309e-14 | 4.3e-10 | 4.3e-10 | ok | **pass** |

H2 summary: 5 of 5 pass; complete = True; kill = False.

### P-DD and E-CIRC: drift–diffusion at V_D = ±0.1 and ±1 mV (relative error vs the P1 network G_DD on the same mesh)

| Run | G·1 mV (A/cm) | predicted | DEVSIM drain ±1 mV | ±0.1 mV | Richardson | P-DD | circuit ±1 mV | circuit Richardson | weighting Richardson | E-CIRC |
|---|---:|---|---:|---:|---:|---|---:|---:|---:|---|
| AR-1 B1 L1 | 3.5e-15 | > 1e-2 | -3.46e+01 | -1.65e+03 | -1.66e+03 | confirmed | +2.49e-04 | -7.2e-10 | -6.7e-10 | yes |
| AR-1 B1 L2 | 3.5e-15 | > 1e-2 | -5.70e+01 | -1.56e+03 | -1.57e+03 | confirmed | +2.49e-04 | -1.1e-09 | -1.2e-09 | yes |
| AR-1 B0 L1 | 3.5e-15 | > 1e-2 | -3.47e+01 | -1.65e+03 | -1.67e+03 | confirmed | +2.49e-04 | -7.2e-10 | -7.2e-10 | yes |
| AR-2 B1 L1 | 1.6e-13 | > 1e-2 | -1.73e+00 | -3.69e+01 | -3.72e+01 | confirmed | +2.49e-04 | -2.0e-10 | -1.9e-10 | yes |
| AR-2 B1 L2 | 1.6e-13 | > 1e-2 | -2.22e+00 | -3.50e+01 | -3.53e+01 | confirmed | +2.49e-04 | -2.1e-10 | -2.1e-10 | yes |
| AR-2 B0 L1 | 1.6e-13 | > 1e-2 | -1.73e+00 | -3.69e+01 | -3.72e+01 | confirmed | +2.49e-04 | -2.0e-10 | -1.9e-10 | yes |
| AR-3 B1 L1 | 2.1e-11 | > 1e-2 | +1.01e-01 | -1.27e+00 | -1.29e+00 | confirmed | +2.49e-04 | -1.9e-10 | -1.9e-10 | yes |
| AR-3 B1 L2 | 2.1e-11 | > 1e-2 | +1.49e-02 | -1.26e+00 | -1.27e+00 | confirmed | +2.49e-04 | -1.9e-10 | -1.9e-10 | yes |
| AR-3 B0 L1 | 2.1e-11 | > 1e-2 | +1.01e-01 | -1.27e+00 | -1.29e+00 | confirmed | +2.49e-04 | -1.9e-10 | -1.9e-10 | yes |
| AR-4 B1 L1 | 2.6e-12 | > 1e-2 | -1.10e+00 | -1.48e+00 | -1.48e+00 | confirmed | +2.46e-04 | -1.9e-10 | -1.9e-10 | yes |
| AR-4 B1 L2 | 2.6e-12 | > 1e-2 | -7.20e-01 | -1.45e+00 | -1.46e+00 | confirmed | +2.46e-04 | -1.9e-10 | -1.9e-10 | yes |
| AR-4 B0 L1 | 2.5e-12 | > 1e-2 | -1.10e+00 | -1.48e+00 | -1.48e+00 | confirmed | +2.46e-04 | -1.9e-10 | -1.8e-10 | yes |
| AR-5 B1 L1 | 1.2e-10 | none | +3.57e-02 | -9.70e-01 | -9.80e-01 | — | +2.46e-04 | -1.8e-10 | -1.8e-10 | yes |
| AR-5 B1 L2 | 1.2e-10 | none | +3.85e-02 | -9.15e-01 | -9.24e-01 | — | +2.46e-04 | -1.8e-10 | -1.8e-10 | yes |
| AR-5 B0 L1 | 1.2e-10 | none | +3.58e-02 | -9.70e-01 | -9.80e-01 | — | +2.46e-04 | -1.8e-10 | -1.8e-10 | yes |
| AR-6 B1 L1 | 1.6e-08 | ≤ 1e-2 | +2.53e-04 | +1.54e-03 | +1.55e-03 | confirmed | +2.44e-04 | -1.8e-10 | -1.8e-10 | yes |
| AR-6 B1 L2 | 1.6e-08 | ≤ 1e-2 | +3.78e-04 | -9.46e-04 | -9.60e-04 | confirmed | +2.44e-04 | -1.8e-10 | -1.8e-10 | yes |
| AR-6 B0 L1 | 1.6e-08 | ≤ 1e-2 | +2.54e-04 | +1.54e-03 | +1.56e-03 | confirmed | +2.44e-04 | -1.8e-10 | -1.8e-10 | yes |
| HO-1 B1 L1 | 6.2e-14 | > 1e-2 | -2.92e+00 | -9.51e+01 | -9.60e+01 | confirmed | +2.49e-04 | -2.0e-10 | -2.1e-10 | yes |
| HO-2 B1 L1 | 1.1e-12 | > 1e-2 | -1.11e+00 | -6.18e+00 | -6.23e+00 | confirmed | +2.49e-04 | -1.9e-10 | -1.9e-10 | yes |
| HO-3 B1 L1 | 4.6e-11 | > 1e-2 | -6.38e-02 | -8.96e-01 | -9.04e-01 | confirmed | +2.46e-04 | -1.8e-10 | -1.8e-10 | yes |
| HO-4 B1 L1 | 8.4e-10 | none | +2.79e-03 | +3.46e-02 | +3.49e-02 | — | +2.45e-04 | -1.8e-10 | -1.8e-10 | yes |
| HO-5 B1 L1 | 1.7e-15 | > 1e-2 | -7.01e+01 | -3.39e+03 | -3.42e+03 | confirmed | +2.49e-04 | -1.2e-09 | -8.2e-10 | yes |

P-DD by state: AR-1 confirmed; AR-2 confirmed; AR-3 confirmed; AR-4 confirmed; AR-5 no prediction; AR-6 confirmed; HO-1 confirmed; HO-2 confirmed; HO-3 confirmed; HO-4 no prediction; HO-5 confirmed.  E-CIRC: all = True, worst |rel| = 1.2e-09.

### P-SS: DEVSIM coupled small-signal at level 1 (relative error vs the P1 network G_DD)

| State | body | G_E | G_M | κ_G | E_dwr | predicted > 1e-2 | outcome |
|---|---|---:|---:|---:|---:|---|---|
| AR-1 | B1 | -6.28e+03 | +1.74e+03 | 6.7e+16 | 0.83 | yes | confirmed |
| AR-1 | B0 | -6.51e+03 | +1.64e+03 | 6.5e+16 | 0.81 | no | — |
| AR-2 | B1 | -1.37e+02 | +4.17e+01 | 6.8e+16 | 0.67 | yes | confirmed |
| AR-2 | B0 | -1.33e+02 | +4.25e+01 | 7.0e+16 | 0.75 | no | — |
| AR-3 | B1 | -1.06e+00 | +3.13e-01 | 1.2e+18 | 14.37 | yes | confirmed |
| AR-3 | B0 | -1.05e+00 | +3.12e-01 | 1.4e+18 | 16.00 | no | — |
| AR-4 | B1 | -2.26e+01 | +2.51e+01 | 1.4e+17 | 1.07 | no | — |
| AR-4 | B0 | -2.29e+01 | +2.52e+01 | 1.4e+17 | 1.08 | no | — |
| AR-5 | B1 | -4.61e-01 | +5.53e-01 | 1.2e+17 | 1.00 | no | — |
| AR-5 | B0 | -4.54e-01 | +5.50e-01 | 1.2e+17 | 0.93 | no | — |
| AR-6 | B1 | -3.58e-03 | +4.58e-03 | 4.9e+14 | 0.00 | no | — |
| AR-6 | B0 | -3.57e-03 | +4.75e-03 | 4.9e+14 | 0.00 | no | — |
| HO-1 | B1 | -3.63e+02 | +1.12e+02 | 6.7e+16 | 0.87 | yes | confirmed |
| HO-1 | B0 | -3.71e+02 | +1.05e+02 | 6.5e+16 | 0.89 | no | — |
| HO-2 | B1 | -1.95e+01 | +6.18e+00 | 7.2e+16 | 0.79 | yes | confirmed |
| HO-2 | B0 | -1.99e+01 | +6.28e+00 | 7.0e+16 | 0.88 | no | — |
| HO-3 | B1 | -1.40e+00 | +1.59e+00 | 4.2e+17 | 3.52 | no | — |
| HO-3 | B0 | -1.42e+00 | +1.63e+00 | 4.0e+17 | 3.70 | no | — |
| HO-4 | B1 | -7.53e-02 | +8.09e-02 | 9.9e+15 | 0.09 | no | — |
| HO-4 | B0 | -7.43e-02 | +8.16e-02 | 9.9e+15 | 0.08 | no | — |
| HO-5 | B1 | -7.00e+03 | +7.79e+03 | 1.2e+17 | 1.24 | yes | confirmed |
| HO-5 | B0 | -7.13e+03 | +7.07e+03 | 1.2e+17 | 0.84 | no | — |

P-SS summary: 6 of 6 predicted cases confirmed.

### ψ̂ check (P1, B1, level 1)

| AR-1 | AR-2 | AR-3 | AR-4 | AR-5 | AR-6 | HO-1 | HO-2 | HO-3 | HO-4 | HO-5 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 11.994 (target 12) | 15.993 (target 16) | 21.013 (target 21) | 12.000 (target 12) | 16.000 (target 16) | 21.020 (target 21) | 14.993 (target 15) | 17.993 (target 18) | 15.000 (target 15) | 18.000 (target 18) | 11.994 (target 12) |

Ledger: 84 runs, failures ['P2_AR-4_B1_lev2', 'P2_AR-4_B1_lev1', 'P2_AR-4_B0_lev1'], total wall 4.12 h, max RSS 5957 MB.

### Unsuccessful attempts (kept in the ledger)

| Job | start (UTC) | exit | wall (s) | peak RSS (MB) | cause |
|---|---|---:|---:|---:|---|
| P2_AR-4_B1_lev2 | 2026-11-05T01:48:45+00:00 | -9 | 629.8 | 5570 | killed by the memory cgroup (OOM) |
| P2_AR-4_B1_lev1 | 2026-11-05T02:35:29+00:00 | -9 | 355.3 | 5957 | killed by the memory cgroup (OOM) |
| P2_AR-4_B0_lev1 | 2026-11-05T02:41:25+00:00 | -15 | 60.6 | 2144 | stopped by me (SIGTERM) to prevent a cascade of OOM kills (amendment A9) |

Runs with `MALLOC_MMAP_THRESHOLD_=131072` (amendment A9): 31 of 84 ledger lines.

### Post-hoc: absolute error of DEVSIM post-solve contact currents (max over ±0.1, ±1 mV; vs the circuit-node current)

| Mesh | runs | drain (A/cm) | source (A/cm) | G_DD of those runs (S/cm) |
|---|---:|---|---|---|
| dev1_L1_lev1 | 8 | 1.6e-11 – 1.8e-11 | 1.5e-11 – 1.7e-11 | 3.5e-12 – 2.1e-08 |
| dev1_L1_lev2 | 3 | 2.5e-11 – 2.7e-11 | 2.5e-11 – 2.7e-11 | 3.5e-12 – 2.1e-08 |
| dev1_L2_lev1 | 1 | 1.6e-11 – 1.6e-11 | 1.5e-11 – 1.5e-11 | 1.7e-12 – 1.7e-12 |
| dev2_L2_lev1 | 8 | 3.3e-12 – 1.6e-11 | 2.4e-12 – 1.6e-11 | 2.5e-09 – 1.6e-05 |
| dev2_L2_lev2 | 3 | 3.1e-12 – 1.6e-11 | 4.7e-12 – 1.6e-11 | 2.6e-09 – 1.6e-05 |

### Post-hoc: DEVSIM raw contact currents identical across gate voltages (same mesh, body and bias)

| Mesh | body | V_D | current | runs identical / runs | value (A/cm) | G_DD span of those runs (S/cm) |
|---|---|---:|---|---|---:|---|
| dev1_L1_lev1 | B0 | +1e-04 | I_drain | 3/3 | +1.440079e-11 | 3.52e-12 – 2.11e-08 |
| dev1_L1_lev1 | B0 | +1e-04 | I_source | 3/3 | +1.515618e-11 | 3.52e-12 – 2.11e-08 |
| dev1_L1_lev1 | B0 | -1e-04 | I_drain | 3/3 | +1.556295e-11 | 3.52e-12 – 2.11e-08 |
| dev1_L1_lev1 | B0 | -1e-04 | I_source | 3/3 | +1.515618e-11 | 3.52e-12 – 2.11e-08 |
| dev1_L1_lev1 | B0 | +1e-03 | I_drain | 2/3 | +1.583860e-11 | 3.52e-12 – 1.62e-10 |
| dev1_L1_lev1 | B0 | +1e-03 | I_source | 2/3 | +1.515618e-11 | 3.52e-12 – 1.62e-10 |
| dev1_L1_lev1 | B0 | -1e-03 | I_drain | 2/3 | +1.607576e-11 | 3.52e-12 – 1.62e-10 |
| dev1_L1_lev1 | B0 | -1e-03 | I_source | 2/3 | +1.515618e-11 | 3.52e-12 – 1.62e-10 |
| dev1_L1_lev1 | B1 | +1e-04 | I_drain | 5/5 | +1.440079e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev1 | B1 | +1e-04 | I_source | 5/5 | +1.515618e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev1 | B1 | -1e-04 | I_drain | 5/5 | +1.556295e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev1 | B1 | -1e-04 | I_source | 5/5 | +1.515618e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev1 | B1 | +1e-03 | I_drain | 4/5 | +1.583860e-11 | 3.53e-12 – 1.12e-09 |
| dev1_L1_lev1 | B1 | +1e-03 | I_source | 4/5 | +1.515618e-11 | 3.53e-12 – 1.12e-09 |
| dev1_L1_lev1 | B1 | -1e-03 | I_drain | 4/5 | +1.607576e-11 | 3.53e-12 – 1.12e-09 |
| dev1_L1_lev1 | B1 | -1e-03 | I_source | 4/5 | +1.515618e-11 | 3.53e-12 – 1.12e-09 |
| dev1_L1_lev2 | B1 | +1e-04 | I_drain | 3/3 | +2.405542e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev2 | B1 | +1e-04 | I_source | 3/3 | +2.528419e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev2 | B1 | -1e-04 | I_drain | 3/3 | +2.515295e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev2 | B1 | -1e-04 | I_source | 3/3 | +2.528419e-11 | 3.53e-12 – 2.11e-08 |
| dev1_L1_lev2 | B1 | +1e-03 | I_drain | 2/3 | +2.313104e-11 | 3.53e-12 – 1.62e-10 |
| dev1_L1_lev2 | B1 | +1e-03 | I_source | 2/3 | +2.528419e-11 | 3.53e-12 – 1.62e-10 |
| dev1_L1_lev2 | B1 | -1e-03 | I_drain | 2/3 | +2.352599e-11 | 3.53e-12 – 1.62e-10 |
| dev1_L1_lev2 | B1 | -1e-03 | I_source | 2/3 | +2.528419e-11 | 3.53e-12 – 1.62e-10 |

### Post-hoc: which G_DD form carries the round-off (body B0, cuts through the channel)

| Run | nodes | energy / (−I_S) − 1 | source-half channel cuts vs −I_S (max) | cuts vs energy (range) | max c_e / G |
|---|---:|---:|---:|---|---:|
| AR-1_B0_lev0 | 21924 | +1.05e-10 | 1.0e-11 | [-1.15e-10, -9.84e-11] | 4.7e+17 |
| AR-1_B0_lev1 | 43792 | +4.75e-10 | 2.8e-11 | [-4.82e-10, -4.47e-10] | 5.0e+17 |
| AR-1_B0_lev2 | 87007 | +1.10e-09 | 1.5e-10 | [-1.24e-09, -1.08e-09] | 5.1e+17 |
| AR-2_B0_lev0 | 21924 | -6.53e-11 | 4.0e-11 | [+2.52e-11, +6.67e-11] | 1.0e+16 |
| AR-2_B0_lev1 | 43792 | +1.16e-10 | 1.0e-10 | [-1.22e-10, -1.48e-11] | 1.1e+16 |
| AR-4_B0_lev0 | 69353 | -4.41e-11 | 3.8e-11 | [+5.64e-12, +4.44e-11] | 6.5e+14 |
| AR-4_B0_lev1 | 137970 | -2.50e-11 | 2.5e-11 | [-4.92e-13, +2.47e-11] | 6.7e+14 |
