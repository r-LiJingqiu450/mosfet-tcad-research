### D1 (post hoc): what limits the extended-precision agreement (X2) at the 1e-10 level

R = Richardson value of the DEVSIM contact current / same-run network value − 1. a = coefficient of x² (x = V_D/V_T) from the ±0.1 mV value. A sinh-shaped response (diffusion-limited transport) has a = 1/6 = 0.16667 and leaves R_sinh = −x₁⁴/12000 = -1.866e-10; R − R_sinh is then minus the relative error of the network value. The forms indicator is WP2's relative difference of the two network forms on the same mesh.

| State | a | R, G_DD | R, G_SD | R − R_sinh, G_DD | R − R_sinh, G_SD | WP2 forms indicator G_DD / G_SD | R, circuit / weighting-function current |
|---|---:|---:|---:|---:|---:|---|---|
| AR-1 | 0.16643 | -6.16e-10 | -2.44e-10 | -4.3e-10 | -5.7e-11 | -4.6e-10 / -6.6e-15 | -6.16e-10 / -6.16e-10 |
| AR-2 | 0.16644 | -1.93e-10 | -9.08e-11 | -6.2e-12 | +9.6e-11 | -1.4e-10 / +6.7e-16 | -1.93e-10 / -1.93e-10 |
| AR-3 | 0.16641 | -1.86e-10 | -8.68e-11 | +4.0e-13 | +1.0e-10 | -1.1e-10 / -2.7e-14 | -1.86e-10 / -1.86e-10 |
| AR-4 | 0.16448 | -1.85e-10 | -1.96e-10 | +1.4e-12 | -9.0e-12 | +3.6e-11 / +8.0e-14 | -1.85e-10 / -1.85e-10 |
| AR-5 | 0.16417 | -1.82e-10 | -2.45e-10 | +4.8e-12 | -5.8e-11 | +1.7e-11 / -1.4e-14 | -1.82e-10 / -1.82e-10 |
| AR-6 | 0.16334 | -1.79e-10 | -2.04e-10 | +7.4e-12 | -1.7e-11 | +2.0e-12 / +1.1e-13 | -1.79e-10 / -1.79e-10 |

### D2 (post hoc): WP2 level-1 reference values recomputed on extended-precision equilibria

| State | G_DD shift (P2X) | G_DD shift (AC ext.) | WP2 forms indicator G_DD | G_SD shift (P2X) | WP2 forms indicator G_SD | double-precision reruns |
|---|---:|---:|---:|---:|---:|---|
| AR-1 | -7.7e-11 | -7.7e-11 | -4.6e-10 | +1.0e-11 | -6.6e-15 | +0.0e+00, +0.0e+00 |
| AR-2 | -1.1e-12 | -1.1e-12 | -1.4e-10 | +3.4e-11 | +6.7e-16 | +0.0e+00, +0.0e+00 |
| AR-3 | +2.2e-16 | +2.2e-16 | -1.1e-10 | +1.2e-11 | -2.7e-14 | +0.0e+00, +0.0e+00 |
| AR-4 | +7.1e-13 | +7.1e-13 | +3.6e-11 | -2.5e-11 | +8.0e-14 | +0.0e+00, +0.0e+00 |
| AR-5 | -7.4e-15 | -7.4e-15 | +1.7e-11 | +4.6e-11 | -1.4e-14 | +0.0e+00, +0.0e+00 |
| AR-6 | -4.4e-16 | -4.4e-16 | +2.0e-12 | +2.3e-11 | +1.1e-13 | +0.0e+00, +0.0e+00 |

### D3 (post hoc): extended-precision current errors against float64 resolution

| State | e_D (A/cm) | e_S (A/cm) | ulp of largest current (A/cm) | e_S / ulp |
|---|---:|---:|---:|---:|
| AR-1 | 6.3e-30 | 2.0e-30 | 7.9e-31 | 2.53 |
| AR-2 | 2.5e-29 | 2.3e-29 | 2.5e-29 | 0.91 |
| AR-3 | 0.0e+00 | 1.2e-27 | 3.2e-27 | 0.36 |
| AR-4 | 0.0e+00 | 1.3e-28 | 4.0e-28 | 0.33 |
| AR-5 | 0.0e+00 | 9.1e-27 | 2.6e-26 | 0.35 |
| AR-6 | 0.0e+00 | 5.3e-25 | 3.3e-24 | 0.16 |

### D4 (post hoc): DEVSIM built-in AC (double) against the WP2 small-signal `G_M` (same Jacobian, UMFPACK through `devsim.umfpack`)

| State | G_AC double (S/cm) | WP2 G_M (S/cm) | bitwise equal | relative difference | WP2 κ_G | WP2 E_round estimate |
|---|---:|---:|---|---:|---:|---:|
| AR-1 | 6.145795e-09 | 6.145795e-09 | yes | +0.0e+00 | 6.7e+16 | 2.4e-01 |
| AR-2 | 6.910470e-09 | 6.910470e-09 | yes | +0.0e+00 | 6.8e+16 | 2.4e-01 |
| AR-3 | 2.776593e-08 | 2.776593e-08 | yes | +0.0e+00 | 1.2e+18 | 4.1e+00 |
| AR-4 | 6.652069e-08 | 6.652069e-08 | yes | +0.0e+00 | 1.4e+17 | 2.6e-01 |
| AR-5 | 1.868963e-07 | 1.869280e-07 | **no** | -1.7e-04 | 1.2e+17 | 2.2e-01 |
| AR-6 | 1.604885e-05 | 1.604885e-05 | yes | +0.0e+00 | 4.9e+14 | 8.8e-04 |

### D5 (post hoc): G_SD from the DEVSIM source contact current (Richardson) vs the same-run network G_SD

| State | double (WP2 P2) | extended (WP4 P2X) |
|---|---:|---:|
| AR-1 | -1.00e+00 | -2.44e-10 |
| AR-2 | -1.00e+00 | -9.08e-11 |
| AR-3 | -1.01e+00 | -8.68e-11 |
| AR-4 | -1.00e+00 | -1.96e-10 |
| AR-5 | -9.40e-01 | -2.45e-10 |
| AR-6 | +6.36e-04 | -2.04e-10 |

### D6 (post hoc): cost (ledger wall time and peak RSS)

| State | P2 double | P2 extended | ratio | Newton iterations double / extended | AC double (solve) | AC extended (solve) |
|---|---:|---:|---:|---:|---:|---:|
| AR-1 | 55 s, 529 MB | 598 s, 484 MB | 10.9× | 43 / 43 | 20 s (1.6 s), 541 MB | 190 s (10.9 s), 578 MB |
| AR-2 | 55 s, 531 MB | 598 s, 483 MB | 10.9× | 43 / 43 | 20 s (1.6 s), 524 MB | 190 s (10.8 s), 578 MB |
| AR-3 | 55 s, 505 MB | 597 s, 485 MB | 10.8× | 43 / 43 | 20 s (1.6 s), 521 MB | 190 s (10.9 s), 578 MB |
| AR-4 | 379 s, 1221 MB | 5961 s, 1428 MB | 15.7× | 44 / 44 | 133 s (11.3 s), 1600 MB | 1959 s (106.4 s), 1808 MB |
| AR-5 | 382 s, 1221 MB | 5966 s, 1428 MB | 15.6× | 44 / 44 | 132 s (11.1 s), 1600 MB | 1957 s (106.7 s), 1808 MB |
| AR-6 | 383 s, 1221 MB | 5956 s, 1428 MB | 15.6× | 44 / 44 | 132 s (11.3 s), 1601 MB | 1955 s (106.8 s), 1807 MB |

### D7 (post hoc): P2X_AR-4, killed first attempt vs rerun

The first attempt's log has 582 lines; 582 of them are identical to the rerun's first 582 lines (attempt 1 is a prefix of the rerun). The rerun adds 3 lines: the summary, the output path and the exit code.
