### WP4-A: H3, DEVSIM default drain contact current scored with `tools/score.py`

| Tier | result |
|---|---|
| A-R | fail 6 |
| B-R | fail 6 |
| HO | fail 5 |

H3 kill (every external-path case passes trivially): False.

### WP4-B1: DEVSIM extended precision, post-solve contact current (level 1, body B1; paired with the WP2 double-precision run)

| State | ψ̂ | e_D double (A/cm) | e_D extended | e_S extended | reduction | DEVSIM drain G, Richardson: double → extended | circuit Richardson (ext) | X1 | X2 | X3 | wall double → extended |
|---|---:|---:|---:|---:|---:|---|---:|---|---|---|---|
| AR-1 | 12 | 1.6e-11 | 6.3e-30 | 2.0e-30 | 2.5e+18 | -1.66e+03 → -6.16e-10 | -6.2e-10 | yes | yes | yes | 55 s → 597 s |
| AR-2 | 16 | 1.6e-11 | 2.5e-29 | 2.3e-29 | 6.4e+17 | -3.72e+01 → -1.93e-10 | -1.9e-10 | yes | yes | yes | 54 s → 597 s |
| AR-3 | 21 | 1.8e-11 | 0.0e+00 | 1.2e-27 | inf | -1.29e+00 → -1.86e-10 | -1.9e-10 | yes | yes | yes | 55 s → 596 s |
| AR-4 | 12 | 7.1e-12 | 0.0e+00 | 1.3e-28 | inf | -1.48e+00 → -1.85e-10 | -1.9e-10 | yes | yes | yes | 379 s → 5960 s |
| AR-5 | 16 | 1.6e-11 | 0.0e+00 | 9.1e-27 | inf | -9.80e-01 → -1.82e-10 | -1.8e-10 | yes | yes | yes | 382 s → 5965 s |
| AR-6 | 21 | 3.3e-12 | 0.0e+00 | 5.3e-25 | inf | +1.55e-03 → -1.79e-10 | -1.8e-10 | yes | yes | yes | 383 s → 5955 s |

B1 summary: X1 6/6, X2 6/6, X3 6/6.
Network G_DD from the extended equilibrium vs the double equilibrium: |rel| ≤ 7.7e-11; final Newton updates ≤ 2.2e-34 V.

### WP4-B2: DEVSIM built-in AC analysis at 1e-6 Hz (level 1, body B1), G_AC vs the same-run network G_DD

| State | ψ̂ | double | extended | WP2 P3 `G_M` (UMFPACK) | predicted > 1e-2 | Y1 (double) | Y2 (extended) |
|---|---:|---:|---:|---:|---|---|---|
| AR-1 | 12 | +1.74e+03 | +3.26e+03 | +1.74e+03 | yes | confirmed | confirmed |
| AR-2 | 16 | +4.17e+01 | +6.74e+01 | +4.17e+01 | yes | confirmed | confirmed |
| AR-3 | 21 | +3.13e-01 | +5.61e-01 | +3.13e-01 | yes | confirmed | confirmed |
| AR-4 | 12 | +2.51e+01 | +2.95e+01 | +2.51e+01 | no | — | — |
| AR-5 | 16 | +5.53e-01 | +5.87e-01 | +5.53e-01 | no | — | — |
| AR-6 | 21 | +4.58e-03 | +4.27e-03 | +4.58e-03 | no | — | — |

B2 summary: Y1 3/3, Y2 3/3.
