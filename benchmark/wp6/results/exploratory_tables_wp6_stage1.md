### Exploratory description, stage 1 (post hoc; decides no outcome)

**Per task**

| Task | state | class | G_ref (S/cm) | outcomes | max \|rel. error\| | wall (s) |
|---|---|---|---:|---|---:|---|
| T01 | AR-1 | low | 3.5326e-12 | CS CS CS | 1.1e-04 | 254, 219, 246 |
| T02 | AR-2 | low | 1.6192e-10 | CS CS CS | 2.6e-05 | 282, 158, 363 |
| T03 | AR-3 | low | 2.1147e-08 | CS CS CS | 1.1e-05 | 307, 293, 234 |
| T04 | AR-4 | low | 2.5505e-09 | CS CS CS | 2.8e-06 | 2488, 528, 458 |
| T05 | AR-5 | moderate | 1.2035e-07 | CS CS CS | 1.5e-05 | 869, 766, 1029 |
| T06 | AR-6 | moderate | 1.5976e-05 | CS CS NR | 2.1e-05 | 2658, 586, 2702 |
| T07 | HO-1 | low | 6.1776e-11 | CS CS CS | 9.9e-05 | 264, 233, 553 |
| T08 | HO-2 | low | 1.1211e-09 | CS CS CS | 1.7e-05 | 804, 264, 282 |
| T09 | HO-3 | low | 4.5692e-08 | CS CS CS | 5.4e-05 | 386, 890, 560 |
| T10 | HO-4 | moderate | 8.4012e-07 | CS CS CS | 2.1e-05 | 587, 1271, 671 |
| T11 | HO-5 | low | 1.7171e-12 | CS CS CS | 3.4e-05 | 154, 243, 200 |

Reported values: 32; median |rel. error| 1.3e-05, largest 1.1e-04; 31 within 1e-4 and 32 within 1e-3.

**By task class**

| Class | runs | outcomes | median wall (s) | median turns | cost (USD) |
|---|---:|---|---:|---:|---:|
| low | 24 | {'CS': 24} | 282 | 8 | 13.21 |
| moderate | 9 | {'CS': 8, 'NR': 1} | 869 | 8 | 4.89 |

**Paths the agents' code touched** (runs whose tool inputs contain the call; mechanical count)

| Path | runs |
|---|---:|
| DEVSIM contact current (get_contact_current) | 33/33 |
| AC small-signal solve (type='ac') | 30/33 |
| circuit current V1.I (get_circuit_node_value) | 33/33 |
| drain bias steps (circuit_alter) | 33/33 |
| extended precision (extended_solver/model/equation) | 6/33 |
| Jacobian export (get_matrix_and_rhs) | 4/33 |
| Richardson / extrapolation (by name) | 1/33 |
| electron density read back (Electrons) | 10/33 |

**Final method** (analyst coding of each run's reported method)

| Category | runs | ids |
|---|---:|---|
| A: central differences of a converged DC current, double precision | 25 | W6-001, W6-002, W6-004, W6-007, W6-008, W6-009, W6-011, W6-013, W6-015, W6-016, W6-017, W6-018, W6-019, W6-020, W6-021, W6-022, W6-023, W6-025, W6-027, W6-028, W6-029, W6-030, W6-031, W6-032, W6-033 |
| B: central differences of a converged DC current, extended precision | 2 | W6-010, W6-014 |
| C: exact linearisation about the equilibrium: Scharfetter-Gummel resistor network on the mesh edges | 5 | W6-003, W6-005, W6-006, W6-012, W6-024 |
| D: no report | 1 | W6-026 |

**Refused tool calls** (by the mechanism that refused them)

| Mechanism | calls |
|---|---:|
| classifier: Classifier unavailable | 9 |
| workspace guard | 7 |
| Claude Code's command check | 3 |

Guard refusals by target: {'own background-task output': 5, 'code text matched the home-directory rule': 2}. Runs with any refusal: W6-004, W6-005, W6-006, W6-024, W6-027, W6-033.
