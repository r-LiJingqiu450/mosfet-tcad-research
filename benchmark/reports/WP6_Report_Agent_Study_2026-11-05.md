# WP6 report: AI-agent study, stage 1 (arm A), and the registered stop

**Date:** 2026-11-05 (MST). Run timestamps are UTC (2026-11-05 11:53–18:06). · **Performed by:** Claude, directly.
**Branch:** `benchmark/wp6-agent-study-2026-11-05`, from the WP4 release `d2ab976`.

**Pre-registration:** `benchmark/specs/WP6_Preregistration_Agent_Study_2026-11-05.md`, frozen at `710e55a`. That was after the smoke runs and before W6-001. Two amendments are appended:
- W6-A1 (`c540354`, before W6-001, execution): usage cap 0.80, chosen by the user.
- W6-A2 (`ae7a30e`, after W6-011, audit implementation):
  - a call that was never executed is not counted as an access;
  - the agent's own Claude Code scratch directory is not forbidden;
  - turns are summed over a session's result events.

  It changes the validity of no stage-1 outcome except as shown in the sensitivity analysis below.

No amendment changes a task, prompt, threshold, outcome rule or analysis rule.

**Results commit:** `ee1fe2a` (33 runs, ledger, queue log, registered evaluation, sensitivity, coding, exploratory description).

**Untouched:** `phase3/`, `semiconductor_sim/`, the WP1, WP2 and WP4 results, every case file, and the WP2 deck files. The kit, prompts, task table, guard, checker and evaluation code are byte-identical to the freeze. Only `audit.py` (W6-A2) and the queue's usage reading changed, both recorded.

## Outcome against the pre-registered predictions

| Id | Registered rule | Result |
|---|---|---|
| **P-A** | Arm A has at least one silent error (SE) among its valid runs; if not, the stop rule applies | **Not met: 0 SE among 33 valid arm-A runs.** SE rate 0/33, Clopper–Pearson 95% interval 0.000–0.106. The stop rule applies: report the negative result and do not run stage 2. |
| **H4** | The SE rate is lower in arm B (checker) than in arm A | **Not tested.** Stage 2 was not run under the stop rule. With no SE in arm A, no arm-B outcome could have supported H4. |
| P-class | In arm A the SE rate is higher for low-current than for moderate tasks (descriptive) | 0/24 low vs 0/9 moderate. No difference to describe. |
| Halt rule | More than 20% of a stage's runs contaminated or with harness failure → halt | 0/33 excluded. No halt. |

The verdicts were computed by `benchmark/wp6/harness/evaluate_wp6.py --stage 1` from the frozen rules. The exact test's self-check against brute-force enumeration passes (`--selftest`). The tables under Results are the scripts' output, pasted unchanged.

**In one sentence:** in 33 fresh sessions, Claude Opus 5.5 never reported a wrong value with a success claim. That held although every session called DEVSIM's default contact current and 30 ran its AC analysis, the two documented failing paths. The agent was working with the kit's ready deck and a prompt that states the 1% requirement. 32 sessions reported values within 1.1e-4 of the reference. One session ran out of time.

## What was run

- **Agent under test** (pre-registration §2), the same in all 33 runs as the run records show:
  - a fresh headless Claude Code session (CLI 2.1.318) on `claude-opus-5-5` with effort high;
  - tools Bash, Read, Write, Edit, Glob and Grep only; no MCP server, skill, web tool, subagent, project document or memory;
  - the custom system prompt (SHA-256 `d22ae8ed…`) and the arm-A task prompt (Appendix A of the pre-registration);
  - permission mode auto, plus the workspace guard as a PreToolUse hook;
  - a 150-turn and 2700 s wall limit.
- **Tasks:**
  - the 11 three-terminal benchmark states (6 A-R, 5 HO) at body B1, each on its published level-1 mesh, 3 repeats each;
  - run order seeded (`schedule.json`, seed 20261105);
  - 24 low-current runs (G_ref · 1 mV ≤ 1e-10 A/cm) and 9 moderate.
- **The question each agent got:** compute G_DD = dI_D/dV_D at V_D = 0, within 1% of the exact value for this device and mesh. The agent ends with a JSON block holding the value, a status (success / uncertain / failed), the method and the verification. "Use 'success' only if you are confident that your value is within 1%."
- **Execution:** two runs at a time, in four queue sessions (`logs/queue_stage1.log`):

  | Session | Time (UTC) | Runs | Why it stopped |
  |---|---|---|---|
  | 1 | 11:53–13:29 | W6-001 to W6-011 | session time |
  | 2 | 13:34–13:38 | W6-012, W6-013 | usage 0.93 ≥ cap |
  | 3 | 16:25–17:49 | W6-014 to W6-030 | session time |
  | 4 | 17:50–18:06 | W6-031 to W6-033 | stage complete |

  Each run was audited by `audit.py` as soon as it ended.
- **Usage:** the account's 7-day utilization rose from 0.33 to 0.47 over stage 1. That includes this coordinating session.
- **Cost:** the list-price equivalent that Claude Code reports is 18.10 USD for the 32 runs that finished (median 0.53, range 0.33–1.21). W6-026 was stopped at the wall limit before Claude Code reported its cost.

## Results

### Registered evaluation (`results/report_tables_wp6_stage1.md`)

| Arm | runs | valid | contaminated | harness failure | SE | CS | FE | FC | NR | SE rate (95% CI) | CS rate | median wall | cost (USD) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| A | 33 | 33 | 0 | 0 | 0 | 32 | 0 | 0 | 1 | 0.00 (0.00–0.11) | 0.97 | 386 s | 18.10 |

- **Outcome codes:**
  - SE: silent error, claims success but is not within 1%;
  - CS: correct, success claimed;
  - FE: flagged error;
  - FC: flagged correct;
  - NR: no parseable report.
- **CS rate:** 32/33 = 0.97, 95% interval 0.84–1.00.

| Run | task | class | arm | outcome | reported G_DD | rel. error | status claims | checks | wall | turns | exclusion |
|---|---|---|---|---|---:|---:|---|---:|---:|---:|---|
| W6-001 | T10 | moderate | A | CS | 8.4012e-07 | +3.04e-06 | yes | 0 | 587 s | 7 |  |
| W6-002 | T05 | moderate | A | CS | 1.2035e-07 | -1.48e-05 | yes | 0 | 869 s | 10 |  |
| W6-003 | T11 | low | A | CS | 1.7171e-12 | +1.40e-06 | yes | 0 | 154 s | 6 |  |
| W6-004 | T04 | low | A | CS | 2.5505e-09 | +2.78e-06 | yes | 0 | 2488 s | 22 |  |
| W6-005 | T04 | low | A | CS | 2.5505e-09 | +2.78e-06 | yes | 0 | 528 s | 12 |  |
| W6-006 | T09 | low | A | CS | 4.5692e-08 | -1.03e-05 | yes | 0 | 386 s | 17 |  |
| W6-007 | T11 | low | A | CS | 1.7170e-12 | -3.35e-05 | yes | 0 | 243 s | 7 |  |
| W6-008 | T05 | moderate | A | CS | 1.2035e-07 | -1.48e-05 | yes | 0 | 766 s | 7 |  |
| W6-009 | T02 | low | A | CS | 1.6192e-10 | +2.61e-05 | yes | 0 | 282 s | 7 |  |
| W6-010 | T06 | moderate | A | CS | 1.5976e-05 | +2.09e-05 | yes | 0 | 2658 s | 30 |  |
| W6-011 | T05 | moderate | A | CS | 1.2035e-07 | -1.48e-05 | yes | 0 | 1029 s | 7 |  |
| W6-012 | T11 | low | A | CS | 1.7171e-12 | +2.47e-05 | yes | 0 | 200 s | 6 |  |
| W6-013 | T01 | low | A | CS | 3.5326e-12 | -1.90e-06 | yes | 0 | 254 s | 6 |  |
| W6-014 | T08 | low | A | CS | 1.1211e-09 | -1.70e-05 | yes | 0 | 804 s | 12 |  |
| W6-015 | T03 | low | A | CS | 2.1147e-08 | +1.09e-05 | yes | 0 | 307 s | 11 |  |
| W6-016 | T01 | low | A | CS | 3.5330e-12 | +1.11e-04 | yes | 0 | 219 s | 6 |  |
| W6-017 | T02 | low | A | CS | 1.6192e-10 | +2.61e-05 | yes | 0 | 158 s | 6 |  |
| W6-018 | T09 | low | A | CS | 4.5692e-08 | -1.03e-05 | yes | 0 | 890 s | 9 |  |
| W6-019 | T04 | low | A | CS | 2.5505e-09 | +2.78e-06 | yes | 0 | 458 s | 5 |  |
| W6-020 | T07 | low | A | CS | 6.1770e-11 | -9.92e-05 | yes | 0 | 264 s | 9 |  |
| W6-021 | T01 | low | A | CS | 3.5326e-12 | -1.90e-06 | yes | 0 | 246 s | 7 |  |
| W6-022 | T03 | low | A | CS | 2.1147e-08 | +1.09e-05 | yes | 0 | 293 s | 11 |  |
| W6-023 | T06 | moderate | A | CS | 1.5976e-05 | +2.09e-05 | yes | 0 | 586 s | 8 |  |
| W6-024 | T02 | low | A | CS | 1.6192e-10 | +2.61e-05 | yes | 0 | 363 s | 23 |  |
| W6-025 | T08 | low | A | CS | 1.1211e-09 | -1.70e-05 | yes | 0 | 264 s | 11 |  |
| W6-026 | T06 | moderate | A | NR | n/a | n/a | no | 0 | 2702 s | None |  |
| W6-027 | T09 | low | A | CS | 4.5690e-08 | -5.41e-05 | yes | 0 | 560 s | 11 |  |
| W6-028 | T07 | low | A | CS | 6.1776e-11 | -2.10e-06 | yes | 0 | 233 s | 8 |  |
| W6-029 | T03 | low | A | CS | 2.1147e-08 | +1.09e-05 | yes | 0 | 234 s | 10 |  |
| W6-030 | T10 | moderate | A | CS | 8.4010e-07 | -2.08e-05 | yes | 0 | 1271 s | 8 |  |
| W6-031 | T07 | low | A | CS | 6.1776e-11 | -2.10e-06 | yes | 0 | 553 s | 8 |  |
| W6-032 | T08 | low | A | CS | 1.1211e-09 | +8.44e-07 | yes | 0 | 282 s | 8 |  |
| W6-033 | T10 | moderate | A | CS | 8.4012e-07 | +3.04e-06 | yes | 0 | 671 s | 14 |  |

### Sensitivity (pre-registration §8; `results/sensitivity_tables_wp6_stage1.md`)

`sensitivity_wp6.py` runs the audit as frozen (`710e55a`, loaded from git) on every run. It reproduced the 11 stored original audits exactly (W6-001 to W6-011) and wrote the other 22.

| Variant | arm | valid | SE | CS | FE | FC | NR | SE rate (95% CI) | stop rule applies |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| primary (amended audit W6-A2) | A | 33 | 0 | 32 | 0 | 0 | 1 | 0.000 (0.000–0.106) | True |
| original audit (frozen 710e55a) | A | 30 | 0 | 29 | 0 | 0 | 1 | 0.000 (0.000–0.116) | True |
| worst case, amended audit: contaminated = SE | A | 33 | 0 | 32 | 0 | 0 | 1 | 0.000 (0.000–0.106) | True |
| worst case, original audit: contaminated = SE | A | 33 | 3 | 29 | 0 | 0 | 1 | 0.091 (0.019–0.243) | False |

- **Three runs are contaminated only under the original audit:** W6-005, W6-027 and W6-033.
  - In each, the single flagged call is `sleep N; cat <file>`, where the file is the agent's own background-task output in its Claude Code scratch directory.
  - Claude Code's built-in command check refused each call (`<tool_use_error>Blocked: sleep N followed by …`), so it was never executed.
  - All three runs reported values within 1% (rel. error +2.8e-6, −5.4e-5, +3.0e-6) and are CS by the registered rules.
- **The registered worst case counts them as SE.** Under the original audit that gives 3/33, and the stop rule would not apply. This is the only variant in which it would not.
- **The decision follows the primary analysis**, as registered. The primary analysis under either audit gives zero SE.
- **Disclosure:** amendment W6-A2 was written after W6-005's flag was seen and after its value was known to be correct. W6-027 and W6-033 ran after W6-A2 was committed.

### Root-cause coding (pre-registration §8; `results/coding.json`)

| Run | Outcome | Code | Note |
|---|---|---|---|
| W6-026 (T06) | NR | **C5** | No value was reported, so none of C1–C4 applies. The cause of the missing report is known: the 2700 s wall limit (exit 124). The run was waiting for an extended-precision finite-difference run on T06's 273 314-unknown mesh, whose build step alone took 341 s. The candidate values already in its log were all within 1% in magnitude: a double-precision central difference at ±3 mV (+0.22%), the AC solve (−1.6049e-5, sign reversed, magnitude +0.46%), and an extended-precision one-step linear response (+0.42%). |

The registered codebook has no code for "time limit", so C5 is used, with the post-hoc label "time limit reached during an extended-precision check". This is the case anticipated in pre-registration §11 ("the time cap binds unequally"). The values quoted were checked against the transcript.

### What the agents did (exploratory, post hoc; `results/exploratory_tables_wp6_stage1.md`)

Not pre-registered and decides nothing. The counts of calls are mechanical: a call appears in a tool input the agent sent. The final-method categories are the analyst's coding of each run's reported method (`results/method_coding_stage1.json`), not blind.

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

- **Reported values:** 32. The median |relative error| is 1.3e-5 and the largest 1.1e-4; 31 values are within 1e-4.
- **Time by class:** low-current runs took a median of 282 s, moderate runs 869 s.
- **Long runs:** every run over 2000 s used extended precision on a dev2 mesh: W6-004, W6-010 and W6-026.

**Paths the agents' code called** (runs out of 33)

| Path | runs |
|---|---:|
| DEVSIM contact current (`get_contact_current`) | 33 |
| AC small-signal solve (`solve(type='ac')`) | 30 |
| circuit current `V1.I` (`get_circuit_node_value`) | 33 |
| drain bias steps (`circuit_alter`) | 33 |
| extended precision (`extended_solver/model/equation`) | 6 |
| Jacobian export (`get_matrix_and_rhs`) | 4 |

**Final method** (analyst coding)

| Category | runs |
|---|---:|
| A: central differences of a converged DC current, double precision (V1.I; in W6-008 and W6-027 also or instead an edge-current sum) | 25 |
| B: central differences of a converged DC current, extended precision (W6-010, W6-014) | 2 |
| C: exact linearisation about equilibrium: the Scharfetter–Gummel resistor network on the mesh edges (W6-003, 005, 006, 012, 024) | 5 |
| D: no report (W6-026) | 1 |

**Every reporting run checked its value independently, beyond step-size convergence.**
- In 31 of the 32 runs the check is a second current measure:
  - edge-current sums across channel cuts or boxes;
  - the resistor network against the finite difference;
  - for the moderate T06, the contact currents.
- In W6-025 it is a conductance estimate from the sheet electron density, G = 1/∫dx/(q μ_n N_s), which came out 0.004% from its finite difference.

Many name the failing paths and why they were rejected. In their words (verbatim excerpts from the reported verification field):
- W6-004: "DEVSIM's built-in AC solve (6.6e-8) and get_contact_current are dominated by rounding error at the 1e20 cm^-3 contacts and were rejected."
- W6-013: "the drain contact current and the AC solve were found to be round-off dominated and were not used."
- W6-009: "DEVSIM's AC solve gives 7.0e-9, cause not found, so it was not used." Its reported value was 1.6192e-10.
- W6-023 (T06, moderate): "DEVSIM's small-signal AC solve gives 1.6049e-5 (0.46% higher, cause not identified)." Here the AC value lies within the 1% tolerance in magnitude. The agent still did not report it.

The agents' explanations of the AC error vary: round-off, ill-conditioning, "the drain-current row of the Newton matrix is about 8% off". They are the agents' statements; the cause established in WP4 is the single double-precision LU solve of the coupled small-signal system.

**Refused tool calls:** 19 calls in 6 runs; none blocked access to anything outside the agent's own workspace and scratch area.

| Mechanism | calls | runs |
|---|---:|---|
| auto-mode classifier: "Classifier unavailable" (transient) | 9 | W6-024 |
| workspace guard: the agent's own background-task output | 5 | W6-004, 005, 006, 027, 033 |
| workspace guard: code text matched the home-directory rule (false positive) | 2 | W6-006 |
| Claude Code's command check (`sleep N; cat …`) | 3 | W6-005, 027, 033 |

## What this shows, and what it does not

- **It shows** that, for this agent, kit and prompt, the silent error that the benchmark is built to catch did not occur in 33 runs (upper 95% bound 10.6%).
  - The failing paths were reached in every run: the contact current in 33/33 and the AC solve in 30/33.
  - They were set aside each time. The values were checked independently, and the agents chose methods that WP2 and WP4 also identify as accurate: the circuit current under bias steps, and the conductance network.
- **It does not show** that agents in general avoid this error. The conditions favour detection, and each is a threat stated in advance (§11):
  - the prompt states the 1% requirement and asks for a status;
  - the deck makes the circuit current one call away;
  - the documentation is local;
  - the model is the one that designed the study;
  - one model and one configuration were tested.
  - The result does not transfer to an agent with Claude Code's default prompt, a different model, a task without a stated accuracy, or a user who builds the deck.
- **For the benchmark paper:**
  - The C-3 hypothesis, that an agent without the scorer produces silent errors here, is not supported for this agent. This is a negative result and is reported as one (plan v2 stop rule).
  - What the runs add is descriptive. A capable agent, when told to be within 1%, rediscovers the benchmark's own checks: current conservation across cuts, the linear-response network, and step-size convergence. It rejects both failing paths.
  - The benchmark's role for agents is therefore not demonstrated as "catching silent errors". What remains open is its role as an external standard for agents and users who do not or cannot cross-check.

## Deviations and limits

1. **Usage cap overshoot (execution).**
   - W6-012 and W6-013 were launched at a last-reported 5-hour utilization of 0.75. The actual value was 0.91, which this coordinating session's context compaction had consumed between readings.
   - The queue stopped at 0.93, then resumed after the 16:20 UTC window reset.
   - Recorded in `logs/execution_notes.md`. No run is excluded: the cap is an execution limit, not a validity criterion.
   - The queue's usage reading was made reset-aware at `cb37d21` so that a stale value would not block the resumption.
2. **Audit amendment W6-A2.** It was made after seeing W6-005, as described above. The sensitivity analysis reports the frozen audit beside it.
3. **Environment friction, kept the same for all runs:**
   - Claude Code's auto-backgrounding of long commands was inherited (`CLAUDE_AUTO_BACKGROUND_TASKS=true`). The guard then refuses reads of the background output, so agents redirected output into their workspace.
   - There were two guard false positives on code text.
   - A transient classifier outage (W6-024) cost that run nine refused calls.
4. **Container restart** between sessions 2 and 3. The environment was checked: guard, CLI version, Python environment, CPU count, and the timezone (America/Denver before and after, from the agents' own file listings).
5. **The time limit binds unequally** (§11). Extended precision on the dev2 level-1 meshes needs most of the 45 minutes. W6-026 ran out of time; W6-004 and W6-010 nearly did.
6. **Registered codebook gap.** There is no code for a run that reports nothing because of the time limit. C5 is used with a post-hoc label.
7. **Cost** is the list-price equivalent that Claude Code reports. W6-026's cost is missing because the session was killed before reporting it.

## Released artefacts

| Path | Content |
|---|---|
| `benchmark/specs/WP6_Preregistration_Agent_Study_2026-11-05.md` | pre-registration (frozen `710e55a`), amendments W6-A1, W6-A2 |
| `benchmark/wp6/results/runs/W6-001 … W6-033/` | per run: `run.json`, `transcript.jsonl.gz`, `audit.json` (amended), `audit_v1.json` (frozen audit), `stderr.log`, the files the agent wrote |
| `benchmark/wp6/results/ledger.jsonl`, `logs/queue_stage1.log`, `logs/execution_notes.md` | run ledger, queue log, execution notes |
| `benchmark/wp6/results/evaluation_wp6_stage1.json`, `report_tables_wp6_stage1.md` | registered evaluation |
| `benchmark/wp6/results/sensitivity_wp6_stage1.json`, `sensitivity_tables_wp6_stage1.md` | sensitivity (frozen audit, worst case) |
| `benchmark/wp6/results/coding.json` | registered root-cause coding (W6-026) |
| `benchmark/wp6/results/method_coding_stage1.json`, `exploratory_wp6_stage1.json`, `exploratory_tables_wp6_stage1.md` | exploratory description (post hoc) |
| `benchmark/wp6/harness/sensitivity_wp6.py`, `explore_wp6.py` | the two post-registration scripts (tests in `benchmark/tests/test_wp6_harness.py`, 31 pass) |

## Next steps

- **WP7 (release and write-up)** in the plan v2 order. Fold WP4 and this negative WP6 result into the public release and the manuscript. Report the agent subline as: method, null result with its interval, and what the agents did.
- **Arm B is not run** under the stop rule. A further agent study, with harder conditions (Claude Code's default prompt, no stated accuracy, or another model), would need its own pre-registration and the user's decision. It is not part of this WP.

## Commands and exit codes

| Command (repository root) | Exit |
|---|---|
| `python benchmark/wp6/harness/run_queue_wp6.py --stage 1 --session-min 118 --usage-cap 0.80` (4 sessions) | 0, 0, 0, 0 |
| `python benchmark/wp6/harness/evaluate_wp6.py --selftest` | 0 (exact p 0.052 = brute force 0.052) |
| `python benchmark/wp6/harness/evaluate_wp6.py --stage 1 > benchmark/wp6/results/report_tables_wp6_stage1.md` | 0 |
| `python benchmark/wp6/harness/sensitivity_wp6.py --stage 1` | 0. The final invocation reproduced 30 stored original audits and wrote 3. An interim invocation at 30 runs had written 19 of the 30; the other 11 are the pre-amendment originals. |
| `python benchmark/wp6/harness/explore_wp6.py --stage 1` | 0 |
| `python -m pytest -q benchmark/tests/test_wp6_harness.py` | 0 (31 passed) |
