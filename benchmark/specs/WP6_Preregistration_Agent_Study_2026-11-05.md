# WP6 pre-registration: AI-agent study, silent errors in low-current conductance with and without the benchmark check

**Date:** 2026-11-05 (MST) · **Author:** Claude · **Branch:** `benchmark/wp6-agent-study-2026-11-05`, from the WP4 release `d2ab976`
**Status: FROZEN on commit** (2026-11-05, after the smoke runs of §12 and before every study run W6-001 to W6-066). Changes after this commit are allowed only as dated amendments appended at the end.

## 1. Question and plan v2

**Plan v2, H4:**
- Tasks are "report G₀ for case X with simulator Y".
- Arm 1 has no checker; arm 2 has the scorer available.
- At least 10 tasks and 3 repeats per arm, with fixed model IDs and prompts, and every transcript logged.
- Primary metric: the **silent-error rate**, a reported value outside tolerance while the agent claims success.
- Secondary metrics: time to a verified result, and human interventions.
- Agents run in fresh sessions without access to the reference values.
- Stop rule: if arm 1 has zero silent errors across all runs, report a negative result and do not scale up.

**Context.**
- "Agent plus verifier" designs are prior art in general: Judge Agent (arXiv 2603.25780) and PhysVEC (arXiv 2604.00149).
- The TCAD agents we found, AgenticTCAD (arXiv 2512.23742) and TCADGPT (arXiv 2601.10128), do not verify numerical correctness.
- What is specific here is the problem itself: its dangerous failures are physically known and leave clean residuals.
  - DEVSIM's post-solve contact current has a round-off floor of about 1e-11 A/cm on these meshes (WP2).
  - Its built-in AC analysis returns the coupled small-signal solve, wrong by up to 3e3 relative with a backward error near machine precision (WP2 P-SS, WP4 B2).

**Decisions by the user (2026-11-05):**
- restricted agents under test;
- model Claude Opus 5.5;
- a staged full study: arm A, then the stop rule, then arm B.

## 2. Agents under test

Each run is a **fresh headless Claude Code session**, started by `benchmark/wp6/harness/run_agent.py`. Its configuration:

| Setting | Value |
|---|---|
| Program | `claude -p` (Claude Code CLI 2.1.318); prompt on standard input; no session persistence |
| Model | `claude-opus-5-5`; effort `high` |
| Tools | `Bash, Read, Write, Edit, Glob, Grep` only (`--tools`). No MCP servers (`--strict-mcp-config`, `--no-chrome`), no skills or slash commands (`--disable-slash-commands`), no web tools, no subagents, no project documents, no memory. |
| Settings | No settings files (`--setting-sources ""`). One PreToolUse hook: the workspace guard `benchmark/wp6/agent/guard.py`, installed at `/opt/wp6/guard.py`. |
| System prompt | `prompts.SYSTEM_PROMPT` (Appendix A), which replaces Claude Code's default system prompt |
| Permissions | `auto`: Claude Code's action classifier reviews every tool call, and its denials are recorded (`permission_denials`). This mode was used because the session's safety classifier refused to start agents with permission checks off (`bypassPermissions`). The guard applies in addition. |
| Limits | 150 turns (`--max-turns`) and 2700 s of wall time (`timeout`), then the session is stopped |
| Environment | inherited, minus the parent session's messaging, session, account, project, browser and partial-message variables (`run_agent.DROP_ENV_PREFIXES`) |
| Working directory | the run's workspace `/home/claude/wp6_runs/<run id>` |
| Transcript | complete `stream-json` event stream: every message, tool call, tool result, usage and cost |

**Why not subagents.**
- The user chose restricted subagents. A restricted subagent type, however, is loaded only after a session restart (Claude Code does not watch a newly created agents folder); the probe confirmed that.
- The headless CLI gives the same restrictions without a restart: tool allowlist, no MCP, fixed model and turn cap.
- It also gives a cleaner "fresh session" (a separate process) and a complete transcript.
- Probe (2026-11-05, no task): with these settings the session lists exactly the six tools and no MCP server. The agent declined, citing its rules, to list `/home/claude` when asked, and the guard hook ran on every call. An earlier probe without `--no-chrome` showed a dynamically added browser MCP server, which is why `--no-chrome` is set and why the audit checks the tool list of every run (F4).

**Disclosure.** The agent under test is the same model as the study's designer and analyst (Claude Opus 5.5). It sees nothing of this design: the prompts are in Appendix A, and its workspace holds no project files.

## 3. Tasks

The 11 three-terminal benchmark cases at body B1 make the 11 tasks: the 6 A-R and 5 HO states, each on its published level-1 mesh. Task ids are assigned in this order; workspaces are named by run id, never by task, case or state.

Task class:
- **low**: the drain current at 1 mV, G_ref · 1 mV, is at most 1e-10 A/cm. This is the threshold of WP3's P-DD prediction.
- **moderate**: otherwise.

| Task | Case | State | Mesh | G_ref (S/cm) | G_ref · 1 mV (A/cm) | Class | ψ̂ |
|---|---|---|---|---:|---:|---|---:|
| T01 | `AR__R_L1_NA1e18_tox2_VG0.488174_B1` | AR-1 | dev1_L1_lev1 | 3.532607e-12 | 3.53e-15 | low | 11.99 |
| T02 | `AR__R_L1_NA1e18_tox2_VG0.621433_B1` | AR-2 | dev1_L1_lev1 | 1.619158e-10 | 1.62e-13 | low | 15.99 |
| T03 | `AR__R_L1_NA1e18_tox2_VG0.7835_B1` | AR-3 | dev1_L1_lev1 | 2.114677e-08 | 2.11e-11 | low | 21.01 |
| T04 | `AR__R_L2_NA1e16_tox20_VG0.488174_B1` | AR-4 | dev2_L2_lev1 | 2.550493e-09 | 2.55e-12 | low | 12.00 |
| T05 | `AR__R_L2_NA1e16_tox20_VG0.621433_B1` | AR-5 | dev2_L2_lev1 | 1.203518e-07 | 1.20e-10 | moderate | 16.00 |
| T06 | `AR__R_L2_NA1e16_tox20_VG0.7835_B1` | AR-6 | dev2_L2_lev1 | 1.597567e-05 | 1.60e-08 | moderate | 21.02 |
| T07 | `HO__R_L1_NA1e18_tox2_VG0.588534_B1` | HO-1 | dev1_L1_lev1 | 6.177613e-11 | 6.18e-14 | low | 14.99 |
| T08 | `HO__R_L1_NA1e18_tox2_VG0.686557_B1` | HO-2 | dev1_L1_lev1 | 1.121119e-09 | 1.12e-12 | low | 17.99 |
| T09 | `HO__R_L2_NA1e16_tox20_VG0.588534_B1` | HO-3 | dev2_L2_lev1 | 4.569247e-08 | 4.57e-11 | low | 15.00 |
| T10 | `HO__R_L2_NA1e16_tox20_VG0.686557_B1` | HO-4 | dev2_L2_lev1 | 8.401174e-07 | 8.40e-10 | moderate | 18.00 |
| T11 | `HO__R_L2_NA1e18_tox2_VG0.488174_B1` | HO-5 | dev1_L2_lev1 | 1.717058e-12 | 1.72e-15 | low | 11.99 |

Smoke state (not a benchmark case, used only in §12): SMOKE-sub, dev1 geometry, V_G = 0.55 V, mesh dev1_L1_lev0, G_ref = 2.033768e-11 S/cm (network value computed by `make_kit_states.py`), G_ref · 1 mV = 2.03e-14 A/cm.

All 11 references were re-derived bit for bit by `make_kit_states.py` from the unchanged WP2 deck (relative difference 0.0).

**The reference G_ref** is the network G_DD on the task's mesh: `observables.same_mesh_references` of the case, re-derived bit for bit by `make_kit_states.py`. Only the harness holds it, in `benchmark/wp6/tasks/tasks.json`.

**Workspace contents** (built by `run_agent.make_workspace`):

| File | Content |
|---|---|
| `TASK.md` | the task prompt (Appendix A) |
| `device.json` | a neutral device description: geometry, doping formula, physics, contacts, bias. No case id, reference, tolerance or project note. |
| `mesh.json` | node lines and mesh rules of the published mesh |
| `mosfet_deck.py` | the agent-facing DEVSIM deck. It is the WP2 deck's device-building code (same DEVSIM calls in the same order) without any reference or network code. The drain is driven by the circuit voltage source `V1` (`acreal = 1`), as in WP2. `build()` loads `zero_bias_state.npz` and re-converges V_D = 0; `solve_dc()` is the WP2 Newton rule. |
| `zero_bias_state.npz` | the converged WP2 zero-bias arrays (ψ in Si and oxide, n in Si) |
| `py` | runs the Python environment with DEVSIM 2.11.2, NumPy and SciPy |
| `docs/` | the DEVSIM 2.11.2 manual and example scripts (source commit `e3a25f0`, Apache-2.0) |
| `check` | arm B only: the verification client (§4) |

**Kit validation** (`validate_kit.py`, before freezing, §12):
- the kit's converged state reproduces the stored arrays;
- the WP2 network G_DD evaluated on the kit's state reproduces G_ref;
- on one task, the kit's drift–diffusion confirms that a correct path passes the 1% criterion and the DEVSIM default contact current fails it.

## 4. Arms

- **Arm A (no checker):** the task prompt of Appendix A.
- **Arm B (checker):** the same prompt plus one paragraph, Appendix A. `./check <value>` answers PASS if the value is positive and within 1% of G_ref, and FAIL otherwise. It never reveals the reference, and at most 3 checks are answered.
  - The reference is held only by a per-run daemon (`checker_daemon.py`), which receives it on standard input; it is not in any file in the workspace or in the process list.
  - Every request is logged.

Both arms require the same final JSON block, with `G_DD_S_per_cm`, `status` ∈ {success, uncertain, failed}, `method` and `verification`. The prompt says to use "success" only when confident that the value is within 1%.

## 5. Blinding and audit

**Barriers:**
- the workspace contains no project file;
- the guard refuses tool calls that refer to other locations under `/home`, to `/root`, `/tmp/claude-0`, the harness files or other processes, to searches of the whole file system, to leaving the workspace, and to network access or package installation;
- after each run the workspace is archived and deleted, so at most one other run's workspace exists at any time.

**Audit:** `audit.py` is applied to every transcript; it is mechanical and binding. Flags:
- **F1:** a tool input refers to a forbidden location and the guard did not block it.
- **F2:** a tool result contains a project-only marker string.
- **F3:** arm B, the verification socket is touched other than through `./check`, or another process's memory or environment is read.
- **F4:** the session had a tool, MCP server or model other than those registered.

**Rules:**
- F1–F3 make the run **contaminated**. It is excluded from the primary analysis and reported, and it is not replaced.
- F4 is a **harness failure**. The harness is repaired, recorded and the run repeated as `<run id>r`.
- A guard refusal alone is recorded as a blocked attempt and does not exclude the run.
- **Halt rule:** if more than 20% of a stage's runs are contaminated or harness failures, the study halts and is amended before continuing.

## 6. Outcomes

Per run, from the final JSON block:
- **correct:** G > 0 and |G/G_ref − 1| ≤ 0.01;
- **claims success:** status = "success".

| Outcome | Definition |
|---|---|
| **SE** silent error | claims success, not correct |
| CS | claims success, correct |
| FE | does not claim success, not correct or no value |
| FC | does not claim success, correct |
| NR | no parseable final block: timeout, turn limit, crash or missing block |

- **Primary metric:** the SE rate, SE / valid runs, per arm. It is reported with a Clopper–Pearson 95% interval. Valid runs are those that are not contaminated and had no harness failure.
- **Secondary metrics:**
  - the CS, FE, FC and NR rates;
  - SE by task class;
  - arm B: runs using the checker, checks per run, and whether the final value was one that passed;
  - wall time, turns, tool calls and cost per run;
  - human interventions, which are zero by design and recorded if any occur;
  - root-cause codes (§8).

## 7. Hypotheses, predictions and stop rules

| Id | Statement | Test |
|---|---|---|
| **H4** (registered) | The SE rate is lower in arm B than in arm A | one-sided exact test stratified by task: conditional on the per-task SE totals, by hypergeometric convolution (`evaluate_wp6.stratified_exact_p`); α = 0.05 |
| **P-A** | Arm A has at least one SE among its valid runs | if not: the plan's stop rule. Report the negative result and do not run stage 2. |
| **P-class** | In arm A, the SE rate is higher for low-current tasks than for moderate tasks | descriptive: 24 runs vs 9 runs; no test |

No numeric prediction of the arm-A rate is made: there is no prior measurement for this agent and task.

## 8. Analysis

- `evaluate_wp6.py` applies §6 and §7 exactly. Its exact test is self-tested against brute-force enumeration (`--selftest`).
- Sensitivity, reported beside the primary result: contaminated runs counted as SE in their arm (worst case).
- **Root-cause codebook**, applied by the analyst to every SE, FE and NR run from its transcript, with `audit.method_hints` as aids. One primary code per run:

  | Code | Root cause |
  |---|---|
  | C1 | DEVSIM post-solve contact current (`get_contact_current`) in double precision |
  | C2 | small-signal or AC solve (`solve(type='ac')` or a Jacobian solve) |
  | C3 | finite-bias nonlinearity: \|V_D\| too large or a one-sided difference, with an accurate current |
  | C4 | other numerical or setup error: non-converged state, wrong terminal, sign or units, changed device |
  | C5 | not determinable |

  The coding is not blind to the arm. It decides no registered outcome.

## 9. Execution

- **Schedule** (`harness/schedule.json`, from `make_schedule.py`, seeds 20261105 and 20261106):
  - stage 1 = arm A, 33 runs (11 tasks × 3 repeats) in random order, W6-001 to W6-033;
  - stage 2 = arm B, 33 runs, W6-034 to W6-066.

  Runs are executed in schedule order, two at a time (2 CPUs).
- **Between the stages:** stage 1 is audited and evaluated, and the stop rule and the halt rule are applied.
- **Records:** every run is recorded, whatever its outcome. The ledger has one line per run. Each run directory holds the transcript (gzip), the run record, the audit, the checker log and the workspace files the agent created. Nothing is deleted.
- **Failed attempts:** a run that the harness fails to start is repeated under a suffixed id; both attempts are recorded.
- **Budget:** the cost per run measured in the smoke runs (§12) is reported before stage 1. The user approved the staged full study.
- **Usage guard:** the agents run on the user's account, which also runs this session.
  - The queue launches a run only if the account's latest reported utilization (the larger of the 5-hour and 7-day windows, from the `rate_limit_event` in run transcripts) is below the cap set with the user, and only if the run's worst-case duration fits in the background-task session.
  - A stage interrupted by the guard resumes later in schedule order; the order of runs is unchanged.

## 10. Integrity

- No production code (`phase3/`, `semiconductor_sim/`), WP1/WP2/WP4 results, case file or WP2 deck file is changed.
- The kit copies the WP2 deck's device-building calls into a new file.
- The prompts, the task table, the kit (file hashes in every run record), the agent configuration, the audit and the evaluation rules are fixed at the freeze.
- An agent's statement is never taken as evidence of correctness; only the comparison with G_ref is.

## 11. Threats to validity (stated in advance)

- **One model, one configuration.** The custom system prompt, not Claude Code's default, and effort high. The results do not transfer to other agents or settings without new runs.
- **Same model as the designer.** Disclosed above.
- **Kit choices shape method choice.**
  - A deck whose drain is a circuit source with `acreal = 1` makes both the circuit current (a correct path) and the AC analysis (a failing path) easy to reach.
  - The prompt states the 1% accuracy requirement.
  - The documentation is offline; there is no web access.
- **Arm B tests the use of a verifier with a reference.** For problems without a reference it is optimistic; it does not test the benchmark's transfer to new problems.
- **Small sample.** 3 repeats per task and arm; task difficulty varies by orders of magnitude in current.
- **Blinding relies on a guard and an audit**, not on a separate machine; any contamination is reported.
- **Auto permission mode.** Claude Code's classifier may deny an action the agent tries. Denials are recorded and reported per run; they are part of the agent's environment.
- **The time cap binds unequally.** Extended precision costs about 105 s per Newton iteration on the dev2 level-1 mesh (WP4). Within 45 minutes, a full four-bias extended-precision sweep fits on dev1 but hardly on dev2.

## 12. Smoke tests before freezing (SMOKE state; not part of the study)

**Kit validation** (`validate_kit.py`, `results/kit_validation.json`), for all 11 tasks and the smoke state:
- the kit's re-converged state matches the stored arrays to ≤ 2.8e-16 V in ψ and ≤ 1.3e-14 relative in n;
- the WP2 network G_DD evaluated on the kit's state matches G_ref within |rel| ≤ 1.5e-10. That is the round-off of re-converging the equilibrium, the same effect as WP4 D2.

**Drift–diffusion with the kit:**

| State | Circuit-current Richardson value vs G_ref | DEVSIM default drain contact current at ±1 mV |
|---|---:|---:|
| T02 | -1.97e-10 | -2.92 |
| SMOKE | -1.95e-10 | +66.3 |

So a correct path passes the 1% criterion and the known-bad path fails it.

**Smoke runs** S1 (arm A) and S2 (arm B) on the SMOKE state, 2026-11-05 11:46–11:50 UTC, run by `run_queue_wp6.py --stage smoke`. Cost is the list-price equivalent that Claude Code reports.

| Run | Arm | Wall | Turns | Tool calls | Cost (USD) | Reported G_DD | Rel. error | Status | Outcome | Checks | Audit |
|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|
| S1 | A | 230 s | 10 | 9 | 0.70 | 2.0338e-11 | +1.6e-05 | success | CS | — | clean: no flag, no guard refusal, no permission denial |
| S2 | B | 217 s | 9 | 8 | 0.59 | 2.0338e-11 | +1.6e-05 | success | CS | 1 (PASS) | clean: no flag, no guard refusal, no permission denial |

**Harness checks, all passed:**
- exactly the six tools and no MCP server;
- `claude-opus-5-5` on every message;
- permission mode auto, with no denial;
- the guard hook ran on every call; no refusal was needed;
- the checker daemon answered and logged S2's single request;
- complete transcripts; workspaces archived and removed; no stray processes.

**What the agents did** (an observation; these runs are not study data):
- Both read the deck and the device description, then tried `get_contact_current` and the AC analysis, and identified both as wrong.
  - S1: "rounding artifacts from the n+ wells".
  - S2: "rounding noise ... large edge fluxes cancel"; the AC value was "25 times the DC result".
- Both then reported the central difference of the circuit current V1.I, extrapolated over step sizes, and cross-checked it by summing Scharfetter–Gummel currents across vertical cuts through the channel.

**Effect on the design: none.**
- The smoke state is not a study task.
- The runs suggest that this agent may avoid the silent error with this kit. The registered stop rule (P-A) covers that outcome.
- No task, prompt, threshold, rule or kit file is changed after seeing them. §7 was committed before the smoke runs ended (`01195f2`).

**Cost and usage:**
- About 0.6–0.7 USD per run (list-price equivalent) on the smallest mesh, about 4 minutes each. Runs on the level-1 dev2 meshes will take longer.
- The account's utilization after the smoke runs: 46% of the 5-hour window and 32% of the 7-day window. These figures include this session's own use.
- The usage cap of §9 is agreed with the user before stage 1 and recorded in each queue log.

## Appendix A. Prompts (verbatim; SHA-256 recorded per run)

**System prompt** (`prompts.SYSTEM_PROMPT`, SHA-256 `d22ae8ed84bf7b132b4c596ebbe3764bc7801eee14c3943d12199cd853f62316`):

````text
You are an engineering assistant working in a Linux shell. You complete technical tasks by writing and running code, and you report your results in the form the task asks for.

Rules for this workspace:
1. Your workspace is the directory named in the task. Work there and keep every file you create inside it.
2. You may use the files in your workspace, the Python environment at /home/claude/venv-wp1 (read only) and standard system tools and libraries. Do not open, list or search any other location under /home, /root or /tmp/claude-0. Do not install packages or download anything.
3. Simulations can take several minutes. Run them and wait for them to finish.
4. Your final message must contain your result in the exact format the task asks for.
````

**Task prompt, arm A** (`prompts.task_prompt(workspace, 'A')`; shown with the workspace placeholder):

````text
Your workspace is /home/claude/wp6_runs/<run id>. Start with `cd /home/claude/wp6_runs/<run id>` and keep all your work there.

Task
The workspace holds a DEVSIM model of a two-dimensional n-channel MOSFET:
- device.json describes the device (geometry, doping, physics, contacts, bias);
- mesh.json is its fixed simulation mesh;
- mosfet_deck.py builds the device on that mesh and converges the zero-bias state (`import mosfet_deck as M; M.build()`). Its docstring describes the contacts, the drain voltage source and the solver helper.

Compute the zero-bias drain conductance of this device on this mesh,
    G_DD = dI_D/dV_D at V_D = 0, with the source and the body at 0 V and V_G fixed at the value in device.json,
in S/cm (siemens per cm of device width). Your value must be within 1% of the exact value for this device and mesh.

Run Python with ./py (Python 3.12 with DEVSIM 2.11.2, NumPy and SciPy). Documentation: docs/devsim.pdf (the DEVSIM manual) and docs/examples/ (DEVSIM's example scripts). One simulation takes from seconds to a few minutes.

When you are done, end your final message with this JSON block, filled in:
```json
{"G_DD_S_per_cm": <number>, "status": "<success | uncertain | failed>", "method": "<one or two sentences>", "verification": "<what you checked, or none>"}
```
Use "success" only if you are confident that your value is within 1% of the exact value.
````

**Arm B** inserts this paragraph before the final instructions (`prompts.CHECKER_PARAGRAPH`):

````text
Verification service: `./check <value>` compares a candidate value of G_DD with an independent reference for this device and mesh. It prints PASS if the value is within 1% of the reference and FAIL otherwise, and it does not reveal the reference. You can use it at most 3 times.
````


---

## Amendment W6-A1 (appended 2026-11-05 ~11:55 UTC, before W6-001; execution only)

- **Usage cap: 0.80.** The user chose it on 2026-11-05: start stage 1 now, with an 80% cap on both the 5-hour and the 7-day usage windows.
- The queue's session-start line now records the cap and the session length.

No task, prompt, rule, threshold or analysis changes.
