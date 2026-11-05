# WP6: AI-agent study (silent errors in low-current conductance)

Pre-registration: `benchmark/specs/WP6_Preregistration_Agent_Study_2026-11-05.md`.

| Path | Role |
|---|---|
| `agent/guard.py` | PreToolUse hook of the agent under test: refuses tool calls outside its workspace (installed to `/opt/wp6/guard.py` by `harness/install_guard.py`) |
| `kit/` | what an agent sees: `mosfet_deck.py` (DEVSIM deck), `py` (Python wrapper), `check` (arm-B client), `docs/` (DEVSIM 2.11.2 manual and examples, Apache-2.0) |
| `tasks/T01..T11, SMOKE/` | per task: `device.json` (neutral description), `mesh.json` (published mesh), `zero_bias_state.npz` (converged WP2 arrays) |
| `tasks/tasks.json` | private task table: case, state, mesh, reference G_DD (never copied into a workspace) |
| `harness/make_kit_states.py` | regenerates `tasks/` from the WP2 deck and checks each reference bit for bit |
| `harness/validate_kit.py` | checks that the agent-facing deck reproduces the state and the reference; optional drift–diffusion check |
| `harness/prompts.py` | the fixed system and task prompts |
| `harness/make_schedule.py` → `harness/schedule.json` | seeded run order: smoke S1/S2, stage 1 (arm A) W6-001..033, stage 2 (arm B) W6-034..066 |
| `harness/run_agent.py` | one run: workspace, checker daemon (arm B), headless `claude -p` session, archive |
| `harness/checker_daemon.py` | arm-B verification service (PASS/FAIL within 1%, at most 3 checks) |
| `harness/run_queue_wp6.py` | runs a stage two at a time within a session and a usage cap; audits each run |
| `harness/audit.py` | mechanical transcript audit (contamination flags F1–F4, final block, usage) |
| `harness/evaluate_wp6.py` | registered evaluation (outcomes, rates, exact stratified test, stop and halt rules) |
| `results/` | ledger, `runs/<id>/` (run record, gzip transcript, audit, checker log, files the agent wrote), evaluations |

Reproduce (repository root, Python environment with DEVSIM 2.11.2):

```bash
python benchmark/wp6/harness/make_kit_states.py      # tasks/ from the WP2 deck
python benchmark/wp6/harness/validate_kit.py --dd T02
python benchmark/wp6/harness/install_guard.py
python benchmark/wp6/harness/run_queue_wp6.py --stage smoke
python benchmark/wp6/harness/run_queue_wp6.py --stage 1     # repeat until the stage is complete
python benchmark/wp6/harness/evaluate_wp6.py --stage 1
```
