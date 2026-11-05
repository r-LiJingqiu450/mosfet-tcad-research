# WP6 execution notes (dated, append-only)

These notes record execution events that the machine logs (queue logs, ledger, run records) do not explain by themselves.
They are summarized in the WP6 report under deviations.

- **2026-11-05 13:34–13:38 UTC: usage cap overshot (W6-012, W6-013).**
  - The stage-1 queue launched W6-012 and W6-013 at a last-reported 5-hour utilization of 0.75, below the 0.80 cap (W6-A1).
  - The actual value at launch was already 0.91: the first usage events inside both runs report 0.91.
  - The jump from 0.75 (13:28) to 0.91 (13:34) came from the coordinating session (its context compaction and work), not from agent runs. The queue sees usage only through values reported inside run transcripts, so it had no reading in between.
  - Both runs finished normally (exit 0; outcome recorded by the registered rules). The queue then stopped at 0.93.
  - The cap is an execution limit agreed with the user, not a validity criterion, so no run is excluded.
  - Consequence: from here on the coordinating session keeps its own work light while runs are launched, and resumes the queue only after the 5-hour window resets (16:20 UTC).

- **2026-11-05 16:24 UTC: container restart before queue session 3.** The cloud machine had restarted while idle (uptime 0 min at 16:24). Checked before relaunching:
  - guard installed (`/opt/wp6/guard.py`), Claude Code 2.1.318, the Python environment, 2 CPUs;
  - the timezone was America/Denver both before (file times in the agents' `ls -la` output of S1 to W6-013 are UTC−7) and after.
  - Unchanged environment otherwise as far as the run records show (same tools, model, permission mode and session settings in every init event).
- **2026-11-05 17:04 UTC: transient permission-classifier outage in W6-024.** Nine tool calls were refused with "Classifier unavailable" (Claude Code's auto-mode classifier gave no verdict). The agent repeated them and finished (CS). Recorded per run; not a harness failure under the registered F4 definition (tools, MCP servers, model).
- **2026-11-05 18:06 UTC: stage 1 complete** (33 of 33 runs; queue sessions 1–4). Registered evaluation: 0 silent errors among 33 valid arm-A runs, so the stop rule (P-A) applies: stage 2 (arm B) is not run. Halt rule: 0 of 33 runs excluded.
