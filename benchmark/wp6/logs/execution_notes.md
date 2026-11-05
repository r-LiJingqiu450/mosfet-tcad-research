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
