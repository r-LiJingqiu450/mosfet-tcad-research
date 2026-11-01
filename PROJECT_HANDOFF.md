# Active research handoff — 2026-11-01

- Branch: `phase3/direction-reset-2026-11-01`, based on `fde679b561bddf35663933890d7ace617bf9b928`.
- Active question: can a finite n+ well replace Ward's semi-infinite well using a reproducible `h_w/W_dep` saturation criterion?
- Decision checkpoint: the preregistered `h_w/W_dep >= 3` and 1% criterion fails at `h_w=1.0 um`: the 0.50→1.00 um `Delta L_net` change is 1.47115%. The result rejects that specific replacement rule; it does not establish a universal correction law.
- Evidence: `phase3/reports/Direction_Reset_Well_Depth_Report_2026-11-01.md`; raw data under `results/phase3/direction_reset/well_depth/`, DEVSIM cross-check under `results/phase3/direction_reset/well_depth_devsim/`, logs under `logs/phase3/direction_reset/`.
- Independent numerical check: at `h_w=1.0 um`, network vs DEVSIM conductance differences are `1.57e-7` and `3.34e-7`; the condition-number proxy is `3.46e10`, so this is numerical-path consistency under a shared model, not independent physical validation.
- Old finite-bias/two-length/Ward route: preserved as historical evidence, paused as the default topic. Do not resume it without a materially new question.
- Next investment gate: only test the finite-depth correction idea after level refinement, at least one independent finite-volume/discretization path, and multiple `N_A,V_G` conditions can falsify or support a predictive scaling relation. Otherwise retain the negative limitation and select a different direction.
