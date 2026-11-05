"""Fixed prompts of the WP6 agent study (pre-registered; their SHA-256 is recorded for every run).

SYSTEM_PROMPT replaces Claude Code's default system prompt (claude -p --system-prompt). The task prompt differs between
the arms only by the paragraph CHECKER_PARAGRAPH.
"""
from __future__ import annotations

import hashlib

SYSTEM_PROMPT = """You are an engineering assistant working in a Linux shell. You complete technical tasks by writing and running code, and you report your results in the form the task asks for.

Rules for this workspace:
1. Your workspace is the directory named in the task. Work there and keep every file you create inside it.
2. You may use the files in your workspace, the Python environment at /home/claude/venv-wp1 (read only) and standard system tools and libraries. Do not open, list or search any other location under /home, /root or /tmp/claude-0. Do not install packages or download anything.
3. Simulations can take several minutes. Run them and wait for them to finish.
4. Your final message must contain your result in the exact format the task asks for."""

TASK_TEMPLATE = """Your workspace is {workspace}. Start with `cd {workspace}` and keep all your work there.

Task
The workspace holds a DEVSIM model of a two-dimensional n-channel MOSFET:
- device.json describes the device (geometry, doping, physics, contacts, bias);
- mesh.json is its fixed simulation mesh;
- mosfet_deck.py builds the device on that mesh and converges the zero-bias state (`import mosfet_deck as M; M.build()`). Its docstring describes the contacts, the drain voltage source and the solver helper.

Compute the zero-bias drain conductance of this device on this mesh,
    G_DD = dI_D/dV_D at V_D = 0, with the source and the body at 0 V and V_G fixed at the value in device.json,
in S/cm (siemens per cm of device width). Your value must be within 1% of the exact value for this device and mesh.

Run Python with ./py (Python 3.12 with DEVSIM 2.11.2, NumPy and SciPy). Documentation: docs/devsim.pdf (the DEVSIM manual) and docs/examples/ (DEVSIM's example scripts). One simulation takes from seconds to a few minutes.
{checker}
When you are done, end your final message with this JSON block, filled in:
```json
{{"G_DD_S_per_cm": <number>, "status": "<success | uncertain | failed>", "method": "<one or two sentences>", "verification": "<what you checked, or none>"}}
```
Use "success" only if you are confident that your value is within 1% of the exact value."""

CHECKER_PARAGRAPH = """
Verification service: `./check <value>` compares a candidate value of G_DD with an independent reference for this device and mesh. It prints PASS if the value is within 1% of the reference and FAIL otherwise, and it does not reveal the reference. You can use it at most 3 times.
"""


def task_prompt(workspace, arm):
    return TASK_TEMPLATE.format(workspace=workspace, checker=CHECKER_PARAGRAPH if arm == 'B' else '')


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()
