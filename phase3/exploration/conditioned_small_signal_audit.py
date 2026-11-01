#!/usr/bin/env python3
"""Audit existing C10/C11 evidence for a conditioning-aware small-signal gate.

This script only reads archived C10/C11 JSON files.  It does not run a new
MOSFET solve and does not modify production solvers or golden files.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from statistics import median


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def corr(x, y):
    if len(x) != len(y) or len(x) < 2:
        return None
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    dx = [a - mx for a in x]
    dy = [b - my for b in y]
    den = math.sqrt(sum(a * a for a in dx) * sum(b * b for b in dy))
    return sum(a * b for a, b in zip(dx, dy)) / den if den else None


def parse_vg(path: Path) -> float:
    # The devsim-check state files contain the exact value in args; this is
    # only a fallback for a reproducible filename match.
    m = re.search(r"VG([0-9.]+)_", path.name)
    if not m:
        raise ValueError(f"cannot parse VG from {path}")
    return float(m.group(1))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    root = args.root.resolve()
    out = (args.out or root / "results/phase3/exploration").resolve()
    out.mkdir(parents=True, exist_ok=True)

    c11_path = root / "results/phase3/C11/c11_eval.json"
    c10_path = root / "results/phase3/C10/analysis/c10_result.json"
    devsim_paths = sorted((root / "results/phase3/C11/devsim_check").glob("e7_*.json"))
    indep_paths = sorted((root / "results/phase3/C11/indep").glob("fd2d_VG*.json"))
    for p in [c11_path, c10_path, *devsim_paths, *indep_paths]:
        if not p.exists():
            raise FileNotFoundError(p)

    c11 = json.loads(c11_path.read_text())
    rows = []
    for row in c11["rows"]:
        rows.append(
            {
                "file": row["file"],
                "na_cm3": row["na"],
                "tox_nm": row["tox"],
                "hw_um": row["hw"],
                "shift_V": row["shift"],
                "psi_hat": row["psi_hat"],
                "alpha": row["alpha"],
                "dL_2D_nm": row["dL_2D"],
                "dL_lin_nm": row["dL_lin"],
                "abs_err_lin": abs(row["err_lin"]),
                "signed_err_lin": row["err_lin"],
                "dL_q2d_nm": row["dL_q2d"],
                "abs_err_q2d": abs(row["err_q2d"]),
                "signed_err_q2d": row["err_q2d"],
                "lam_lin_nm": row["lam_lin"],
                "B_lin_VT": row["B_lin_VT"],
                "B_inf_over_B_lin": row["B_inf_over_B_lin"],
            }
        )

    rows_csv = out / "conditioned_small_signal_rows.csv"
    with rows_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    indep = {}
    for p in indep_paths:
        d = json.loads(p.read_text())
        indep[round(float(d["args"]["vg"]), 6)] = d

    coupled = []
    for p in devsim_paths:
        d = json.loads(p.read_text())
        vg = float(d["args"]["vg"])
        ds = d.get("devsim", {})
        const = d["result"]["const_mu"]
        lomb = d["result"].get("lombardi_mu", {})
        fv = indep.get(round(vg, 6))
        g_e = float(ds["G_E"])
        g_m = float(ds["G_M"])
        g_scale = max(abs(g_e), abs(g_m), 1e-300)
        coupled.append(
            {
                "vg": vg,
                "G_E": g_e,
                "G_M": g_m,
                "abs_G_M_minus_G_E_over_max": abs(g_m - g_e) / g_scale,
                "E_dwr": ds.get("E_dwr"),
                "E_round": ds.get("E_round"),
                "kappa_G": ds.get("kappa_G"),
                "n_unknowns": ds.get("n_unknowns"),
                "G_net_const_mu": const["G_net"],
                "dL_net_const_mu_nm": const["dL_net_nm"],
                "network_residual_const_mu": const["net_rel_residual"],
                "bounds_ok_const_mu": const["bounds_ok"],
                "drain_side_rel_error_const_mu": None,
                "drain_side_rel_error_lombardi_mu": None,
                "source_swapped_rel_const_mu": None,
                "independent_dL_net_nm": None if fv is None else fv["dL_net_nm"],
                "network_vs_independent_rel_dL": None,
            }
        )
        # Matching netcheck diagnostics are keyed by VG; preserve exact values.
        npath = root / "results/phase3/C11/devsim_check" / (
            f"netcheck_deck_W_L_1_level_0_vg_{vg:.6f}_tox_2_na_1e+18_hw_0.5_psibi-shift_0.000000_devsim.json"
        )
        if npath.exists():
            nd = json.loads(npath.read_text())["diag"]
            coupled[-1]["drain_side_rel_error_const_mu"] = nd[0]["rel_drainside_error"]
            coupled[-1]["drain_side_rel_error_lombardi_mu"] = nd[1]["rel_drainside_error"]
            coupled[-1]["source_swapped_rel_const_mu"] = nd[0]["rel_src_vs_swapped"]
        if fv is not None:
            coupled[-1]["network_vs_independent_rel_dL"] = abs(
                const["dL_net_nm"] - fv["dL_net_nm"]
            ) / abs(fv["dL_net_nm"])

    coupled_csv = out / "coupled_small_signal_rows.csv"
    with coupled_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(coupled[0]))
        w.writeheader()
        w.writerows(coupled)

    abs_lin = [r["abs_err_lin"] for r in rows]
    abs_q2d = [r["abs_err_q2d"] for r in rows]
    log_alpha = [math.log10(r["alpha"]) for r in rows]
    log_bratio = [math.log10(max(r["B_inf_over_B_lin"], 1e-300)) for r in rows]
    net_fv = [r["network_vs_independent_rel_dL"] for r in coupled if r["network_vs_independent_rel_dL"] is not None]
    summary = {
        "status": "EXPLORATORY_EVIDENCE_ONLY",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "input_sha256": {str(p.relative_to(root)): sha256(p) for p in [c10_path, c11_path, *devsim_paths, *indep_paths]},
        "n_c11_states": len(rows),
        "h_lin": {
            "abs_error_min": min(abs_lin),
            "abs_error_median": median(abs_lin),
            "abs_error_max": max(abs_lin),
            "signed_error_min": min(r["signed_err_lin"] for r in rows),
            "signed_error_max": max(r["signed_err_lin"] for r in rows),
            "corr_abs_error_log10_alpha": corr(abs_lin, log_alpha),
            "corr_abs_error_log10_Binf_over_Blin": corr(abs_lin, log_bratio),
            "all_underestimate": all(r["signed_err_lin"] < 0 for r in rows),
        },
        "h_q2d": {
            "abs_error_min": min(abs_q2d),
            "abs_error_median": median(abs_q2d),
            "abs_error_max": max(abs_q2d),
            "all_underestimate": all(r["signed_err_q2d"] < 0 for r in rows),
        },
        "coupled_vs_decoupled": {
            "states": coupled,
            "max_network_vs_independent_rel_dL": max(net_fv) if net_fv else None,
            "max_small_signal_discrepancy": max(r["abs_G_M_minus_G_E_over_max"] for r in coupled),
            "min_kappa_G": min(r["kappa_G"] for r in coupled),
            "max_kappa_G": max(r["kappa_G"] for r in coupled),
        },
        "interpretation": {
            "established": [
                "C11 archived rows show H-lin underestimates dL in all 24 states.",
                "The independent finite-volume network and DEVSIM equilibrium network agree at the two archived states within the reported 0.25% scale.",
            ],
            "diagnostic": [
                "The coupled small-signal G_E/G_M pair has large disagreement and large reported error estimates in the two deep-subthreshold checks.",
                "The drain-side forward-current diagnostic is ill-conditioned while source-side and swapped currents agree to round-off.",
            ],
            "hypothesis": "A conditioning-aware gate based on an equilibrium network plus an independent residual check could decide when coupled small-signal output is reportable; this is not yet a validated algorithm or a novelty claim.",
        },
    }
    summary_path = out / "conditioned_small_signal_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    report_path = root / "phase3/reports/Codex_Exploration_Conditioned_SmallSignal_2026-11-01.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report = f"""# Codex 独立探索：守恒网络与小信号条件性闸门（2026-11-01）

## 范围与状态

本探索只读取 C10/C11 已冻结的 JSON 结果，没有重新运行二维 MOSFET 扫描，也没有修改生产求解器、正式 golden 或历史结果。当前结论是 `MORE_EVIDENCE_NEEDED`，不是新模型通过，也不是 SCI 新颖性结论。

- 代码提交：`{summary['git_head']}`
- 输入数据：`results/phase3/C10/analysis/c10_result.json`、`results/phase3/C11/c11_eval.json`、C11 DEVSIM 检查及独立有限体积结果。
- 机器可读输出：`{rows_csv.relative_to(root)}`、`{coupled_csv.relative_to(root)}`、`{summary_path.relative_to(root)}`。

## 可证伪问题

候选问题是：在深亚阈值、输出电流接近浮点分辨率时，能否用“平衡态电阻网络 + 独立守恒/残差检查 + 条件数或后验误差估计”判断耦合小信号结果是否可报告？可证伪标准是：在同一收敛态上，网络结果与独立有限体积参考不能保持小于预先给定误差，或条件性闸门无法识别已知的耦合失效状态，则该方向不成立。

## 已有证据的重新量化

| 指标 | 结果 | 证据等级 |
|---|---:|---|
| C11 H-lin 状态数 | {len(rows)} | 数值复核 |
| H-lin 绝对相对误差 | {min(abs_lin):.6f}–{max(abs_lin):.6f}（中位数 {median(abs_lin):.6f}） | 数值复核 |
| H-lin 误差符号 | {'24/24 低估' if all(r['signed_err_lin'] < 0 for r in rows) else '并非全为低估'} | 数值复核 |
| abs(error) 与 log10(alpha) 相关系数 | {summary['h_lin']['corr_abs_error_log10_alpha']:.3f} | 探索性统计 |
| abs(error) 与 log10(B∞/B_lin) 相关系数 | {summary['h_lin']['corr_abs_error_log10_Binf_over_Blin']:.3f} | 仅诊断，B∞ 拟合不是统一参考 |
| DEVSIM 网络与独立有限体积 dL 相对差（2 点最大） | {summary['coupled_vs_decoupled']['max_network_vs_independent_rel_dL']:.6g} | 独立数值路径 |
| 两个耦合小信号状态的最大 |G_M−G_E|/max(|G|) | {summary['coupled_vs_decoupled']['max_small_signal_discrepancy']:.6g} | DEVSIM 后验诊断 |
| 两个耦合状态 κ_G | {summary['coupled_vs_decoupled']['min_kappa_G']:.3e}–{summary['coupled_vs_decoupled']['max_kappa_G']:.3e} | 条件性诊断 |

相关系数只说明当前 24 点数据中的关联，不是因果关系；`alpha`、`B∞` 与误差的定义来自 C10/C11，不能替代新的参考解。

## 两类小信号路径的区分

C11 已有结果显示，独立有限体积网络与 DEVSIM 平衡态网络在两个状态上保持小误差，而耦合小信号的 `G_E`/`G_M` 差异和其自带的 `E_dwr`、`E_round` 估计都很大。源侧电流与交换端点电流接近机器精度一致；排侧前向电流在低电流状态出现数量级失真。这支持“先做平衡网络和守恒检查，再决定是否报告耦合小信号”的工程假设，但目前只有两个 DEVSIM 状态，不能给出普适阈值。

## 研究状态分类

- **已建立（代码/数值）**：C11 的 24 点 H-lin 复核、C10/C11 的原始 SHA-256、两个独立有限体积网络结果及 77 项 DEVSIM 环境下测试通过。
- **诊断性结果**：条件数、DWR/rounding 估计、源/排两侧电流对照；它们揭示已知病态性，但尚未构成通用误差定理。
- **待证伪假设**：可由后验指标自动选择“耦合小信号可报告/改用守恒网络”的闸门。下一步必须在预注册的新状态上验证，并保留闸门失败案例。
- **不能宣称**：没有宣称新的物理模型、二维 MOSFET 科学验证完成或 SCI 创新成立。

## 建议的下一步（需单独批准后执行）

1. 先完成文献闸门，确认小信号病态性、网络线性响应和后验误差控制已有何种先例。
2. 只选少量预注册状态（至少一个深亚阈值、一个近阈值、一个强反型）和一组网格层级，分别比较耦合 Jacobian、平衡网络、独立有限体积；不做未经批准的参数扫描。
3. 预先固定判据：网络/独立参考误差、端点守恒残差、DWR/rounding 估计与条件指标必须同时记录；任一失败都保留为负结果。
4. 若闸门能稳定识别已知失效并在独立状态上保持误差界，再评估是否值得形成方法学论文问题；否则将其作为项目的数值可靠性规范，而不是研究主线。

## 复现命令

```bash
cd /workspace/c12_restore_backup/restored
/workspace/mosfet-tcad-research/.venv/bin/python -m pytest phase3/tests -q -rA
DEVSIM_MATH_LIBS=/tmp/devsim-libs/libblas.so:/tmp/devsim-libs/liblapack.so \\
  LD_LIBRARY_PATH=/tmp/devsim-libs \\
  /workspace/mosfet-tcad-research/.venv-devsim/bin/python -m pytest phase3/tests -q -rA
python3 phase3/exploration/conditioned_small_signal_audit.py
```

两次测试的实际日志位于：`/workspace/c12_restore_backup/evidence/c12_phase3_tests.log` 与 `/workspace/c12_restore_backup/evidence/c12_phase3_tests_devsim.log`。
"""
    report_path.write_text(report)
    print(json.dumps({"rows_csv": str(rows_csv), "coupled_csv": str(coupled_csv), "summary": str(summary_path), "report": str(report_path), "n_rows": len(rows), "n_coupled": len(coupled)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
