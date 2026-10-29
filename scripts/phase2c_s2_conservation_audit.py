#!/usr/bin/env python3
"""S2 flux and conservation audit on fixed converged Newton states."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
import warnings
from collections import Counter
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "semiconductor_sim"))

from common import constants as C  # noqa: E402
from common.bernoulli import bernoulli  # noqa: E402
from common.grid import make_junction_grid  # noqa: E402
from stage2_pn import solver as S  # noqa: E402
from stage2_pn.postprocess import current_density  # noqa: E402


TRACE = []


def _trace(frame, event, arg):
    if event == "call" and frame.f_code.co_name == "newton_solve":
        return _trace
    if event == "line" and frame.f_code.co_name == "newton_solve" and frame.f_lineno == 679:
        TRACE.append({"phi_n": frame.f_locals["phi_n"].copy(),
                      "phi_p": frame.f_locals["phi_p"].copy()})
    return _trace


def _stable_flux(x, phi, phi_n, phi_p):
    h = np.diff(x)
    d = (phi[:-1] - phi[1:]) / C.VT
    bp = bernoulli(d)
    bm = bernoulli(-d)
    fn = C.NI * np.exp(np.clip((phi[1:] - phi_n[:-1]) / C.VT, -745.0, 700.0))
    dn = (phi_n[:-1] - phi_n[1:]) / C.VT
    fp = C.NI * np.exp(np.clip((phi_p[:-1] - phi[1:]) / C.VT, -745.0, 700.0))
    dp = (phi_p[1:] - phi_p[:-1]) / C.VT
    jn = C.Q_E * C.D_N / h * bp * fn * np.expm1(dn)
    jp = -C.Q_E * C.D_P / h * bm * fp * np.expm1(dp)
    return jn, jp


def _decimal_b(z):
    return Decimal(1) if z == 0 else z / (Decimal(1) - (-z).exp())


def _decimal_flux(x, phi, n, p, phi_n, phi_p):
    """High precision SG evaluation from the captured QFP state and returned state."""
    with localcontext() as ctx:
        ctx.prec = 70
        q = Decimal(str(C.Q_E)); dn_const = Decimal(str(C.D_N)); dp_const = Decimal(str(C.D_P))
        vt = Decimal(str(C.VT)); ni = Decimal(str(C.NI))
        jn_qfp, jp_qfp, jn_ret, jp_ret = [], [], [], []
        for i, h_float in enumerate(np.diff(x)):
            h = Decimal(str(h_float))
            ph_i, ph_j = Decimal(str(phi[i])), Decimal(str(phi[i + 1]))
            d = (ph_i - ph_j) / vt
            bp, bm = _decimal_b(d), _decimal_b(-d)
            ni_i = ni * ((ph_i - Decimal(str(phi_n[i]))) / vt).exp()
            ni_j = ni * ((ph_j - Decimal(str(phi_n[i + 1]))) / vt).exp()
            pi_i = ni * ((Decimal(str(phi_p[i])) - ph_i) / vt).exp()
            pi_j = ni * ((Decimal(str(phi_p[i + 1])) - ph_j) / vt).exp()
            jn_qfp.append(q * dn_const / h * (bp * ni_j - bm * ni_i))
            jp_qfp.append(-q * dp_const / h * (bm * pi_j - bp * pi_i))
            nr_i, nr_j = Decimal(str(n[i])), Decimal(str(n[i + 1]))
            pr_i, pr_j = Decimal(str(p[i])), Decimal(str(p[i + 1]))
            jn_ret.append(q * dn_const / h * (bp * nr_j - bm * nr_i))
            jp_ret.append(-q * dp_const / h * (bm * pr_j - bp * pr_i))
        return (np.array([float(v) for v in jn_qfp]), np.array([float(v) for v in jp_qfp]),
                np.array([float(v) for v in jp_qfp]), np.array([float(v) for v in jp_ret]),
                np.array([float(v) for v in jn_ret]), np.array([float(v) for v in jp_ret]))


def _balance(x, jn, jp, U):
    h = np.diff(x)
    hi = 0.5 * (h[:-1] + h[1:])
    qU = C.Q_E * float(np.sum(U[1:-1] * hi))
    eb = jn[1:] - jn[:-1] - C.Q_E * hi * U[1:-1]
    pb = jp[1:] - jp[:-1] + C.Q_E * hi * U[1:-1]
    total = jn + jp
    scale = max(abs(qU), abs(float(np.median(total))), 1e-30)
    return {
        "qU_A_m2": float(qU),
        "Jn_delta_A_m2": float(jn[-1] - jn[0]),
        "Jp_delta_A_m2": float(jp[-1] - jp[0]),
        "Jn_over_qU": float((jn[-1] - jn[0]) / qU) if qU else None,
        "minus_Jp_over_qU": float(-(jp[-1] - jp[0]) / qU) if qU else None,
        "electron_local_inf_A_m2": float(np.max(np.abs(eb))),
        "hole_local_inf_A_m2": float(np.max(np.abs(pb))),
        "electron_local_relative": float(np.max(np.abs(eb)) / scale),
        "hole_local_relative": float(np.max(np.abs(pb)) / scale),
        "total_current_min_A_m2": float(np.min(total)),
        "total_current_max_A_m2": float(np.max(total)),
        "total_current_median_A_m2": float(np.median(total)),
        "total_current_spread_A_m2": float(np.ptp(total)),
        "total_current_spread_relative": float(np.ptp(total) / scale),
    }


def _qfp_residual(x, result, Va):
    Na = C.NA_CM3 * 1e6; Nd = C.ND_CM3 * 1e6
    bc = S.contact_values(Va, Na, Nd)
    phi = result["phi"]
    phin = result["_phi_n"]
    phip = result["_phi_p"]
    hw = np.diff(x); hi = 0.5 * (hw[:-1] + hw[1:])
    Na_arr = np.where(x >= 0.0, 0.0, Na); Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    F, *_ = S._qfp_assemble(phi, phin, phip, x, Na_arr, Nd_arr, hw, hi,
                             *S.contact_carriers(bc, Na, Nd))
    return float(np.linalg.norm(F)), float(np.max(np.abs(F)))


def _run_grid(name, config, selected=(-0.5, -1.0, -2.0)):
    x = make_junction_grid(**config)
    states = {}
    rows = []
    eq_state = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sys.settrace(_trace)
        try:
            eq = S.newton_solve(x, C.NA_CM3 * 1e6, C.ND_CM3 * 1e6, 0.0)
            eq_trace = TRACE[-1]
            eq = {**eq, "_phi_n": eq_trace["phi_n"], "_phi_p": eq_trace["phi_p"]}
            eq_state = eq
            state = eq
            for step in range(1, 21):
                va = -0.1 * step
                r = S.newton_solve(x, C.NA_CM3 * 1e6, C.ND_CM3 * 1e6, va,
                                   phi=state["phi"], n=state["n"], p=state["p"])
                t = TRACE[-1]
                state = {**r, "_phi_n": t["phi_n"], "_phi_p": t["phi_p"]}
                if round(va, 10) in {round(v, 10) for v in selected}:
                    states[round(va, 10)] = state
        finally:
            sys.settrace(None)
    warning_counts = Counter(str(w.message) for w in caught)
    for va, state in sorted(states.items()):
        J, Jt, jn_direct, jp_direct = current_density(state["phi"], state["n"], state["p"], x)
        jn_stable, jp_stable = _stable_flux(x, state["phi"], state["_phi_n"], state["_phi_p"])
        hp_jn_qfp, hp_jp_qfp, _, _, hp_jn_ret, hp_jp_ret = _decimal_flux(
            x, state["phi"], state["n"], state["p"], state["_phi_n"], state["_phi_p"])
        U = S.srh_rate(state["n"], state["p"])[0]
        direct = _balance(x, jn_direct, jp_direct, U)
        stable = _balance(x, jn_stable, jp_stable, U)
        high = _balance(x, hp_jn_qfp, hp_jp_qfp, U)
        qfp_l2, qfp_inf = _qfp_residual(x, state, va)
        rows.append({
            "grid": name, "N": len(x), "V_a": va, "converged": bool(state["converged"]),
            "iters": int(state["iters"]), "qfp_residual_l2": qfp_l2, "qfp_residual_inf": qfp_inf,
            "J_center_direct_A_m2": float(J),
            "direct": direct, "stable": stable, "high_precision_qfp": high,
            "max_direct_vs_stable_Jn": float(np.max(np.abs(jn_direct - jn_stable))),
            "max_direct_vs_stable_Jp": float(np.max(np.abs(jp_direct - jp_stable))),
            "max_stable_vs_high_Jn": float(np.max(np.abs(jn_stable - hp_jn_qfp))),
            "max_stable_vs_high_Jp": float(np.max(np.abs(jp_stable - hp_jp_qfp))),
            "max_direct_vs_high_returned_Jn": float(np.max(np.abs(jn_direct - hp_jn_ret))),
            "max_direct_vs_high_returned_Jp": float(np.max(np.abs(jp_direct - hp_jp_ret))),
            "warnings": dict(warning_counts),
        })
    return rows


def main():
    configs = {
        "coarse": dict(L=2e-6, d0=20e-9, transition=0.15e-6, h_max=40e-9),
        "medium": dict(L=2e-6, d0=5e-9, transition=0.30e-6, h_max=20e-9),
        "production": dict(L=2e-6),
    }
    rows = []
    for name, config in configs.items():
        rows.extend(_run_grid(name, config))
    payload = {
        "metadata": {
            "command": "MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/phase2c_s2_conservation_audit.py",
            "python": sys.version, "platform": platform.platform(),
            "solver_sha256": hashlib.sha256((ROOT / "semiconductor_sim/stage2_pn/solver.py").read_bytes()).hexdigest(),
            "official_golden_touched": False,
        },
        "formulae": {
            "direct": "Jn=qDn/h[B(d)n_j-B(-d)n_i], Jp=-qDp/h[B(-d)p_j-B(d)p_i]",
            "stable_qfp": "factor Bernoulli terms and evaluate quasi-Fermi differences with expm1",
            "high_precision": "Decimal(70 digits) SG evaluation from captured final QFP state",
            "local_electron": "Jn[i+1]-Jn[i]-q*h_i*U[i]",
            "local_hole": "Jp[i+1]-Jp[i]+q*h_i*U[i]",
        },
        "selected_bias_V": [-0.5, -1.0, -2.0],
        "records": rows,
    }
    out = ROOT / "verification/phase2c_s2_conservation.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fields = ["grid", "N", "V_a", "converged", "iters", "qfp_residual_l2", "qfp_residual_inf",
              "J_center_direct_A_m2", "direct", "stable", "high_precision_qfp",
              "max_direct_vs_stable_Jn", "max_direct_vs_stable_Jp", "max_stable_vs_high_Jn",
              "max_stable_vs_high_Jp", "max_direct_vs_high_returned_Jn", "max_direct_vs_high_returned_Jp"]
    with (ROOT / "verification/phase2c_s2_conservation.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, dict) else value
                             for key, value in row.items() if key in fields})
    print(json.dumps({"records": len(rows), "all_converged": all(r["converged"] for r in rows),
                      "warnings": sum(len(r["warnings"]) for r in rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
