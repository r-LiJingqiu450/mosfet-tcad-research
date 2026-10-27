# -*- coding: utf-8 -*-
"""
单截面纵向 Poisson（准平衡）求解器 —— 第三阶段 MOSFET 的核心引擎。

系列化说明：三对角求解统一走 common.thomas.thomas_solve_banded
（scipy LAPACK 带状优先、纯 Python 回退）；全部物理常数取自
common.constants（全系列唯一参数源）。

物理模型（电荷片 / 渐变沟道近似 GCA 的标准基础）：

    沿沟道 y 分 Ny 个截面；每个截面在"垂直于界面"的 x 方向做准平衡
    （quasi-equilibrium）求解：载流子随纵向能带弯曲按 Boltzmann 分布，
    电子准费米势 phi_n(y) 在该截面内取常数（由源漏/沟道决定），
    空穴与衬底平衡（phi_p = 0）：

        d^2 phi/dx^2 = -(q/eps_si)(p - n - Na)
        p = p0 exp(-phi/V_T)
        n = n0 exp(phi/V_T) exp(-phi_n/V_T)

    边界：phi(0)=phi_s（栅界面），phi(x_max)=0（衬底体内参考）。
    栅压由栅极混合边界反算：V_G = V_FB + phi_s - Q_s/C_ox。

约化模型的范围：
    源漏以理想电子库边界给出，不显式求解源漏结。GCA 是计算模型的
    选择；低电子浓度不意味着完整 DD 方程必然无解。phi_n=0 时退化为 S1。

约定（与第一阶段 mos_capacitor.py 一致）：
    phi 相对衬底体内（phi_bulk = 0），x 指向衬底。
"""

import warnings

import numpy as np

from common import constants as C
from common.thomas import thomas_solve_banded


# ================================================================================
# Poisson 残差与 Newton
# ================================================================================

def _carriers(phi, phi_n, p0, n0, density_clip_limit=None):
    """Optional historical clipping is for paired diagnostics only, off in production."""
    effective = (phi if density_clip_limit is None
                 else np.clip(phi, -density_clip_limit, density_clip_limit))
    with np.errstate(over='raise', invalid='raise'):
        try:
            return (p0*np.exp(-effective/C.VT),
                    n0*np.exp((effective-phi_n)/C.VT))
        except FloatingPointError as exc:
            raise ValueError('Carrier exponent outside floating-point range') from exc


def _residual(phi, phi_n, p0, n0, Na_arr, aa, bb0, cc, density_clip_limit=None):
    """
    泊松残差 F_i（内部节点）：
        F = lap + (q/eps)(p - n - Na)  (strong-form residual)
    采用第一阶段相同的非均匀二阶中心差分 + 电荷权重。

    p = p0 exp(-phi/V_T),  n = n0 exp(phi/V_T) exp(-phi_n/V_T)。
    """
    p, n = _carriers(phi, phi_n, p0, n0, density_clip_limit)
    lap = aa * phi[:-2] + bb0 * phi[1:-1] + cc * phi[2:]
    return lap + (C.Q_E / C.EPS_SI) * (p[1:-1] - n[1:-1] - Na_arr[1:-1])


def solve_poisson_slice(phi_s, phi_n_ch, x, Na=C.NA, phi_init=None,
                        tol=1e-9, max_iter=200, clamp=0.3, warn=True,
                        damped=True, max_backtrack=12, density_clip_limit=None):
    """Poisson Newton with joint un-damped correction and scaled residual tests.

    A small accepted damped step alone is never considered convergence.
    The line search cannot silently accept a failed Armijo trial. Diagnostics
    report a componentwise backward error, not a physical model error bound.
    Cross-grid initial guesses require explicit coordinate interpolation by caller.
    """
    x = np.asarray(x, float)
    if (x.ndim != 1 or len(x) < 3 or np.any(np.diff(x) <= 0)
            or not np.all(np.isfinite(x)) or Na <= 0 or tol <= 0
            or max_iter < 1 or clamp <= 0 or max_backtrack < 1
            or not np.isfinite(phi_s+phi_n_ch)):
        raise ValueError('Invalid Poisson inputs')
    hL, hR = np.diff(x)[:-1], np.diff(x)[1:]
    aa = 2/(hL*(hL+hR)); bb0 = -2/(hL*hR); cc = 2/(hR*(hL+hR))
    p0 = .5*(Na+np.sqrt(Na**2+4*C.NI**2)); n0 = C.NI**2/p0
    Na_arr = np.full_like(x, Na); phi_n = np.full_like(x, phi_n_ch)
    if phi_init is None:
        ld = np.sqrt(C.EPS_SI*C.VT/(C.Q_E*Na))
        if phi_s > 0:
            w = min(np.sqrt(2*C.EPS_SI*phi_s/(C.Q_E*Na)),x[-1])
            phi = np.where(x<w,phi_s*(1-x/w)**2,0.)
        else:
            phi = phi_s*np.exp(-x/ld)
    else:
        phi = np.asarray(phi_init,float).copy()
        if phi.shape != x.shape or not np.all(np.isfinite(phi)):
            raise ValueError('Initial potential must be on the supplied x grid')
    phi[0], phi[-1] = phi_s, 0.
    row_scale = np.abs(bb0)*C.VT + C.Q_E*Na/C.EPS_SI
    potential_min=float(np.min(phi)); potential_max=float(np.max(phi))
    clip_evaluations=0
    def diagnostics(v):
        nonlocal potential_min,potential_max,clip_evaluations
        potential_min=min(potential_min,float(np.min(v)))
        potential_max=max(potential_max,float(np.max(v)))
        if density_clip_limit is not None:
            clip_evaluations+=int(np.count_nonzero(np.abs(v)>density_clip_limit))
        p,n = _carriers(v,phi_n,p0,n0,density_clip_limit)
        f = _residual(v,phi_n,p0,n0,Na_arr,aa,bb0,cc,density_clip_limit)
        diag = bb0-C.Q_E/C.EPS_SI/C.VT*(p[1:-1]+n[1:-1])
        mag = (np.abs(aa*v[:-2])+np.abs(bb0*v[1:-1])+np.abs(cc*v[2:])
               +C.Q_E/C.EPS_SI*(p[1:-1]+n[1:-1]+Na))
        merit = float(np.linalg.norm(f/row_scale))
        floor = float(64*np.finfo(float).eps*np.linalg.norm(mag/row_scale))
        residual_v = float(np.max(np.abs(f)/np.abs(diag)))
        berr = float(np.max(np.abs(f)/np.maximum(mag,1.)))
        return f,diag,merit,floor,residual_v,berr
    converged=False; reason='max_iterations'; backtracks=0; correction=float('inf')
    for niter in range(1,max_iter+1):
        F,diag,merit,floor,residual_v,berr=diagnostics(phi)
        low=aa.copy(); low[0]=0.; up=cc.copy(); up[-1]=0.
        delta=thomas_solve_banded(low,diag,up,-F)
        correction=float(np.max(np.abs(delta)))
        if correction <= tol and residual_v <= tol:
            converged=True; reason='correction_and_residual'; break
        delta *= min(1.,clamp/max(correction,np.finfo(float).tiny))
        accepted=False
        for k in range(max_backtrack if damped else 1):
            alpha=2.**(-k)
            trial=phi.copy(); trial[1:-1]+=alpha*delta
            trial_diag=diagnostics(trial)
            if (not damped or trial_diag[2] <= (1-2e-4*alpha)*merit
                    or (merit <= floor and trial_diag[2] <= 2*max(floor,1e-30))):
                phi=trial; accepted=True; backtracks+=k; break
        if not accepted:
            reason='line_search_failed'; break
    F,diag,merit,floor,residual_v,berr=diagnostics(phi)
    if not converged:
        low=aa.copy(); low[0]=0.; up=cc.copy(); up[-1]=0.
        correction=float(np.max(np.abs(thomas_solve_banded(low,diag,up,-F))))
        if correction <= tol and residual_v <= tol:
            converged=True; reason='correction_and_residual'
    if not converged and warn:
        warnings.warn(f'Poisson failed: {reason}, phi_s={phi_s}, phi_n={phi_n_ch}',
                      RuntimeWarning,stacklevel=2)
    p,n=_carriers(phi,phi_n,p0,n0,density_clip_limit)
    return dict(phi=phi,n=n,p=p,iters=niter,converged=converged,
                status=reason,correction_inf=correction,residual_volts=residual_v,
                backward_error=berr,scaled_residual=merit,roundoff_floor=floor,
                backtracks=backtracks,potential_min=potential_min,
                potential_max=potential_max,clip_evaluations=clip_evaluations,
                density_clip_limit=density_clip_limit)


# ================================================================================
# 掺杂剖面（原 doping.py，并入本模块：最小闭环仅一处使用）
# ================================================================================

def doping_profile(x, Na):
    """
    返回 (Na_arr, Nd_arr)（m^-3），长度与 x 一致。

    最小闭环为均匀 p 衬底：Na_arr = Na 常数，Nd_arr = 0。
    """
    Na_arr = np.full_like(x, Na, dtype=float)
    Nd_arr = np.zeros_like(x, dtype=float)
    return Na_arr, Nd_arr


# ================================================================================
# 统一入口（兼容 gate_boundary 调用）与后处理：电荷 / 面密度
# ================================================================================

def solve_slice(phi_s, phi_n_ch, x, Na=C.NA, state=None,
                tol=1e-9, max_iter=200, clamp=0.3, verbose=False,
                warn=True):
    """
    单截面求解统一入口。最小闭环用准平衡 Poisson。

    state 可选：上一个 (phi, n, p) 解，作为初值加速延续收敛。
    """
    phi_init = state[0] if state is not None else None
    res = solve_poisson_slice(phi_s, phi_n_ch, x, Na=Na,
                              phi_init=phi_init, tol=tol,
                              max_iter=max_iter, clamp=clamp,
                              warn=warn)
    return res


# ================================================================================
# 后处理：电荷 / 面密度
# ================================================================================

def surface_charge(res, x):
    """
    半导体总面电荷 Q_s = eps_si dphi/dx|_0（三点一阶前向差分，同第一阶段）。
    返回 (Q_s, E_s)。
    """
    phi = res["phi"]
    h0 = x[1] - x[0]
    dphi_dx0 = (-3.0 * phi[0] + 4.0 * phi[1] - phi[2]) / (2.0 * h0)
    Q_s = C.EPS_SI * dphi_dx0
    E_s = -dphi_dx0
    return Q_s, E_s


def electron_sheet_density(res, x):
    """电子面密度 N_s = ∫ n dx（m^-2），梯形积分。"""
    n = res["n"]
    hw = np.diff(x)
    return float(np.sum(0.5 * (n[:-1] + n[1:]) * hw))


def hole_sheet_density(res, x):
    """空穴面密度 ∫ p dx（m^-2）。"""
    p = res["p"]
    hw = np.diff(x)
    return float(np.sum(0.5 * (p[:-1] + p[1:]) * hw))


def total_charge_density(res, x, Na=C.NA):
    """总空间电荷面密度 Q_s = q ∫ (p - n - Na) dx（C/m^2），供交叉校验。"""
    p = res["p"]
    n = res["n"]
    Na_arr = np.full_like(x, Na, dtype=float)
    rho = C.Q_E * (p - n - Na_arr)
    hw = np.diff(x)
    return float(np.sum(0.5 * (rho[:-1] + rho[1:]) * hw))
