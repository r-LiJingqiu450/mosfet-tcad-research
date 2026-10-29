# Phase 2C：半导体数值守恒与后处理可信性审计

日期：2026-10-29（Asia/Shanghai）  
分支：`phase2c/conservation-postprocess-audit-2026-10-29`  
起点：`phase2b/verification-candidates-2026-10-28`，`e1598b9597e99fecefcf3a057ca86a9258f969b4`  
最终状态：**PARTIAL**

本轮没有修改生产物理模型、Newton 停止标准、官方 golden 或 `main`。审计固定同一批已收敛 S2 解，分别评估现有双精度 SG 通量、使用最终 QFP 变量的稳定 `expm1` 形式、70 位 Decimal 高精度参考和控制体积残差。结果显示：现有后处理在低反偏区会出现严重浮点相消；稳定公式与高精度 QFP 参考几乎一致，但局部守恒误差仍由求解器停止误差和低载流子行缩放主导。因此只修改后处理公式不能证明整个生产解达到物理守恒，本轮保留生产实现并交付误差证据。

## 基线保护

* 当前分支从指定 Phase 2B HEAD 创建，工作区初始干净。
* `main` 为 `c64c463ba35fb3c2f177c9558217a20387046795`，未修改。
* 官方 S3 golden SHA-256：`17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828`。
* 官方 S1/S2 golden 仍不存在；候选文件仍在 `verification/candidates/phase2b_v1_0/`，未移动到官方目录。

## A：S2 离散关系与误差定义

`_newton_system` 和 `current_density` 使用的现有 SG 形式为

\[
J_n=\frac{qD_n}{h}[B(d)n_{i+1}-B(-d)n_i],\qquad
J_p=-\frac{qD_p}{h}[B(-d)p_{i+1}-B(d)p_i],
\]

其中 `d=(φ_i-φ_{i+1})/V_T`，`B(d)=d/(1-e^{-d})`。内部节点控制体积残差是

\[
F_n=J_{n,i+1/2}-J_{n,i-1/2}-q h_i U_i,
\]
\[
F_p=J_{p,i+1/2}-J_{p,i-1/2}+q h_i U_i.
\]

因此全域积分关系应为

\[
J_n(x_R)-J_n(x_L)=q\sum_i h_iU_i,
\quad
J_p(x_R)-J_p(x_L)=-q\sum_i h_iU_i.
\]

审计同时记录：

* 局部电子/空穴残差的无穷范数及相对于 `max(|q∫U dx|, |J|)` 的比例；
* 两端电子/空穴电流差与 `q∫U dx` 的比值；
* 所有截面 `J_t=J_n+J_p` 的最小值、最大值、中位数和全域 spread；
* 现有双精度表达式与稳定表达式的最大差异；
* 稳定表达式与 Decimal 70 位 QFP 参考的差异。

稳定表达式把 Bernoulli 公因子和准费米势差分开，用 `expm1` 计算。例如电子通量写成

\[
J_n=\frac{qD_n}{h}B(d)n_i
e^{(\phi_{i+1}-\phi_{n,i})/V_T}
\operatorname{expm1}[(\phi_{n,i}-\phi_{n,i+1})/V_T],
\]

空穴通量使用对应的 `B(-d)` 和准费米势差。高精度参考从 Newton 最终 QFP 数组计算载流子，再用 Decimal 直接评估 SG 表达式；它不是新的物理模型。

## S2 结果

数据文件：

* `verification/phase2c_s2_conservation.json`
* `verification/phase2c_s2_conservation.csv`
* `verification/plots/phase2c_s2_conservation.png`
* `verification/plots/phase2c_s2_formula_error.png`

每个网格固定从平衡态以 `0.1 V` 步长推进到 `-0.5, -1.0, -2.0 V`，9 个目标状态全部满足生产 Newton 的 `converged=True`。下表中的第一列比值是现有双精度通量，第二列是最终 QFP 稳定形式，第三列是 Decimal QFP 参考；`spread` 是现有双精度总电流全域 spread 的相对值；`ΔJn` 是现有通量与稳定形式的最大绝对差。

| 网格 | Va | 现有 ΔJn/q∫U | 稳定 ΔJn/q∫U | 高精度 QFP | 局部电子误差 | 局部空穴误差 | 全域 spread | max ΔJn |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| coarse | -2.0 | 0.999981 | 0.995921 | 0.995921 | 1.29e-3 | 7.67e-4 | 2.34e-3 | 1.39e-2 |
| coarse | -1.0 | 1.002240 | 1.022335 | 1.022335 | 1.11e-2 | 4.92e-3 | 1.11e-2 | 4.08e-2 |
| coarse | -0.5 | 1.073941 | 0.998758 | 0.998758 | 1.18e-1 | 2.98e-2 | 1.18e-1 | 2.28e-1 |
| medium | -2.0 | 0.992654 | 0.993595 | 0.993595 | 2.28e-2 | 7.86e-3 | 2.28e-2 | 2.72e-2 |
| medium | -1.0 | 1.199972 | 1.248139 | 1.248139 | 1.74e-1 | 6.29e-2 | 2.66e-1 | 2.02e-1 |
| medium | -0.5 | 0.820118 | 0.498839 | 0.498839 | 1.16 | 0.347 | 1.44 | 1.58 |
| production | -2.0 | 1.054499 | 1.022558 | 1.022558 | 0.239 | 6.85e-2 | 0.286 | 0.512 |
| production | -1.0 | 1.235444 | 2.176000 | 2.176000 | 2.07 | 0.466 | 2.74 | 2.88 |
| production | -0.5 | 1.000000 | 1.000000 | 1.000000 | 6.96 | 2.19 | 6.96 | 8.24 |

这组数据区分了四种误差：

1. **浮点表达式误差。** 稳定形式和 70 位 QFP 参考的电子通量最大差异约 `1e-9 A/m²`，空穴通量约 `1e-6–1e-5 A/m²`；现有通量与高精度“返回载流子”参考的差异可达到数 `A/m²`，证明低电流区确有大数相消和载流子重构误差。
2. **后处理误差。** 现有 `current_density` 只接收 `(φ,n,p)`，丢失 Newton 内部的 QFP 数组；重新从双精度 `n,p` 反推 QFP 会把舍入噪声放大。稳定公式改善了公式评估，但并未让所有局部残差归零。
3. **网格误差。** coarse/medium/production 的比值和 spread 不呈简单单调关系，说明低电流区域不能用一次网格比较替代收敛研究。
4. **求解器停止误差。** QFP 缩放残差可很小，但低载流子行的反缩放会产生 `O(1 A/m²)` 的原始连续性残差；因此“Newton converged”与“器件电流守恒达到 2%”是两个不同条件。

根因明确但不适合在本轮做生产修复：将 `current_density` 改为接受 QFP 变量可以改善表达式评估，却不能消除低载流子行的停止误差；提高预算或放宽容差会改变数值验收标准。生产实现保留，问题转为 Phase 2D 候选：增加 QFP 状态接口并以独立局部守恒阈值回归保护。

## B：S1 表面电荷和阈值口径

数据文件：

* `verification/phase2c_s1_postprocess.json`
* `verification/phase2c_s1_postprocess.csv`
* `verification/plots/phase2c_s1_charge_methods.png`
* `verification/plots/phase2c_s1_threshold_methods.png`

三种定义为

\[
Q_{3pt}=\epsilon\frac{-3\phi_0+4\phi_1-\phi_2}{2h_0},
\]
\[
Q_{flux}=\epsilon\frac{\phi_1-\phi_0}{h_0},
\quad
Q_{int}=q\sum_{i\in interior}(p_i-n_i-N_A)h_i.
\]

离散控制体积关系是 `Q_flux - Q_right_flux = Q_int`，而非直接把接触边界节点的浓度用普通梯形积分并入空间电荷。审计中 `Q_flux` 与 `Q_int` 的相对差约 `1e-14–1e-13`；三点导数在强边界层下会偏离积分定义。

### 代表点的后处理偏差

| N | φs | 三点与积分相对差 | 通量与积分相对差 |
|---:|---:|---:|---:|
| 201 | -0.5 V | 41.4% | 3.8e-16 |
| 201 | 0.3 V | 1.35% | 2.4e-15 |
| 201 | 0.85 V | 15.4% | 2.2e-14 |
| 801 | -0.5 V | 39.9% | 4.7e-16 |
| 801 | 0.3 V | 0.333% | 5.1e-14 |
| 801 | 0.85 V | 9.00% | 4.4e-14 |
| 3201 | -0.5 V | 37.9% | 6.2e-16 |
| 3201 | 0.3 V | 0.0831% | 7.0e-13 |
| 3201 | 0.85 V | 3.07% | 2.7e-14 |

`φ_s=-0.5 V` 的通量和积分电荷随 N 仍明显变化，说明该强积累点在当前有限域和网格范围内尚未达到可接受的网格收敛；不能把它用作已收敛的物理参考。

### N=201 阈值口径

阈值扫描均取 `φ_s` 样本中最接近 `2ψ_F` 的 `0.715 V`，然后使用不同 `Q_s` 定义映射到 `V_G`：

| 方法 | V_T (N=201) | W (N=201) |
|---|---:|---:|
| 现有三点后处理 | 0.9962624863 V | 303.099 nm |
| 边界通量 | 0.9929703052 V | 299.551 nm |
| 体积积分 | 0.9929703034 V | 299.551 nm |

三点与通量阈值差约 `3.29 mV`；到 `N=801` 降至约 `1.05 mV`，`N=3201` 仍约 `0.286 mV`。这三条序列属于不同后处理口径，不能混合成一条“网格收敛曲线”。平带电容也受定义影响：N=201 三点为 `0.593714 C/Cox`，通量/积分约 `0.580110 C/Cox`。

本轮没有修改 `postprocess.py`。三点导数是生产历史口径，通量/积分是控制体积自洽口径；后续若更换生产口径，必须建立新版本候选 baseline 并记录差异。

## C：测试与完整性

实际命令：

```text
MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/phase2c_s2_conservation_audit.py
exit 0; 9 records, all converged, 0 warnings

MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/phase2c_s1_postprocess_audit.py
exit 0; 5 grids and requested surface-potential points written

MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/phase2c_make_plots.py
exit 0; four audit plots available across Phase 2C S1/S2 plot files

MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest semiconductor_sim -q
exit 0; 94 passed, 4 skipped, 4 warnings
```

4 条 warning 仍是 Phase 2A 已知的 deprecated Gummel 对照 warning。没有关闭、过滤或删除它们。审计脚本和测试均通过 `py_compile`，所有差异和命令输出分别保存在 `verification/phase2c_*`。

正式/候选基线完整性：

* 官方 `stage3.json` SHA-256 仍为 `17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828`；
* 官方 S1/S2 golden 仍缺失；
* Phase 2B 候选 `stage1.json` 为 `93244365320f005eb61ae242f5416cbb4cf9a9af92b1e63917370e336adc0bd0`；
* Phase 2B 候选 `stage2.json` 为 `7883dfdd80fe5fc3b768ed3e21088f6b212682dd74048c033a935dd1162afeac`。

## 科学结论与后续进入条件

代码验证已经证明离散 SG 残差、QFP Newton 终止和控制体积关系的实现路径可追踪；数值收敛验证显示三套网格上的 Newton 返回收敛；器件物理验证仍受低反偏电流后处理、QFP 状态丢失、低载流子行缩放和 S1 强积累网格未收敛限制。当前不能声称“全域电流守恒已通过”、不能把三点表面电荷当成唯一物理真值，也不能据此宣称 SCI 创新或实验精度。

进入后续 S3 高保真对标研究前，至少需要：

1. 为生产 S2 暴露或保存最终 QFP 状态，使用稳定通量和高精度对照重新定义端点/局部守恒验收；
2. 明确低载流子行的原始残差阈值，证明 solver stop error 小于物理电流误差预算；
3. 选择并冻结 S1 表面电荷口径，完成 `φ_s=-0.5 V` 的域长和网格收敛；
4. 为新口径建立独立版本候选 baseline，保持官方 golden 不变。
