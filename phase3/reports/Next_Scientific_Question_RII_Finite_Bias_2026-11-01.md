# R-II 有限漏偏压与长度仿射性：独立探索报告

## 决策

本轮没有提出新的求解器或新的提取算法。研究问题是：**C8/C9 已讨论的双长度端部项抵消，在解析 n⁺ 源漏区（R-II）和有限漏偏压下是否仍保持可检验的长度结构？** 这是对既有零偏结果的物理边界测试，而不是把双长度公式重新包装成创新。

结果为 `MORE_EVIDENCE_NEEDED`。在 `V_G=1.32, 1.60 V`、`L=1,2,4 µm`、`V_D=0.1–20 mV` 的小规模数据中，

\[
R(L,V_D)=\left[I(L,V_D)/V_D\right]^{-1}
\]

对 L 的三点拟合保持近似仿射；仿射斜率随 `V_D` 增大而变化。用低阶多项式把斜率外推到 `V_D=0`，得到的斜率分别与独立零偏网络值相差 `4.9×10⁻⁷`（1.60 V）和 `1.6×10⁻⁶`（1.32 V）。这支持一个**有限偏压下的条件性数值观察**，但不证明新的物理定律，也不足以宣称 SCI 原创性。

## 与最近工作的关系

C8/C9 已明确讨论双长度差分和其提取含义；因此本轮三长度试验的实质增量是：把问题移到解析 n⁺ 源漏区的 R-II 参考，并直接观察有限 `V_D` 的响应。最近的项目文献审查还指出 Ward 1988 §6.4 已研究线性区有限漏偏压的端部效应。当前结果应理解为对该物理问题的数值核查候选，而非与 Ward 理论竞争的新公式。

## 模型和可复现设置

- DEVSIM 2.11.2，Python 3.12.14，NumPy 2.3.5，SciPy 1.17.0；BLAS/LAPACK 通过 `DEVSIM_MATH_LIBS='libblas.so.3:liblapack.so.3'` 加载。
- R-II：解析 n⁺ 源漏区，`N_A=10¹⁸ cm⁻³`，侧向 erfc 梯度 `σ=15 nm`，`t_ox=2 nm`，`h_w=0.5 µm`，level 0，恒定 `μ=1350 cm²/(V s)`。
- 长度 `L={1,2,4} µm`；电势边界保持平衡 `PsiLat`（`s_couple=0`），漏端电子准费米边界继续到 `V_D`；源端为 0 V。
- 每个漏偏压从零开始分 5 步 Newton continuation；生产收敛判据是 `max(|Δψ|, V_T·relative electron update) ≤ 10⁻¹⁰ V`，并额外迭代一次。
- 每个状态检查漏/源端电流、`V1.I`、所有内部 x 截面的 SG 电流，而不是只检查中心截面。

运行命令（每个长度一条，实际退出码均为 0）：

```bash
DEVSIM_MATH_LIBS='libblas.so.3:liblapack.so.3' .venv-devsim/bin/python -u \
  phase3/next_question/physical_r2_finite_bias.py \
  --L {1,2,4} --level 0 --vg {1.32,1.60} \
  --vd 0.0001 0.001 0.005 0.020 --out results/phase3/next_question/physical_r2_finite_bias
python phase3/next_question/analyze_physical_r2.py \
  --out-dir results/phase3/next_question/physical_r2_finite_bias/analysis
```

实际记录分别保存在 `logs/phase3/next_question/r2_*.log`；原始 JSON 在 `results/phase3/next_question/physical_r2_finite_bias/`。

## 数值结果

### 长度仿射性

CSV：`results/phase3/next_question/physical_r2_finite_bias/analysis/physical_r2_finite_bias_rows.csv`。

| `V_G` | `V_D` | 仿射斜率 `b` | 最大三点相对残差 | `L=1–2` 与 `L=2–4` 斜率相对差 |
|---:|---:|---:|---:|---:|
| 1.60 V | 0.1 mV | 2088.0445 | 1.06×10⁻¹¹ | 3.8×10⁻¹¹ |
| 1.60 V | 1 mV | 2092.6693 | 1.09×10⁻⁹ | 3.8×10⁻⁹ |
| 1.60 V | 5 mV | 2113.4496 | 2.75×10⁻⁸ | 9.7×10⁻⁸ |
| 1.60 V | 20 mV | 2194.7970 | 4.65×10⁻⁷ | 1.6×10⁻⁶ |
| 1.32 V | 0.1 mV | 32141.4567 | 3.37×10⁻¹¹ | 1.1×10⁻¹⁰ |
| 1.32 V | 5 mV | 34320.7708 | 8.04×10⁻⁸ | 2.6×10⁻⁷ |
| 1.32 V | 20 mV | 41955.2025 | 1.06×10⁻⁶ | 3.4×10⁻⁶ |

单位为 `cm/(S/cm)`；这里的斜率只表示 `R` 对 `L` 的数值斜率，不直接等同于实验迁移率。

### 零偏极限

独立零偏 R-II 网络的斜率为 `2087.5317902`（`V_G=1.60 V`）和 `32098.5435682`（`V_G=1.32 V`）。对有限偏压斜率使用二次外推（仅作有限数据诊断）：

| `V_G` | 外推 `b(0)` | 相对零偏网络值 |
|---:|---:|---:|
| 1.60 V | 2087.5328062 | `+4.87×10⁻⁷` |
| 1.32 V | 32098.5941 | `+1.57×10⁻⁶` |

这说明有限偏压数据连续地回到零偏基线；它不能单独证明外推阶数，也不能替代更小偏压、网格加密和独立离散路径。

### 守恒和收敛

- `V_G=1.60 V`：全域截面绝对电流幅值离散最大相对展宽为 `6.18×10⁻⁷`（0.1 mV，L=4 µm），在 20 mV 降到 `3.63×10⁻⁹`；端点 KCL 最大相对误差为 `4.80×10⁻⁹`（0.1 mV，L=1 µm）。
- `V_G=1.32 V`：0.1 mV 的电流很小，截面展宽最高 `7.44×10⁻⁶`、KCL 最高 `3.78×10⁻⁷`；5 mV 时分别为 `2.55×10⁻⁷` 和 `5.32×10⁻⁹`，20 mV 时分别为 `8.33×10⁻⁸` 和 `2.45×10⁻⁹`。因此近阈值 0.1 mV 点只能作结构性观察，不能作高精度守恒证据。
- 所有已完成的 6 个长度/栅压运行退出码为 0，Newton 每一步均有收敛记录；日志保留 DEVSIM 的真实迭代输出。未隐藏 warning，也未把失败调用当成成功。

### R-II 电势边界规格敏感性

为避免把 `s_couple=0` 默认为唯一物理规格，另外在 `V_G=1.60 V`、相同三长度和 `V_D=5,20 mV` 下测试了 `s_couple=1`（漏端电势随外加电压平移）。“平衡电势”与“随偏压平移”两种规格的最大 `G_eff` 相对差分别为 `1.12×10⁻⁶` 和 `5.09×10⁻⁶`；三长度斜率差为 `3.3×10⁻⁷` 和 `6.7×10⁻⁶`（绝对斜率单位）。这些差异低于本试验的电流/守恒不确定度尺度，表现为测试窗口内的弱敏感性，但不能选择唯一物理边界，也不能证明两者数学等价。原始数据和汇总在 `results/phase3/next_question/physical_r2_boundary_sensitivity/`，求解日志在 `logs/phase3/next_question/boundary_sensitivity/`。

## 科学解释和限制

当前证据最支持的解释是：在这一窄小模型窗口，有限偏压改变了有效内区电阻斜率 `b(V_D)`，但端部层的主要贡献仍可由一个近似与 L 无关的截距吸收。`b(V_D)` 的偏压依赖本身是器件响应，不应被误报为端部项消除后的“真实迁移率”。

尚未建立的内容包括：独立有限体积/有限元重算、网格 level 1–3、非恒定迁移率、Fermi–Dirac 或量子修正、更多栅压与较大偏压、以及和 Ward §6.4 系数的定量对照。`s_couple` 敏感性测试只覆盖一个栅压和两个偏压，不能替代物理接触推导。特别是 `V_G=1.32 V, V_D=0.1 mV` 的守恒误差已经达到 `10⁻⁶–10⁻⁷`，不能支撑更强结论。

因此，本轮不修改生产 S1/S2/S3、不更新任何 golden、不启动 Package E，也不提出新的物理模型。若继续投资，最有判别力的下一步是固定一个近阈值点，做 level 0–2 的独立有限体积电流和截面守恒交叉验证，并把 `b(V_D)` 与 Ward 的有限偏压渐近项逐项对应；若对应不上，方向应保留为验证性结果而不是论文主线。

## 产物索引

- 求解器探针：[`phase3/next_question/physical_r2_finite_bias.py`](../next_question/physical_r2_finite_bias.py)
- 后处理：[`phase3/next_question/analyze_physical_r2.py`](../next_question/analyze_physical_r2.py)
- 边界规格分析：[`phase3/next_question/analyze_boundary_sensitivity.py`](../next_question/analyze_boundary_sensitivity.py)
- 原始 JSON：`results/phase3/next_question/physical_r2_finite_bias/`
- 误差表：`results/phase3/next_question/physical_r2_finite_bias/analysis/physical_r2_finite_bias_rows.csv`
- 摘要：`results/phase3/next_question/physical_r2_finite_bias/analysis/physical_r2_finite_bias_summary.json`
- 图：`results/phase3/next_question/physical_r2_finite_bias/analysis/physical_r2_finite_bias.png`
- 日志：`logs/phase3/next_question/r2_*.log`
