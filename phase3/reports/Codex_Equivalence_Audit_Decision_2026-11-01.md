# 独立等价性审计与科研方向决策（2026-11-01）

## 决策

本轮将“守恒网络 + 条件性小信号闸门”从潜在论文方向降级为**数值可靠性审计协议**。当前证据支持低电流条件下的实用检查，但没有显示新的数学算法、独立误差定理或已建立的 SCI 级原创贡献。状态：`MORE_EVIDENCE_NEEDED`；小信号可靠性不再作为当前首选发表主线。

此前探索报告中把 `endprobe` 的原始 `G_net` 与 C11 的源端 `dL_2D` 放在同一表中，命名上造成了混淆。原始 `G_net` 是排端节点电流的直接求和，而 C11 的 `dL_2D` 使用了源端的 cancellation-free 电流。该问题没有被删除或改写；本报告给出纠正后的独立审计。

## 比较量的数学定义

在固定的零漏压平衡态、恒定迁移率和 SG 离散下，边电导为

\[
 c_e = q\,\mu\,g_e\,\bar n_e,
\]

其中 `g_e` 是几何耦合因子，`\bar n_e` 是平衡载流子权重。对源端电势 `\varphi=0`、漏端电势 `\varphi=1` 的离散网络，令 `L` 为对称加权图 Laplacian，内部节点满足

\[
 (L\varphi)_i=0.
\]

三种在精确算术中等价的网络量是：

1. 源端端点电流 `G_src = -Σ_{i∈source}(Lφ)_i`；
2. 漏端端点电流 `G_drain = Σ_{i∈drain}(Lφ)_i`；
3. 正定边能量 `G_energy = Σ_e c_e(φ_i−φ_j)^2 = φᵀLφ`。

有限精度下，`G_drain` 可能由大电流项相消而失真；边能量没有端点大数相减，源端和边能量是本审计的稳定网络量。

`G_E` 是行列平衡后的 SuperLU 解中 `-x[V1.I]`，`G_M` 是同一导出 Newton 线性系统的 UMFPACK 解中 `-x[V1.I]`。它们不是两个不同物理模型，而是同一全耦合 Jacobian 输出量的两种线性代数路径。只有在矩阵装配、变量缩放、边界条件和输出符号都对应时，才应与零偏 SG 网络电导比较。

独立有限体积结果采用独立网格和 Newton 实现，并通过同一类零偏网络计算 `dL`，因此它是数值路径独立的参考；它不是实验数据，也不是商业 TCAD 验证。

## 纠正后的新计算

配置：`W` 端模型，`L=1 µm`，`N_A=10^18 cm⁻³`，`t_ox=2 nm`，`h_w=0.5 µm`，level 0，恒定迁移率，`V_D→0`。每个状态使用新运行的 `endprobe.py --devsim`，并由包装脚本重新计算边能量。

| `V_G` (V) | `G_energy/G_src−1` | `G_drain/G_src−1` | `G_E/G_src−1` | `G_M/G_src−1` | `ΔL_energy` 与独立 FV |
|---:|---:|---:|---:|---:|---:|
| 0.488174 | `+2.60e−11` | `+1.1913e+02` | `−3.5823e+02` | `+3.6065e+02` | `+0.0430%` |
| 0.621433 | `−5.66e−11` | `+1.4516e+00` | `−7.5387e+00` | `+1.5506e+01` | `−0.1566%` |
| 0.783500 | `+1.78e−11` | `+1.3512e−03` | `−4.5868e−02` | `+6.6864e−02` | `−0.1487%` |

全部三点的边能量与源端电流相对差小于 `5.7e−11`。独立有限体积的 `ΔL` 差异小于 `0.157%`。相反，排端原始求和在深亚阈值达到约 `119×` 相对误差。`G_E` 与 `G_M` 在低电流点也严重偏离稳定网络量；它们的 equilibrated backward residual 仍约 `10^−16`，说明这是典型的前向条件性问题，而不是简单的残差不收敛。

本轮没有把三点拟合成普适阈值。它们只证明：

- 源端/边能量网络量在这三个固定配置下可重复；
- 排端电流后处理存在已知的浮点相消风险；
- 全耦合小信号输出的误差随状态条件性显著增大，但不能据此定义通用切换阈值。

机器可读结果：

- `results/phase3/exploration/equivalence_audit/equivalence_rows.csv`
- `results/phase3/exploration/equivalence_audit/equivalence_summary.json`
- `results/phase3/exploration/equivalence_energy_final/`
- `results/phase3/exploration/equivalence_na1e18/fd2d_VG0.621433.json`

## 误差来源判别

**浮点后处理误差。** `G_drain` 与 `G_src` 使用同一网络解，却相差最高约 `119×`；边能量仍与 `G_src` 一致。这直接定位了端点电流求和的相消，而不是物理响应差异。

**线性求解前向误差。** `G_E` 和 `G_M` 都来自同一导出 Jacobian，且 equilibrated residual 很小；它们对源端网络量的偏差达到数百倍。`κ_G` 在原 C11 两点约 `2.4×10^15–3.8×10^16`，与低电流输出对舍入扰动的敏感性一致。条件数和 DWR 是诊断估计，不是严格真误差界。

**模型/边界条件差异。** 网络与独立 FV 都是零偏、恒定迁移率、冻结平衡态、相同 W 端边界。该窗口内它们的 `ΔL` 一致，因此当前三点没有证据支持“网络/FV 差异来自不同物理量”。`G_E/G_M` 仍然是全 Jacobian 输出；若研究有限漏压、场依赖迁移率、复 AC 导纳或不同端接，它们不能自动视为同一个量。

**求解器停止误差。** 网络内部残差和边能量/源端一致性分别记录；当前差异远大于 `10^−11` 的网络恒等式误差，说明主要问题不在网络内部线性方程的停止判据，而在输出量的条件性与端点后处理。

## 文献与新颖性判断

本仓库已有文献筛查提供了足够的定位：

- Gummel (1964), DOI `10.1109/T-ED.1964.1536099`：Gummel 迭代的原始来源；说明固定点/全 Newton 对比属于成熟数值背景。
- Scharfetter–Gummel (1969), DOI `10.1109/T-ED.1969.16754`：指数拟合通量及其离散背景；零偏网络线性化是该离散的标准极限分析。
- Selberherr (1984), DOI `10.1007/978-3-7091-8752-4`；Markowich、Ringhofer、Schmeiser (1990), DOI `10.1007/978-3-7091-0650-2`：半导体漂移扩散、Poisson 边界和数值解的经典参考。
- Roy (2005), DOI `10.1016/j.jcp.2004.10.036`；Oberkampf–Trucano (2002), DOI `10.1016/S0376-0421(02)00005-2`；Celik et al. (2008), DOI `10.1115/1.2960953`：代码验证、解验证、网格不确定度和报告规范。
- C1b 现有报告引用了 Higham–Mary (2019, *SIAM J. Sci. Comput.* 41(5)) 的概率舍入误差框架，并明确将 Hager–Higham/Skeel 条件估计标为估计量而非误差上界。

这些来源共同说明：SG 线性化、能量/端点守恒、后向误差、条件数和 DWR 都属于成熟工具。当前贡献是把它们组合到本项目的低电流 MOSFET 输出审计中；这可以形成内部可靠性规范或验证性短文，但没有足够证据宣称新的数值算法。

本环境对 Crossref API 的请求被代理拒绝（HTTP 403），因此没有把未经独立核实的搜索片段加入参考文献。上列 DOI 来自仓库已有筛查报告和代码引用。

## 科研方向决定

我决定**停止把条件性小信号闸门作为当前发表主线**，保留以下工程规则：

1. 固定零偏网络量时优先报告源端电流和正定边能量；排端端点和必须通过相消敏感性检查。
2. 全 Jacobian 小信号输出必须同时报告平衡网络对照、后向残差、条件性估计和独立路径差异。
3. 不从两个或三个状态外推通用阈值；任何切换规则必须在预注册的新状态上验证。

下一条候选研究线应回到具有明确物理问题的 GCA/有限漏压或独立二维参考比较，并先完成文献闸门。当前不启动新的大规模二维扫描，也不宣称 SCI 原创性。

## 复现与验证

```bash
cd /workspace/c12_restore_backup/restored
DEVSIM_MATH_LIBS=/tmp/devsim-libs/libblas.so:/tmp/devsim-libs/liblapack.so \
  LD_LIBRARY_PATH=/tmp/devsim-libs \
  /workspace/mosfet-tcad-research/.venv-devsim/bin/python \
  phase3/exploration/equivalence_energy.py \
  --deck W --L 1 --level 0 --vg 0.488174 \
  --tox 2 --na 1e18 --hw 0.5 --devsim \
  --out results/phase3/exploration/equivalence_energy_final

/workspace/mosfet-tcad-research/.venv/bin/python \
  phase3/c11/fd2d_indep.py --na 1e18 --tox 2 --vg 0.621433 \
  --hw 0.5 --L 1 \
  --out results/phase3/exploration/equivalence_na1e18/fd2d_VG0.621433.json

/workspace/mosfet-tcad-research/.venv/bin/python \
  phase3/exploration/equivalence_audit.py
```

三次 `endprobe` 运行退出码均为 0；独立有限体积运行退出码为 0。旧的 preliminary 输出和新纠正输出均保留。
