# Phase 3 选题重置：有限 n⁺ 井深的最小决定性试验

**状态：`PARTIAL / MORE_EVIDENCE_NEEDED`**  
**分支：** `phase3/direction-reset-2026-11-01`  
**基线：** `fde679b561bddf35663933890d7ace617bf9b928`  
**本轮提交：** 见 Git 日志（本报告提交后更新）

## 科学决定

上一轮有限漏偏压、双长度和 Ward 对照已经完成为探索性检查，但没有形成可辩护的原创方法。本轮比较了三个科学上不同的方向：

1. SG 离散熵/自由能证书：相关结构和守恒分析已有成熟先例，适合作为可靠性工具，暂不构成器件物理贡献。
2. 耦合 Newton 条件数/阻尼切换：线搜索、continuation 和预条件已有广泛先例；当前没有可检验的新物理预测。
3. 有限 n⁺ 井深对 R-II 端修正的影响：Ward 1988 §6.4/§6.5 使用半无限井假设，而其 p.144 深度估计依赖缺失的 p.142。项目 C6 的旧 0.5→1.0 μm 比较还受深度网格吸附影响。

选择第 3 项，因为它改变的是物理边界规格，而不是继续增加同一条旧路线的网格、长度或偏压。可证伪的问题是：有限井深能否用耗尽宽度 `W_dep` 归一化，并在某个明确比例后替代半无限井？

最近的文献证据是仓库中的 Ward 审计：M. J. Ward, *Asymptotic Methods in Semiconductor Device Modeling*, Caltech thesis (1988), DOI `10.7907/08TK-NM84`；Ward–Odeh–Cohen, *SIAM J. Appl. Math.* 50 (1990), DOI `10.1137/0150066`。审计只把全文页码和公式作为 E1 证据；它没有把本试验的数值结果当成已有文献结论。

## 预注册规格与方法

固定 F0/R-II 规格：`N_A=10^16 cm^-3`、`t_ox=20 nm`、突变 n⁺、`L=2 μm`、level 0、`V_G=1.25 V`、恒定迁移率 `1350 cm²/(V s)`，`σ=0`。解析耗尽尺度为

\[
W_{dep}=0.3039293463\;\mu\mathrm m.
\]

井深断点为 `h_w={0.10, 0.25, 0.50, 1.00} μm`。`well_depth_probe.py` 在横向网格不变的条件下插入精确的深度断点，并通过 equilibrium/Newton、网络分解和可选 DEVSIM 路径记录 `ΔL_net`、`G_net`、网络残差和收敛信息。该脚本没有改动生产 S1/S2/S3 求解器或任何 golden 文件。

预注册判据为：在 `h_w/W_dep ≥ 3` 后，相邻深度的 `ΔL_net` 相对变化不超过 1%，且网络残差低于 `10^-10`。这只是本试验的验收门槛，不是普适误差界。

## 原始结果

| `h_w` (μm) | `h_w/W_dep` | `ΔL_net` (nm) | `G_net` (S/cm) | Newton 步数 | 网络残差 | KCL 相对残差 | 相邻 `ΔL` 变化 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 0.329024 | 98.6930761 | 1.9279023974e-1 | 14 | 3.78e-15 | 3.32e-8 | — |
| 0.25 | 0.822560 | 137.2289367 | 2.0140126003e-1 | 17 | 2.96e-13 | 1.33e-8 | 28.0814% |
| 0.50 | 1.645119 | 153.3487890 | 2.0523584709e-1 | 18 | 6.44e-14 | 1.93e-8 | 10.5119% |
| 1.00 | 3.290238 | 155.6384568 | 2.0579238756e-1 | 18 | 4.36e-14 | 3.32e-8 | **1.47115%** |

四个命令均返回退出码 0；每个原始 JSON 同时记录网格节点规模、网格间距、Newton 末步和壁钟时间。完整逐项数据在 `results/phase3/direction_reset/well_depth/well_depth_rows.csv`，原始求解在同目录四个 `well_depth_*.json` 文件。

**决定性结果：** `h_w/W_dep=3.29` 时，0.50→1.00 μm 的 `ΔL_net` 变化为 1.47115%，超过预注册的 1% 门槛；网络残差本身已远低于门槛。因此本试验拒绝了“在 3 个耗尽宽度处即可安全替代半无限井”的具体主张。结果仍显示单调趋于饱和，但尚未达到可用于边界替换的误差标准。

## 独立 DEVSIM 交叉检查

在 `h_w=1.00 μm` 的相同规格下，使用 `devsim==2.11.2` 的耦合路径独立重算。网络电导为 `0.2057923875644292 S/cm`，DEVSIM 两种提取分别为 `0.20579232307998246` 和 `0.20579245631365636 S/cm`，相对差为 `1.57e-7` 和 `3.34e-7`。DEVSIM 输出的条件数代理为 `3.46e10`，双重加权残差和舍入估计分别为 `3.57e-7` 与 `5.52e-8`。

这支持“在该离散规格下两条实现路径数值相容”的判断，但不等于独立物理验证；模型、边界和材料仍相同，且系统明显病态。DEVSIM 原始输出保存在 `logs/phase3/direction_reset/well_depth/hw1.00_devsim.log`，机器可读结果在 `results/phase3/direction_reset/well_depth_devsim/`。

## 代码、命令和测试证据

新增的隔离代码为 `phase3/direction_reset/well_depth_probe.py` 与 `phase3/direction_reset/analyze_well_depth.py`。代表命令如下：

```text
DEVSIM_MATH_LIBS='libblas.so.3:liblapack.so.3' .venv-devsim/bin/python -u phase3/direction_reset/well_depth_probe.py --L 2 --level 0 --vg 1.25 --hw 1.0 --na 1e16 --tox 20 --sigma 0 --out results/phase3/direction_reset/well_depth_devsim
MPLCONFIGDIR=/tmp/mplconfig python phase3/direction_reset/analyze_well_depth.py
python -m py_compile phase3/direction_reset/well_depth_probe.py phase3/direction_reset/analyze_well_depth.py
.venv-devsim/bin/python -m pytest -q
```

退出码分别为 `0`、`0`、`0`、`0`。当前仓库实际测试为 **94 passed, 4 skipped, 4 warnings**（27.93 s）；warning 是已有 S2 退化测试中的 Armijo/Newton/Gummel 诊断，没有被隐藏。完整日志为 `logs/phase3/direction_reset/full_pytest.log`。测试统计与较早 G-0 报告中的 100/171 项不同，是当前检出树实际收集到的测试数量，本报告不把历史统计冒充本次结果。

## 解释、限制和下一步门槛

这四点只覆盖一个 `N_A`、一个栅压、一个主网格 level 和一个恒定迁移率设置；没有证明跨掺杂、跨偏压或跨网格的相似律。`ΔL_net` 的趋于饱和可以是有限井物理、边界条件和离散误差的组合。DEVSIM 交叉检查也共享同一物理规格，不能排除模型假设错误。

因此当前没有足够证据宣称原创算法或 SCI 贡献。最强的可继续候选是“有限井深修正/适用性图”，但只有在以下独立结果出现后才值得投入：level 加密仍保留趋势；至少一个独立有限体积或不同离散路径复现趋势；多个 `N_A,V_G` 工况能给出可预测而非事后拟合的尺度关系；并且文献核查确认该判据不是已有结果的直接重述。若这些条件失败，应该把本轮结论保留为对半无限井假设的负面限制，而不是继续弱化阈值。

本轮没有实施生产代码修复，没有改变正式 golden、候选 baseline 或 `main`。旧有限偏压路线的结果保留在历史文件中，但不再作为默认研究主题。

## 证据索引

- 预注册选择和门槛：`phase3/direction_reset/Direction_Reset_Decision_2026-11-01.md`。
- 原始求解：`results/phase3/direction_reset/well_depth/*.json`。
- 汇总 CSV/JSON/图：`results/phase3/direction_reset/well_depth/well_depth_rows.csv`、`well_depth_summary.json`、`well_depth_saturation.png`。
- DEVSIM 日志和 JSON：`logs/phase3/direction_reset/well_depth/hw1.00_devsim.log`、`results/phase3/direction_reset/well_depth_devsim/`。
- 最近文献和假设审计：`phase3/reports/C5_Ward_FullText_Audit_2026-10-31.md`。
