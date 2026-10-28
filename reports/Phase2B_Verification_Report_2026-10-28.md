# Phase 2B 数值验证与候选基线报告

日期：2026-10-28（Asia/Shanghai）  
分支：`phase2b/verification-candidates-2026-10-28`  
基线：`phase2a/s2-convergence-audit-2026-10-28`，`fb7a311303f3b366b556d2a5385def0b9c256399`  
最终判定：**PARTIAL**  
建议 Phase 2C：**有条件进入验证研究，不进入新的物理模型开发**。

本轮建立了 S1/S2 独立数值验证、候选 baseline 和文献筛查。官方 `stage1.json`、`stage2.json` 仍然缺失，没有写入 `semiconductor_sim/results/golden/`；官方 S3 golden 也保持不变。生产 S2 Newton 在指定偏压延续路径上满足其 QFP 残差和步长停止条件，但低反偏电流后处理的独立守恒检查存在网格和大数相消敏感性，因此不能把整个 S2 物理验证标为 PASS。

## A/B：基线保护与 Gummel 分类

* 当前分支从 `fb7a311` 创建，`main` 仍为 `c64c463ba35fb3c2f177c9558217a20387046795`。
* 工作区起始时干净；`main` 没有修改。
* `semiconductor_sim/results/golden/stage3.json` SHA-256 为 `17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828`。
* S1/S2 官方 golden 不存在，本轮没有生成或覆盖。
* Phase 2A 的 Gummel 结果保留为诊断例子：`omega=1` 和 `omega=0.5` 的 30 步外层分别停在 `1.465e-2 V` 和 `7.324e-3 V`，仍标记为 deprecated/nonconvergent。候选 S2 baseline 记录了这四类 warning 的来源，candidate test 也要求 warning 可见。

基线环境为 Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0、Matplotlib 3.10.8、pytest 9.1.1。完整运行命令和退出码见各 `verification/phase2b_*_command.log` 文件。

## C：S1 验证结果

独立脚本为 `scripts/phase2b_s1_audit.py`，输出见 `verification/phase2b_s1_verification.json` 和 `verification/phase2b_s1_metrics.csv`。它计算了离散残差、有限体积内部电荷积分、边界通量、表面三点差分电荷、解析耗尽近似和网格敏感性。

候选配置 `NGRID=201` 的标准数字为：

| 量 | 数值 | 解析/参考 | 相对误差 |
|---|---:|---:|---:|
| `Vt_numeric` | 0.9962624863 V | 0.9963499971 V | 8.78e-5 |
| `w_num` | 303.099 nm | 303.929 nm | 0.273% |
| `Cfb_numeric` | 0.593714 | 0.594739 | 0.172% |
| Newton 平均/最大步数 | 4.05 / 17 | 全点收敛 | 全部 `converged=True` |

网格加密到 `N=201,401,801,1601,3201` 时，`V_G(2ψ_F)` 为 `0.9953359, 0.9960196, 0.9962553, 0.9963246, 0.9963434 V`，相对最密网格的误差逐步下降；相邻观测阶约 `1.64, 1.88, 2.23`。这支持二阶离散趋势，但最密网格仍使用有限长度和边界条件，不能把它作为精确解析解。

残差与电荷检查需要分清定义：

* 内部有限体积电荷积分与离散边界通量的相对误差为约 `1e-14` 到 `1e-13`，这是离散方程的自洽检查。
* `postprocess.py` 的三点表面导数在 `φ_s=-0.5 V` 的相对电荷差约 `39.9%`，在 `φ_s=0.85 V` 约 `9.0%`；这不是 Poisson 残差失败，而是表面导数后处理和强边界层解析度的敏感性。S1 不能把所有 C-V/电荷点都称为独立物理验证。
* 归一化残差 `||F||_2 / [(q/ε)N_a√M]` 在审计点约 `1e-13` 到 `2e-11`。求解器的 `converged` 主要由步长决定，本轮报告同时保留残差，避免把步长停止单独解释为物理精度。

## D：S2 生产 Newton 验证

独立脚本为 `scripts/phase2b_s2_audit.py`，配置了 coarse、medium 和 production 三套结区网格；生产网格为 `make_junction_grid(L=2 µm)`。结果见：

* `verification/phase2b_s2_verification.json`
* `verification/phase2b_s2_metrics.csv`
* `verification/phase2b_s2_convergence_history.csv`
* `verification/phase2b_s2_command.log`

在每套网格上，从平衡态沿 `0.1 V` 步长推进到 `+0.5 V` 和 `-2.0 V`。78 个生产/网格/偏压记录全部返回 `converged=True`，无 warning；每条记录含 QFP 残差、终止步数、最后更新量、总电流、电子/空穴电流、SRH 积分、耗尽宽度和 Newton 历史。

生产网格的代表结果：

| 偏压 | Newton 步 | QFP `||F||₂` | `J` (A/m²) | `q∫U dx` (A/m²) | 电子差/`q∫U` |
|---:|---:|---:|---:|---:|---:|
| 0 V | 8 | 3.84e-10 | 4.04e-9 | 4.04e-9 | 由于平衡态抵消，比例不稳定 |
| +0.5 V | 8 | 3.26e-11 | 1.2142e4 | 1.0879e4 | 0.9999996 |
| -0.5 V | 6 | 4.85e-10 | -1.1749 | -1.17488 | 1.0000004 |
| -1.0 V | 6 | 2.26e-10 | -2.0323 | -2.03225 | 约 2.12，受端点 SG 相消影响 |
| -2.0 V | 5 | 4.99e-11 | -3.4199 | -3.41988 | 1.0205 |

QFP 残差和 `newton_solve` 的 `converged` 是数值求解器判据；电子/空穴守恒比值是独立后处理检查。`-1 V` 附近的端点差分比值不满足 2% 检查，而总电流和 SRH 积分仍接近。这一差异需要在后续工作中改进 SG 电流后处理或采用控制体积端点通量，当前不能写成完整物理验证通过。

网格对比显示，`+0.5 V` 的电流从 coarse 到 production 约为 `11978 -> 12142 A/m²`，`-1 V` 约为 `-2.028 -> -2.032 A/m²`；Newton 都收敛，但低电流区守恒误差更敏感。脚本同时保存原始残差分量，供后续决定是否需要独立端点通量实现。

## E：候选 baseline

候选文件位于 `verification/candidates/phase2b_v1_0/`：

* `stage1.json` SHA-256：`93244365320f005eb61ae242f5416cbb4cf9a9af92b1e63917370e336adc0bd0`
* `stage2.json` SHA-256：`7883dfdd80fe5fc3b768ed3e21088f6b212682dd74048c033a935dd1162afeac`
* `manifest.json`：版本、配置、代码基线、容差、warning、独立证据路径和官方 golden 状态。

这些文件明确标记 `candidate-not-official-golden`，来源是 `semiconductor_sim/tests/golden_helpers.py` 的降配配置。S1 配置是 `NGRID=201`；S2 配置是 `L=2 µm, V=+0.5/-1.0 V`。它们是当前代码的可复现候选回归线，不是历史实验结果。`semiconductor_sim/tests/test_phase2b_candidate_regression.py` 验证候选文件哈希、配置隔离和数值重现，并要求 S2 的 Gummel warning 保持可见。

候选测试命令：

```text
MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest semiconductor_sim/tests/test_phase2b_candidate_regression.py -q
exit 0; 3 passed
```

加入候选回归测试后的全量命令为：

```text
MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest semiconductor_sim -q
exit 0; 94 passed, 4 skipped, 4 warnings
```

## F：文献与 SCI 问题筛查

文献报告为 [Phase2B_Literature_Screening_2026-10-28.md](/workspace/mosfet-tcad-research/reports/Phase2B_Literature_Screening_2026-10-28.md)，机器可读来源为 `verification/literature_screening_sources.json`。记录了 Gummel (1964)、Scharfetter–Gummel (1969)、Pao–Sah (1966)、Tsividis/McAndrew、Roy (2005)、Oberkampf–Trucano、Celik、Selberherr、Markowich 等来源及 DOI/出版者链接。Crossref API 在当前代理中返回 403，因此报告保留了可核验链接，没有把搜索摘要当成已验证全文。

可验证的研究问题包括：GCA 误差分解、阈值分支约定敏感性、CLM 与离散误差的符号分解、1D 到独立 2D/商业 TCAD 的转移性、以及网格/残差/参考解联合不确定度。当前仓库没有测量数据、独立 2D 求解器或已冻结商业 TCAD deck，因此不作原创性或 SCI 分区承诺。

## 测试、风险和 Phase 2C 决策

本轮新增脚本均通过 `py_compile`，S1/S2 审计脚本退出码为 0。加入候选测试后的全量测试为 94 passed、4 skipped、4 warnings；Phase 2A 原始基线仍为 91 passed、4 skipped、4 warnings。S1/S2 官方 golden 继续跳过。生产 S2 Newton 的偏压延续在本轮三套网格上返回收敛，但低反偏端点电流后处理没有达到统一的 2% 独立守恒标准。

因此本轮为 **PARTIAL**：

* 可以进入 Phase 2C 的验证研究：先实现端点控制体积电流/守恒后处理，对 `-1 V` 低电流区重复检查；
* 可以继续使用候选 baseline 做回归，但不得把它们移动到官方 golden 目录；
* 暂缓新物理模型、2D/短沟道结论和论文级精度声明，直到独立参考与电流后处理问题解决。
