# Codex 独立探索：守恒网络与小信号条件性闸门（2026-11-01）

## 范围与状态

本探索只读取 C10/C11 已冻结的 JSON 结果，没有重新运行二维 MOSFET 扫描，也没有修改生产求解器、正式 golden 或历史结果。当前结论是 `MORE_EVIDENCE_NEEDED`，不是新模型通过，也不是 SCI 新颖性结论。

- 代码提交：`e61a443f2ef4c492b63e8b2c9fa827bd3145a6fb`
- 输入数据：`results/phase3/C10/analysis/c10_result.json`、`results/phase3/C11/c11_eval.json`、C11 DEVSIM 检查及独立有限体积结果。
- 机器可读输出：`results/phase3/exploration/conditioned_small_signal_rows.csv`、`results/phase3/exploration/coupled_small_signal_rows.csv`、`results/phase3/exploration/conditioned_small_signal_summary.json`。

## 可证伪问题

候选问题是：在深亚阈值、输出电流接近浮点分辨率时，能否用“平衡态电阻网络 + 独立守恒/残差检查 + 条件数或后验误差估计”判断耦合小信号结果是否可报告？可证伪标准是：在同一收敛态上，网络结果与独立有限体积参考不能保持小于预先给定误差，或条件性闸门无法识别已知的耦合失效状态，则该方向不成立。

## 已有证据的重新量化

| 指标 | 结果 | 证据等级 |
|---|---:|---|
| C11 H-lin 状态数 | 24 | 数值复核 |
| H-lin 绝对相对误差 | 0.133353–0.241581（中位数 0.196701） | 数值复核 |
| H-lin 误差符号 | 24/24 低估 | 数值复核 |
| abs(error) 与 log10(alpha) 相关系数 | 0.834 | 探索性统计 |
| abs(error) 与 log10(B∞/B_lin) 相关系数 | 0.869 | 仅诊断，B∞ 拟合不是统一参考 |
| DEVSIM 网络与独立有限体积 dL 相对差（2 点最大） | 0.0014866 | 独立数值路径 |
| 两个耦合小信号状态的最大 |G_M−G_E|/max(|G|) | 1.98779 | DEVSIM 后验诊断 |
| 两个耦合状态 κ_G | 2.402e+15–3.821e+16 | 条件性诊断 |

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
DEVSIM_MATH_LIBS=/tmp/devsim-libs/libblas.so:/tmp/devsim-libs/liblapack.so \
  LD_LIBRARY_PATH=/tmp/devsim-libs \
  /workspace/mosfet-tcad-research/.venv-devsim/bin/python -m pytest phase3/tests -q -rA
python3 phase3/exploration/conditioned_small_signal_audit.py
```

两次测试的实际日志位于：`/workspace/c12_restore_backup/evidence/c12_phase3_tests.log` 与 `/workspace/c12_restore_backup/evidence/c12_phase3_tests_devsim.log`。
