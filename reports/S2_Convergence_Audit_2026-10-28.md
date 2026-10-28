# S2 收敛问题审计（Phase 2A）

日期：2026-10-28（Asia/Shanghai）  
分支：`phase2a/s2-convergence-audit-2026-10-28`  
基线：`restore/phase1-v1.2`，HEAD `78b38b714aedfb38bbb4e43d762e4744ec3cae11`  
本轮结论：**PARTIAL**

本轮没有修改 S2 生产求解器、物理参数、收敛阈值或 golden 文件。四条 warning 均能在当前代码中复现，并已定位到调用链和数值触发条件。生产路径使用的全耦合 QFP Newton 在平衡态收敛；warning 来自 `run.py` 中保留的 Gummel 方法对照路径，因此不能把测试通过解释为 Gummel 已收敛。由于 Gummel 固定点仍然失稳，本轮不把 warning 隐藏或改写为成功。

## 复现环境与基线保护

代码和环境检查结果：

* 分支为 `phase2a/s2-convergence-audit-2026-10-28`，从恢复基线创建；`main` 仍为 `c64c463ba35fb3c2f177c9558217a20387046795`。
* `semiconductor_sim/stage2_pn/solver.py` 和 `run.py` 在审计前后没有工作区差异。
* Python 3.12.14，NumPy 2.3.5，SciPy 1.17.0，Matplotlib 3.10.8，pytest 9.1.1。
* `semiconductor_sim/results/golden/stage3.json` 审计前后 SHA-256 均为 `17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828`。
* `stage1.json` 和 `stage2.json` 仍不存在；本轮没有生成或冻结它们。

定向测试命令和真实退出码：

```text
MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest semiconductor_sim/tests/test_stage_closure.py::test_s2_equilibrium_degeneration -vv -s
exit 0; 1 passed, 4 warnings

MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest semiconductor_sim/stage2_pn/tests -q
exit 0; 9 passed

MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest semiconductor_sim/tests/test_common.py semiconductor_sim/tests/test_stage_closure.py::test_s2_equilibrium_degeneration -q
exit 0; 27 passed, 4 warnings

MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python -m pytest semiconductor_sim -q
exit 0; 91 passed, 4 skipped, 4 warnings

MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/reproduce_phase1.py --fresh
exit 0; collected 95, passed 91, skipped 4, failures 0, errors 0
```

完整输出分别保存在 `verification/s2_convergence_test.log`、`verification/s2_targeted_tests.log`、`verification/s2_kernel_closure_tests.log`、`verification/s2_full_pytest.log` 和 `verification/s2_reproduce_phase1.log`。fresh 流程产生的时间和图形元数据曾改变工作区中的派生结果；为保护恢复基线，已将这些非本轮交付文件恢复到提交版本，并保留了真实运行日志。测试临时目录也已清理。

## 调用链与四条 warning

触发链为：

```text
test_s2_equilibrium_degeneration
  -> stage2_pn.run.run(L=2 um, V_f=0.5, V_r=-1.0)
  -> run() 的 Gummel 对照循环
  -> gummel_solve(V_a=0, max_outer=30, omega=1.0 / 0.5)
  -> solve_poisson_pn(..., clamp=0.3, max_iter=100, max_backtrack=12)
```

生产平衡解 `solve_equilibrium -> newton_solve` 是另一条路径，不调用这个 Gummel 对照循环。

### 1. Armijo 线搜索 warning

触发位置是 `semiconductor_sim/stage2_pn/solver.py:311-335`。Gummel 第二个外层迭代（`omega=1`）给 Poisson 子问题的初始缩放残差为 `9.649102844680032e-5`。未限幅的线性 Poisson Newton 方向为 `|delta|_inf=1.6708041137440302 V`，生产对照路径先以 `clamp=0.3 V` 截断，之后 12 次回溯到 `alpha=2^-12=2.44140625e-4`，单步移动仅 `7.32421875e-5 V`，仍不满足充分下降条件。

独立探针 `verification/s2_clamp_probe.log` 对照了相同输入：

* `clamp=0.3`：100 步、`converged=False`、最终残差 `1.0284941386771807e-4`；
* `clamp=1`：15 步、`converged=True`、最终残差 `5.670100290714211e-16`；
* `clamp=10`：2 步、`converged=True`、最终残差 `5.087909564253258e-16`；
* `damped=False, clamp=0.3`：7 步、`converged=True`、最终残差 `5.802411478665453e-16`。

因此这条 warning 的直接原因是 Gummel 对照的 Poisson 步限幅与固定载流子电荷不匹配，不是 Armijo 判据本身证明了物理解不存在。

### 2. 内层 Newton 100 步耗尽

它与第一条 warning 属于同一个内层调用的后果。`solve_poisson_pn` 只有在 `dmove < 1e-11 V` 时才返回收敛；第二个外层（`omega=1`）100 步后残差从 `9.649102844680032e-5` 变为 `1.0284941386771807e-4`，每步仍为 `7.32421875e-5 V`，所以不能宣称收敛。`omega=0.5` 第二个外层的对应残差为 `3.463034834173071e-5 -> 3.4694871461366527e-5`。

诊断脚本逐步记录了 60 个 Gummel 外层和 5,733 个 Poisson Newton 步。`omega=1` 有 29 个内层调用未收敛、2,833 个 Armijo 失败步；`omega=0.5` 有 28 个内层调用未收敛、2,782 个失败步。这些计数来自独立镜像诊断，不依赖 pytest warning 去重。

### 3. `omega=1` 的 Gummel 外层 warning

`gummel_solve` 在 `solver.py:854-885` 的停止条件是 `dphi < 1e-9 V`。第 30 步的记录为：

```text
dphi = 1.4648437500008882e-2 V
dn   = 9.156556445963915e-1
dp   = 7.868867769483983e4
inner_iters = 100, inner_converged = False
```

这不是容差边界效应。使用独立的精确线性 Poisson 步（去掉内层限幅和 Armijo，仅用于诊断）时，`omega=1` 在第 2 个外层的 `dphi=1.6708 V`，第 3 个外层达到约 `1.97e12 V`，随后继续放大。结果见 `verification/s2_exact_poisson_probe.log`。

### 4. `omega=0.5` 的 Gummel 外层 warning

第 30 步为：

```text
dphi = 7.3242187500002776e-3 V
dn   = 3.1564443643122336e-1
dp   = 3.2749772842730174e-1
inner_iters = 100, inner_converged = False
```

它比 `omega=1` 小，是欠松弛后的二周期模态减弱结果，但仍比 `1e-9 V` 停止阈值大约 7 个数量级。精确 Poisson 外层探针显示 `omega=0.5` 在第三个外层已经达到约 `1.99e4 V`，所以把默认 omega 简单改小也不是有证据的最小修复。完整逐步数据在 `verification/s2_convergence_history.csv`，结构化摘要在 `verification/s2_convergence_diagnostics.json`。

## 生产全 Newton 对照

独立命令：

```text
MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/diagnose_s2_production_newton.py
exit 0; iters=8, converged=True, warnings=[]
```

`verification/s2_production_newton.log` 记录了缩放残差与步长：初始残差约 `1.21e-4`，前 3 步降至 `7.591e-10`，随后在组装舍入底附近；第 8 步 `dmove=1.075e-16 V`、`||F||=3.845e-10`，满足当前代码的“步长 + 残差地板”判据。该结果证明生产平衡路径达到其明确实现的停止条件，但不证明弃用的 Gummel 对照已收敛。

## 是否修改代码

本轮没有实施生产代码修复。原因是根因已经明确，但候选改变各自对应不同含义：提高 Poisson clamp 会改变 Gummel 对照的内层步长策略；去掉阻尼会暴露外层固定点发散；继续增加 100/30 的预算或放宽容差只会延长或掩盖未收敛。当前 `run.py` 已明确把 Gummel 标为 deprecated、主流程使用全 Newton；保留 warning 能让这个数值风险继续可见。

没有新增“修复回归测试”，因为没有声称修复任何生产行为。新增的 `scripts/diagnose_s2_convergence.py` 是审计用独立 tracing/mirror 脚本，不改变导入的生产模块；它本身的真实运行退出码为 0。

## 文件与 SHA-256

本轮新增文件及用途见 `verification/s2_changed_files_manifest.json`。主要文件为：

* `scripts/diagnose_s2_convergence.py`：逐步镜像 Poisson/Gummel 更新并输出结构化证据；
* `verification/s2_convergence_diagnostics.json`：每个 omega 的 30 步外层摘要和代表性内层记录；
* `verification/s2_convergence_history.csv`：5,794 行（含表头）的逐步历史；
* `verification/s2_convergence_test.log`、`s2_targeted_tests.log`、`s2_kernel_closure_tests.log`、`s2_full_pytest.log`、`s2_reproduce_phase1.log`：真实命令输出和退出码；
* `verification/s2_clamp_probe.log`、`s2_exact_poisson_probe.log`、`s2_production_newton.log`：根因与生产路径对照。

S1/S2 golden 的当前状态仍为缺失；唯一生产 golden `stage3.json` 未被写入，SHA-256 前后相同。没有删除原始数据或重写历史。

## Git 与后续建议

本轮交付应只存在于 `phase2a/s2-convergence-audit-2026-10-28`。`main` 未修改，也没有重复尝试被拒绝的 GraphQL PR 创建。由于 warning 仍存在，状态为 **PARTIAL**，不建议把本轮标为 PASS 或进入 S1/S2 golden 恢复。

下一步应先由人工决定 Gummel 对照的研究目的：若它只用于展示失败模式，可将报告和测试说明明确标为“预期非收敛诊断”；若要把它变成可用求解器，需要单独设计带线搜索的全耦合/Anderson 固定点算法，并用电流守恒、质量作用和独立参考解验证，不能在本轮通过增大预算或隐藏 warning 处理。
