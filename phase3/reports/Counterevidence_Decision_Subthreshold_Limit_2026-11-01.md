# C10–C12 回溯证据与下一研究决策

**日期：** 2026-11-01  
**工作分支：** `phase3/subthreshold-counterevidence-2026-11-01`  
**基线提交：** `bf705f040cc0901a9cbb88ef6e411ea70782500d`  
**本阶段新二维器件仿真：** 0

本报告处理的是对 `F(eta, chi)=Delta L/W_dep` 及其
`F=F_inf+A/chi` 近似的回溯压力测试。它不改写 C8–C12 的历史报告、C10
原始结果或正式 golden。C10、C11、C12 的修正结论优先于早期的正面解释。

## 1. 回溯检查确实复现了反例

输入为：

- 附加文件 `chi_fit_C10_retrospective_check.json`，SHA-256
  `b23885461406fd853f2569a67e7965dfe6c2128b072bfd8db4b8e81d80a65d16`；
- `results/phase3/C10/analysis/c10_result.json`，SHA-256
  `3c3ca80264976168bb86a70ec6c67524b487e7e955e505a831794180546366ee`；
- `results/phase3/subthreshold_limit/subthreshold_limit_summary.json`，SHA-256
  `91d9798881f525933ac3eca1c4333775cc0f623652c7c64076f87d0e42efe07e`。

`phase3/subthreshold_limit/verify_retrospective.py` 读取 C10 的
`all_states`，只选择 `level=0`、请求 `h_w=0.5 um`、无内建电势偏移、且与两
个拟合器件相同的 `N_A/t_ox`。计算定义为

\[
  \chi=\hat\psi_{bi,eff}-\hat\psi,\qquad
  \eta=\frac{\epsilon_{Si}}{\epsilon_{ox}}\frac{t_{ox}}{W_d},\qquad
  F_{obs}=\frac{\Delta L_{source}}{W_d},
\]

\[
  F_{pred}=F_{\infty}+A/\chi,
  \qquad e_F=F_{pred}-F_{obs}.
\]

`Delta L_source` 是 C10 采用的无相消源侧结果；排水侧结果在深亚阈值下
可能违反 Rayleigh 界限，因此没有替代源侧值。脚本选出 10 个状态，并且与
附加 JSON 的每个数值字段逐项相同（最大绝对差为 0）。结果保存在：

- `results/phase3/subthreshold_limit/retrospective/retrospective_check.csv`；
- `results/phase3/subthreshold_limit/retrospective/retrospective_check.json`；
- `logs/phase3/subthreshold_limit/retrospective_check.log`。

代表性结果如下；误差是对 `F_obs` 的相对误差：

| `N_A` (cm⁻³) | `t_ox` (nm) | `V_G` (V) | `chi` | `eta` | `F_obs` | `F_pred` | `e_F/F_obs` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e16 | 20 | 0.488174 | 24.8414 | 0.312882 | 2.012641 | 1.752842 | −12.908% |
| 1e16 | 20 | 0.588534 | 21.8414 | 0.277341 | 1.845218 | 1.716493 | −6.976% |
| 1e16 | 20 | 0.686557 | 18.8414 | 0.251682 | 1.716709 | 1.668569 | −2.804% |
| 1e18 | 2 | 0.488174 | 29.4465 | 0.312882 | 2.174522 | 1.695388 | −22.034% |
| 1e18 | 2 | 0.588534 | 26.4465 | 0.277340 | 2.063164 | 1.687767 | −18.195% |
| 1e18 | 2 | 0.621433 | 25.4465 | 0.267936 | 2.017492 | 1.684828 | −16.489% |
| 1e18 | 2 | 0.686557 | 23.4465 | 0.251682 | 1.934484 | 1.678196 | −13.248% |
| 1e18 | 2 | 0.783500 | 20.4258 | 0.231919 | 1.829223 | 1.665719 | −8.938% |
| 1e18 | 2 | 0.877763 | 17.4465 | 0.216378 | 1.742736 | 1.649181 | −5.368% |
| 1e18 | 2 | 0.971587 | 14.4465 | 0.203512 | 1.667415 | 1.625634 | −2.506% |

全体的最大绝对相对误差为 **22.034%**。源侧电流在 10 个状态中全部满足
界限；排水侧有 9/10 个状态不满足界限，这是 C10 的相消问题，不是拟合误差
的补救手段。

这不是预注册的新实验，也没有重新拟合系数。它证明的是：当前发布的、每个
器件固定一组 `(F_inf,A)` 的 `chi`-only 曲线不能外推到这些 C10 状态。由于
`eta` 随偏压变化，而且 C10 已识别出
`(hat psi, hat psi_bi, t_hat, h_hat, n_face/N_A)` 五个无量纲依赖，失败不等于
所有可能的 `F(eta,chi)` 都被证伪；它足以撤回当前曲线的跨器件预测解释。

## 2. 三个不同层次的命题必须分开

### 2.1 局部经验拟合

`F_inf+A/chi` 是在 C8 的两个点（`alpha=15,50`）上拟合、在同一器件的
`alpha=150` 上作内部回代的经验插值。它没有给出 `F_inf(eta)` 或 `A(eta)` 的
生成规则，也没有处理 C10 的 `psi_hat`、结面屏蔽和几何深度依赖。因此它不应
被称为电子无关极限或普适定律。原有六点和 C10 的差异必须同时保留。

### 2.2 切片电流或表面电荷近似

在固定的二维平衡势和恒定迁移率假设下，C11 给出了两个不同网络的量：

\[
  G_{rows}\le G_{2D}\le G_{slice},
  \qquad
  \Delta L_{rows}\le\Delta L_{2D}\le\Delta L_{slice}.
\]

`G_slice` 等价于把每个横向切片在纵向短接后得到的片电导积分；`G_rows` 是
去掉垂直连结后的行网络。Rayleigh 单调性使它们成为界限，但除非界限重合，
切片积分不是实际二维网络修正。C11 的观测是多数深亚阈值状态的界限宽度不超过
约 1.2%，但器件 2 在 `hat psi=12` 时约 5.4%；这不能被写成普适误差上界。

### 2.3 实际二维网络的控制渐近结果

若令

\[
  \varepsilon_n=\frac{N_s}{N_A W_d}\to0,
\]

要得到受控的电子无关极限，至少还必须固定或显式控制
`eta`、`hat psi`、`hat psi_bi`、`t_hat`、`h_hat`、`n_face/N_A`、结面边界条件、
`L/W_d`、统计与迁移率模型，并证明二维势在结面和栅边角附近随
`epsilon_n` 统一收敛。远通道的片电子比很小，不能单独保证近结区的电子种群和
屏蔽也按同阶变小。当前资料没有这样的收敛证明，也没有数值证据支持
`O(epsilon_n)` 的统一余项；该阶数不应继续使用。

此外，任何积分表达式的误差必须拆开：

1. **静电反馈误差：** 以电子自由势代替含电子势；
2. **垂直电荷剖面误差：** 以表面片电荷或一维剖面代替二维体积分；
3. **电流重分布误差：** 以切片/行网络代替完整二维电阻网络；
4. **求解与浮点误差：** 包括 C10 排水侧 `1-phi` 相消。

只有在给定的二维平衡势上，切片积分才是网络问题的近似/界限；它不是仅由
`eta,chi` 决定的参数自由预测器，也不等同于完整漂移扩散小信号网络。

## 3. 研究方向决策

| 候选方向 | 现有证据 | 本阶段决定 |
|---|---|---|
| 继续 `F_inf+A/chi` 作为跨器件极限 | 10 个 C10 状态误差 −2.51% 到 −22.03%；系数无生成规则 | **关闭为预测主张**，保留为探索性回代 |
| 增加 `F(eta,chi)` 拟合或第三器件 | 可能吸收偏差，但 C10 已有五个无量纲依赖；逐器件重拟合不能叫跨器件预测 | **暂不投资**；没有能区分新知识与再拟合的判别实验 |
| 运行 C11 建议的双栅/无掺杂第三结构 | C12 已证明零偏长沟道公式是 Liang–Taur (2004) 的特例，纯 `Ein` 只是薄膜极限 | **撤回**；不重复已有理论 |
| 低电流端电流的可审计可靠性证书 | C10 观察到排水侧相消；C11 独立有限体积解在两个代表状态复核到 ≤0.25%；Rayleigh 界限能区分有效/失效计算 | **选择为下一问题的筛选方向**，尚非原创性结论 |

下一阶段应问的具体问题是：

> 能否用 Rayleigh 界限、源/漏电流互易性和独立有限体积残差，构造一个在低电流
> 小信号结果被解释之前就能拒绝数值失真的、可复核的端电流验收流程？

这个问题与 `F(eta,chi)` 不同：它针对可验证的数值可信性，而不是继续寻找一个
被 C10–C12 证据削弱的闭式物理极限。它仍有明显先验工作风险，不能预先称为
SCI 创新。进入实施前需做一次针对加权残差、端电流误差估计和低电流小信号求解的
文献门槛，并在已有状态上预注册判别条件。最小可判别测试应同时包含：

- 含已知相消的 C10 排水侧值；
- 无相消源侧值和交换边界值；
- C11 的独立有限体积结果；
- 至少一组未参与条件设计的留出状态。

只有在不把经验阈值冒充数学上界、并能量化误报/漏报后，才值得判断这是否形成
方法学贡献。

## 4. 结论与可复现命令

**科学结论：** 提议的电子无关 `chi` 极限作为当前固定系数模型不经受 C10
反例；可保留的内容是局部经验坐标、C11 的切片/行网络界限，以及 C10–C11
对相消错误和非线性前渐近区的诊断。尚无受控余项、普适 `F(eta,chi)` 或跨器件
预测器。

**状态：** `PARTIAL — LIMIT_REJECTED_AS_GENERAL_MODEL; NEXT_QUESTION_SELECTED_FOR_SCREENING`

复现回溯结果：

```bash
python phase3/subthreshold_limit/verify_retrospective.py
```

该命令退出码为 0，输出 `n=10`、最大绝对相对误差 `22.034023067506247%`、
源侧界限失败 `0`、排水侧界限失败 `9`。本阶段未运行新的二维器件仿真，未修改
生产 S1/S2/S3、正式 golden 或 C10/C11/C12 历史数据。

相关测试实际运行两次。第一次直接调用 `.venv-devsim/bin/python` 时，DEVSIM
找不到默认名称的 BLAS/LAPACK 动态库，导入退出码为 1；诊断日志为
`logs/phase3/subthreshold_limit/devsim_import_without_math_env.log`。随后使用
系统中实际存在的库路径重跑：

```bash
DEVSIM_MATH_LIBS='/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3:/usr/lib/x86_64-linux-gnu/libblas.so.3' \\
  .venv-devsim/bin/python -m pytest -q semiconductor_sim phase3/tests
```

该命令退出码为 0，结果为 **171 passed, 4 skipped, 4 warnings**（28.96 s）。
警告仍是已有的 S2 Newton/ Gummel 非收敛诊断；没有隐藏或降级。完整日志在
`logs/phase3/subthreshold_limit/combined_pytest.log`。

相关证据索引：

- C10 原始分析与数值问题：`phase3/reports/C10_Research_Phase_Report_2026-10-31.md`；
- C11 独立有限体积、Rayleigh 界限和模型失败：`phase3/reports/C11_Research_Phase_Report_2026-10-31.md`；
- C12 对 Liang–Taur 既有理论的文献门槛：`phase3/reports/C12_MultiGate_Literature_Assessment_2026-11-01.md`；
- 逐点机器可读结果：`results/phase3/subthreshold_limit/retrospective/`；
- 当前 handoff：`PROJECT_HANDOFF.md`。
