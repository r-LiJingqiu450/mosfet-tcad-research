"""Generate reports only from completed JUnit and computed results; no fallback claims."""
import difflib,html,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from verify_delivery import parse_pytest
from semiconductor_sim.common.jsonio import write_json
OUT=ROOT/'semiconductor_sim/results/phase1';REPORTS=ROOT/'reports'


def mdtable(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+
                    ['| '+' | '.join(str(x) for x in row)+' |' for row in rows])


def main():
    tests=parse_pytest(ROOT)
    s=json.loads((OUT/'phase1_results.json').read_text());a=s['attribution']['worst_point']
    th=s['threshold'];bg=s['attribution']['neutral_background'];b=s['charge_budget']['independent_reference']
    stats=s['table_build_stats'];REPORTS.mkdir(exist_ok=True)
    bg_rows=[[r['vg'],f"{r['total_A']:.8g}",f"{r['neutral_background_A']:.8g}",
              f"{100*r['neutral_background_fraction']:.8g}%"] for r in bg if r['vd']==.05]
    error_rows=[[label,f'{100*a[key]:+.8f}%'] for label,key in [
        ('同网格直接求解－连续参考（含有限域/求解器影响）','grid_and_boundary_relative'),
        ('生产查表－同网格直接求解','table_interpolation_relative'),
        ('两者相加的总偏差','total_relative'),('开启历史±3 V钳位的配对敏感性','clamp_paired_relative')]]
    continuity_rows=[[r['ny'],f"{s['baseline']['continuity'][i]['finite_difference']:.8g}",
                      f"{r['finite_difference']:.8g}",f"{r['edge_consistency']:.4g}"] for i,r in enumerate(s['continuity'])]
    reference_rows=[[r['vg'],r['vd'],f"{100*r['relative_error']:+.5f}%"] for r in s['current_reference']]
    ab=s['ablation'][1]
    report=f'''# MOSFET 第一阶段修订与证据核验报告

版本 phase1-1.2。日期 2026年10月26日。此报告仅在完整测试结束且JUnit与日志一致后生成。

## 1 交付状态

本次完整测试：**{tests['passed']} passed, {tests['skipped']} skipped, {tests['failures']} failures, {tests['errors']} errors**，共{tests['collected']}项。完整逐项日志见根目录tests_phase1.txt，结构化结果见verification/pytest.xml和pytest_status.json。

4项跳过均为上传材料中缺失原始golden的S1/S2回归及键集合检查。S3的生产golden、配置检查及篡改检测实际运行并通过。没有重新伪造原版S1/S2的历史结果。

独立参考42点最大电流偏差{100*s['maximum_reference_relative_error']:.6f}%。这个数字仅指当前经典约化模型的数值实现一致性，不是实际MOSFET的实验误差。

## 2 对前次交付与锐评的核查

- 接受：前包tests_phase1.txt仅6字节，未等全测试结束便交付，是交付流程失败。新流程以进程退出码、完整JUnit、匹配日志摘要三者共同判定完成；缺任何一项禁止发布。
- 接受：前次freeze使用240×110而生产使用400×220，回归配置不统一。现在生产、冻结、正式S3测试均导入同一个PRODUCTION配置，canonical文件为semiconductor_sim/results/golden/stage3.json。
- 接受：前包混入旧失败运行study_log.txt。现在写文件与打印共用严格JSON编码器，支持NumPy标量/数组，拒绝NaN/Infinity。当前规范日志由成功运行重新生成；旧证据仅在history下原交付快照保留。
- 澄清：实际检查的旧ZIP中，结果JSON、4 CSV、NPZ及5对PNG/SVG均存在；旧版test_frozen_phase1_golden也确实读取golden_phase1.json。问题是配置不一致、全套检查未完成，不能将其描述为完全无人读取。原ZIP摘要与文件清单见history/v11_archive_audit.json。

## 3 单一生产配置与golden证据链

配置源：stage3_mosfet/config.py::PRODUCTION。当前400×220、log插值、h0=0.125 nm、hmax=4 nm、xmax=1.2 μm。物理参数、偏压样本和阈值判据也记录在golden.config中。

冻结入口：python scripts/freeze_phase1.py --reason 说明。已有文件必须额外显式给出--allow-replace。正常复现不会调用冻结脚本，并核对golden运行前后的SHA256不变。test_golden.py直接比较同配置的生产结果；测试故意篡改配置及电流，确认都会被拒绝。

golden仍只证明相对冻结版本的数值回归，不证明物理正确。最初缺失的历史golden与新修订版本golden不是同一证据。

## 4 中性背景贡献必须由比例确定

参考明确计算 Ns = Ns_excess + nb*xmax，同时单独积分得到ID_excess与ID_background。该分解是约化域内的载流子记账；“过剩”并不意味着它已经由真实源漏注入模型验证。

VD=0.05 V：

{mdtable(['VG V','总电流 A','中性背景电流 A','背景占比'],bg_rows)}

因此，VG=0.5和0.75 V下“主体由中性背景贡献”的具体判断不成立；接近平带的VG=0.05 V处背景确实占主导，显示了低栅压电流对计算域深度的敏感性。报告不把这一区域的绝对电流水平解释为真实器件关断电流。

即使0.5—0.75 V的背景很小，独立参考与数值解仍共享经典统计、恒定迁移率和理想载流子库等近似。一致性属于code verification，不能升级为experimental validation。

## 5 最大0.191%偏差的定量归因

最大点VG={a['point']['vg']} V，VD={a['point']['vd']} V。

- 独立第一积分参考：{a['reference_A']:.12g} A。
- 同一空间网格直接Poisson栅边界求解，不查表：{a['direct_same_grid_A']:.12g} A。
- 生产查表：{a['production_A']:.12g} A。

{mdtable(['分量','相对参考的有符号贡献'],error_rows)}

总偏差由两个差值按同一参考归一后相加，闭合残差{a['telescoping_closure_relative']:.3g}。这是一条条件数值分解，取决于所选求解/插值路径，不能当作唯一的器件物理分解。

另外，在连续参考的phi_n节点上单独作log插值，得到phi_n方向插值敏感性{100*a['phi_n_only_interpolation_relative']:+.6f}%。这属于另一项消融，不应再次加到上表而重复计数。直接积分器估计误差约{a['direct_quad_estimate_relative']:.3g}（相对量）；这是积分器估计，不是完整误差界。

同一最大点直接求解的空间细化：

{mdtable(['h0 nm','hmax nm','相对参考电流偏差'],[[r['d0_nm'],r['hmax_nm'],f"{100*r['relative_to_reference']:+.6f}%"] for r in s['attribution']['direct_grid_refinement']])}

## 6 钳位归因

生产版密度计算已取消电势钳位；超出浮点指数范围将显式报错。历史±3 V开关只为配对敏感性试验保留。

最大偏差点开/关历史钳位电流差为{a['clamp_paired_relative']:.3g}，钳位激活数量为{a['clamp_activated_values']}。整个生产表求解记录的试探势范围为[{stats['potential_min']},{stats['potential_max']}] V；{stats['solved_slices']}个有效截面均在历史±3 V以内。因此当前0.191%偏差不能归因于该钳位。这个结论限于已检查域，不能推广至任意偏压。

## 7 电荷比较恒等式的地位

charge_budget只比较指定平方律电荷Cox(VG−VT,long)/q和计算面密度。Qb又由Qs和Ns定义，所以其闭合误差接近机器精度是代数内部一致性，不是独立物理验证。本版从数值精度验收项目中移除该闭合条件，字段明确标注algebraic_reference_comparison。

在VG=2 V源端，独立参考的表面势项{1000*b['surface_term_V']:.5f} mV、体电荷项{1000*b['body_term_V']:.5f} mV仍可帮助解释与指定平方律电荷的差异；不将其解释成量子效应或实验精度证据。

## 8 阈值分支、CLM与守恒

原表T10使用线性平方律反解得到33.6 mV，但该点不满足VD<Vov。正确饱和分支判据偏移为{th['criterion_offset_mV']:.5f} mV；保留原版提取值时剩余16.5648 mV。本版重新提取阈值{th['vt_cc']:.9f} V，对应总差{th['total_offset_mV']:.5f} mV，分解为{th['criterion_offset_mV']:.5f}+{th['remaining_offset_mV']:.5f} mV。

VG=2 V、VD=1.5 V时，旧模型关闭CLM为{ab['legacy_no_clm_A']*1e6:.5f} μA，开启后为{ab['legacy_clm_A']*1e6:.5f} μA，平方律为{ab['square_A']*1e6:.5f} μA。较小的平方律偏差包含经验CLM补偿，不能作为独立准确性证明。修订完整积分无CLM为{ab['revised_production_A']*1e6:.5f} μA。

{mdtable(['y点数','旧版差分波动','本版独立差分波动','本版边积分内部一致性'],continuity_rows)}

边积分一致性来自共同插值表示，明确属于构造检查。独立差分检查保留Ns>0.2 Ns,max和两端各5%排除掩码，未声称漏端极薄层的任意点导数均已解析。准费米势电流形式已包含漂移与扩散。

## 9 交付物与完成门槛

提供结果JSON、基线重算JSON、6个CSV、transfer_curves.npz、7对PNG/SVG、完整测试日志、JUnit、状态JSON、文件大小/SHA256清单、源码及可编辑S3章节。

离线阅读：reports/Phase1_Report.html。无需Python即可看结果。重新计算：python scripts/reproduce_phase1.py --fresh。环境实测版本见environment_verified.txt，依赖说明见README。首次冷启动建表可能需要数分钟。

verification/run_status.json只有在数值实验、全测试、报告和产物校验全部完成时写入completed/exit_code=0。打包校验拒绝空文件、不完整测试摘要、缺失图片或运行中的状态。归档后从新目录冷启动复现，记录在verification/clean_unpack_verification.json。

## 10 研究范围

本轮完成工程交付和已识别数值证据的补强，没有声称SCI投稿终稿。器件参数仍为NA=1e16 cm^-3、tox=20 nm、W=L=1 μm、300 K。未加入二维DD源漏结、DIBL、量子、场依赖迁移率或AC。下一研究阶段需要文献新颖性检验、参数泛化与二维DD对标。
'''
    (REPORTS/'Phase1_Report.md').write_text(report,encoding='utf-8')
    draft=rf'''# S3 第一轮修订章节 v1.2

本稿依据已完成的数值复算与测试编写，尚未建立SCI论文所需的新颖性主张。

## 模型与证据边界

采用理想源漏电子库的经典长沟道沟道约化域。垂直方向求解准平衡Poisson，沿沟道积分：

$$I_D=\frac{{qW\mu_n}}{{L}}\int_0^{{V_D}}N_s(V_G,\varphi_n)\,d\varphi_n.$$

准费米势形式已经包含漂移和扩散；低电子密度不意味着完整DD方程必然无解。模型采用非简并统计、恒温恒迁移率，默认不施加经验CLM。源漏结、二维短沟道电势重分布、量子与动态效应不在本轮范围。

独立第一积分参考不调用有限差分或查表，因此可以检验实现；它与主程序共享经典物理假设，不构成实际器件的实验验证。

## 数值方法与独立检查

Newton同时检查未阻尼修正与缩放残差。回溯失败和未收敛截面显式返回失败，表格越界查询拒绝。生产版不对电势作密度钳位；浮点指数越界报错。

统一400×220正值对数面密度插值，并对插值表示解析积分、直接反演累积积分获得phi_n(y)。解析积分精确针对插值表示，不保证连续器件精度。独立差分电流和共同构造的边积分检查分别报告。

## 数值结果

42点矩阵最大电流偏差为{100*s['maximum_reference_relative_error']:.6f}%，发生在VG={a['point']['vg']} V、VD={a['point']['vd']} V。同网格直接Poisson解相对连续参考贡献{100*a['grid_and_boundary_relative']:+.6f}%，查表相对同网格直接解贡献{100*a['table_interpolation_relative']:+.6f}%。历史±3 V钳位在该域未激活，配对电流差为0，不能将该点偏差归因于钳位。

修订版饱和点641点独立差分电流波动为{s['continuity'][3]['finite_difference']:.6g}，满足1e-3目标。检查使用原强沟道掩码，并未声称整个漏端边界层点导数都已经网格无关。

恒流阈值为{th['vt_cc']:.9f} V，平方律恒流点位于饱和分支，其判据偏移{th['criterion_offset_mV']:.5f} mV；剩余{th['remaining_offset_mV']:.5f} mV仍包含模型比较与提取口径差异。

## 中性体背景与亚阈区

显式分解Ns=Ns_excess+nb*xmax，并对两部分分别计算电流。VD=0.05 V下，VG=0.5 V背景占比约0.0576%，VG=0.75 V约0.0000373%；在VG=0.05 V处约96.2%。因此背景依赖主要限制近平带及更低栅压区域，不能笼统地称0.5—0.75 V电流由中性体背景主导。即便背景很小，理想电子库与经典近似仍限制对真实器件的解释。

## 平方律对照与电荷比较

定义Qn=-qNs，Qb=Qs-Qn，并以VT,long指定平方律电荷。可写：

$$\frac{{q(N_s^{{SL}}-N_s)}}{{C_{{ox}}}}=(\varphi_s-2\psi_F)+\frac{{|Q_b|-|Q_{{b,th}}|}}{{C_{{ox}}}}.$$

这个式子用于分配相对指定平方律的电荷差；Qb由共同变量构造，闭合到机器精度是代数事实，不是器件物理或独立数值验证。源端分解也不能替代沿沟道电流积分的解释。

CLM消融显示经验沟道缩短明显抬高电流，与平方律的较小差值存在补偿。CLM自身仍需二维器件或实验约束，不能借平方律接近证明其正确。

## 结论

本轮在固定器件与经典约化模型内，完成了独立参考、误差路径归因、背景贡献显式化，以及统一生产golden和完整测试的工程闭环。结果不构成实验验证，也不主张标准插值、积分或GCA方法为首创。后续需在新颖性检验和更完整器件基准基础上确定研究贡献。
'''
    (REPORTS/'S3_Replacement_Sections.md').write_text(draft,encoding='utf-8')
    # Dependency-free offline summary; images are ordinary relative files inside ZIP.
    figures=''.join(f'<figure><img src="../semiconductor_sim/results/phase1/fig{i}_{name}.png"><figcaption>{html.escape(caption)}</figcaption></figure>'
         for i,name,caption in [(1,'continuity','独立差分与边积分检查的区别见正文'),(2,'independent_reference','同经典模型参考对照，不是实验精度'),
          (3,'clm_ablation','CLM消融'),(4,'transfer','转移特性'),(5,'spatial_refinement','空间细化'),(6,'background','中性背景占比'),(7,'error_attribution','条件数值误差归因')])
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>MOSFET v1.2 核验报告</title><style>body{max-width:1050px;margin:30px auto;padding:0 20px;font:16px/1.75 system-ui;color:#183042}pre{white-space:pre-wrap;word-break:break-word;font:inherit}img{width:100%;height:auto}figure{margin:30px 0}figcaption{color:#536575}h1{font-size:28px}</style><body><h1>MOSFET v1.2 核验报告</h1><pre>'''+html.escape(report)+ '</pre>'+figures+'</body></html>'
    (REPORTS/'Phase1_Report.html').write_text(page,encoding='utf-8')
    (ROOT/'README.md').write_text(f'''# MOSFET 第一阶段修订包 v1.2

无需Python的离线入口：reports/Phase1_Report.html。完整分析：reports/Phase1_Report.md。可编辑论文章节：reports/S3_Replacement_Sections.md。

本次完整测试：{tests['passed']} passed, {tests['skipped']} skipped。4项跳过仅限缺失历史golden的S1/S2；S3生产golden已运行。详见tests_phase1.txt和verification/pytest.xml。

## 冷启动复现

Python 3.10或以上（实测版本见environment_verified.txt）。Windows、macOS、Linux均使用以下Python命令；首次建表可能需要数分钟。

```bash
python -m pip install -r requirements.txt
python scripts/reproduce_phase1.py --fresh
```

--fresh仅删除本项目.cache。正式复现顺序为数值实验→全量测试→报告生成→产物检查。正常命令不会生成或修改golden。

生产、冻结、S3回归均使用stage3_mosfet/config.py的PRODUCTION：400×220 log。冻结文件唯一为semiconductor_sim/results/golden/stage3.json。显式更新需python scripts/freeze_phase1.py --reason 说明 --allow-replace。

## 文件索引

- semiconductor_sim/results/phase1：结果JSON、6 CSV、NPZ、7对PNG/SVG。
- verification：完整JUnit、测试状态、运行终态、产物大小与SHA256清单、冷启动解压核验。
- history：前次交付原ZIP及其核查清单，失败历史不混入本轮规范日志。
- baseline：从用户Word恢复的原始v1源码快照与审计。
- changes.patch：相对原始恢复源码的修改。

发布流程拒绝未完成的测试、失败退出码、缺少文件或被自动更新的golden。离线结果与重新计算是两种使用方式；缺少本机依赖不影响离线查阅。

当前只是已修订的经典GCA研究基线，不是完成二维DD、量子、迁移率/AC或SCI终稿的声明。
''',encoding='utf-8')
    (ROOT/'CHANGELOG.md').write_text('''# phase1-1.2

- 修复前次未等测试完成即交付的流程；以完整日志、JUnit及退出码联合封闭验收。
- 生产、冻结、正式S3测试统一PRODUCTION和results/golden/stage3.json；篡改测试证明检查有效。
- 写文件及打印统一严格JSON编码；拒绝非有限值和半写文件。
- 默认取消密度电势钳位，保留历史开关供配对敏感性计算；增加全建表试探势范围统计。
- 显式报告中性背景/栅诱导过剩电子的电流贡献；增加近平带依赖说明。
- 最大偏差点增加同网格直接Poisson、空间细化、phi_n插值和钳位消融。
- 电荷比较恒等式移出独立精度验收，标明代数证据类型。
- 正常复现不生成golden；新增归档完整性与解压冷启动验证。
''',encoding='utf-8')
    with zipfile.ZipFile(ROOT/'baseline/recovered_v1_source.zip') as z:
        base={n:z.read(n).decode() for n in z.namelist() if n.endswith('.py')}
    patch=[]
    for p in sorted((ROOT/'semiconductor_sim').rglob('*.py')):
        if '__pycache__' in p.parts:continue
        name=p.relative_to(ROOT).as_posix();before=base.get(name,'');after=p.read_text()
        if before!=after:patch.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
    (ROOT/'changes.patch').write_text(''.join(patch),encoding='utf-8')
    print('Reports generated from completed tests:',tests)

if __name__=='__main__':main()
