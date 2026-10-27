# MOSFET 第一阶段修订包 v1.2

无需Python的离线入口：reports/Phase1_Report.html。完整分析：reports/Phase1_Report.md。可编辑论文章节：reports/S3_Replacement_Sections.md。

本次完整测试：91 passed, 4 skipped。4项跳过仅限缺失历史golden的S1/S2；S3生产golden已运行。详见tests_phase1.txt和verification/pytest.xml。

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
