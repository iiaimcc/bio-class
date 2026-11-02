# design_out（待教师生成）

本目录是课程三任务卡 2–4 的**红灯兜底**位置，当前为空占位。教师必须在预演周按
`../README.md` 第 4 节 runbook 用教师账号实际生成后放入，否则"在线不可用时使用同参数预生成结果"这条方案没有实物。

## 需要放入的文件

| 文件 | 内容 | 来源 |
|---|---|---|
| `bb_<seed>_<n>.pdb` | RFdiffusion backbone，2–4 个（含至少 1 个淘汰候选） | ColabDesign RFdiffusion |
| `params_rfdiffusion.json` | 目标结构、hotspot、长度、种子、设计数 | 同上 |
| `sequence_candidates.csv` | ProteinMPNN 输出序列与分数（≥3 条） | ProteinMPNN |
| `params_proteinmpnn.json` | 固定链、设计链、采样数、种子 | 同上 |
| `refold_<candidate>.pdb` | 回折结构（≥2 条） | ColabFold 或在线 Boltz |
| `refold_summary.csv` | pLDDT/pAE、链映射、界面残基数、与 1YY9 表位重叠 | 学生/教师计算 |
| `run_log.md` | 工具版本、日期、运行链接或截图、失败记录 | — |

## 输入

- 目标结构：`../inputs/EGFR_DIII_1YY9.pdb`（1YY9 A 链 310–480）
- hotspot：`../inputs/EGFR_DIII_hotspots.txt`（22 个 4.5 Å 表位残基）

## 自检

- [ ] 每个文件能说清"谁在什么时候用什么参数生成的"；
- [ ] backbone 不穿入 EGFR、至少 1 个靠近 hotspot；
- [ ] 回折结构的接触残基落在目标表位附近；
- [ ] `refold_summary.csv` 里的数值来自工具输出，不是估计值；
- [ ] 学生提交物模板要求标注"教师预生成"。
