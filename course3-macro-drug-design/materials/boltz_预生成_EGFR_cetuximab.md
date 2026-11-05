# Boltz 预生成结果说明（1YY9 参照复核）

> 本文件是教师参照材料，不是课程三的生成式设计主线。主线是"生成 backbone → 设计序列 → 回折复核"（学生用在线 RFdiffusion/ProteinMPNN，教师兜底产物见 `fallback/expected/design_out/`）。
> 本文件的数值**已由教师在本机实际跑出**（2026-11-05），原始输出在 `fallback/expected/boltz_out/`。

## 1. 运行信息

| 项 | 值 |
|---|---|
| 运行命令 | `boltz predict boltz_egfr_cetuximab_msa.yaml --recycling_steps 3 --diffusion_samples 1 --devices 1 --output_format pdb --out_dir out` |
| MSA 模式 | 离线：登录节点用 ColabFold MMseqs2 服务器（`env` 模式）先生成 `.a3m`，再在 GPU 节点读本地 MSA 运行（计算节点无外网 DNS） |
| 运行日期与环境 | 2026-11-05；`app-boltz`（boltz 2.2.1）+ `asglibc239` 容器；V100S 单卡，预测约 1 分钟 |
| 输出位置 | `materials/fallback/expected/boltz_out/`（含结构、置信度 JSON、叠合图、复算脚本） |

## 2. 基线结果（实测，不是占位）

| 指标 | 数值 | 一句话解读（供学生对照） |
|---|---|---|
| iptm | **0.972** | 模型对"这三条链按 1YY9 的方式结合"很自信；自信 ≠ 亲和力，更不能写成 KD |
| ptm / complex_plddt | 0.978 / 0.981 | 单体折叠与复合物整体置信度都高 |
| 链对 ipTM | A–C 0.966、A–D 0.969 | EGFR 与轻链、重链两个界面都被认为可靠 |
| 与 1YY9 叠合 RMSD（EGFR domain III，148 CA） | **0.321 Å** | 预测结构与实验结构几乎重合 |
| VL / VH 叠合 RMSD（CA） | 0.243 Å / 0.229 Å | 抗体两条链也对得上 |
| 4.5 Å 界面残基（EGFR 侧） | 实验 22 个 / 预测 21 个，重合 21 个 | 界面位置基本一致；实验独有 A411 未进入预测界面 |
| 界面原子数（EGFR 侧 4.5 Å 内） | 实验 77 / 预测 70 | 接触密度接近 |

叠合图：`fallback/expected/boltz_out/align_vs_1YY9_global.png`（全局）、`align_vs_1YY9_interface.png`（界面放大）；数值的机器可读版本：`compare_vs_1yy9.json`；复算脚本：`compare_vs_1yy9.py`、`render_views.pml`。

## 3. 参照结果使用边界

- 1YY9 的真实 PDB 链为 A=EGFR、C=Fab light、D=Fab heavy；本次 Boltz 输入沿用 A/C/D，与 PDB 一致，叠合前只把预测的 A 链编号整体加 309（PDB 310–480 = 输入 1–171）；
- ipTM/pAE/RMSD 只能用于结构复核和模型比较，不能写成 KD、阻断能力或药效；
- 这份结果是**已知序列**的重建（序列来自 UniProt 与 1YY9 SEQRES），难度低于学生设计全新序列，对比时要说明；
- 新 mini-protein binder 的结果记录在 `fallback/expected/design_out/` 与学生的 `refold_review.md`，不要塞入本表；
- 学生使用教师预生成结构时必须标明"教师预生成"并写清工具、日期、参数。
