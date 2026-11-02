# Boltz 预生成结果说明（1YY9参照复核）— 教师填写模板

> 本文件是可选的教师参照材料，不是课程三的生成式设计主线。课程三主线是 RFdiffusion backbone → ProteinMPNN 序列 → 回折复核。
> 若教师实际运行 Boltz 复核 1YY9 参照复合物，在课前填写真实结果；**禁止编造数值**。
> 当前所有数值栏均为待填占位：本课程包**没有**随附 Boltz 运行结果，教师未跑通时不要把这些栏位讲成已有数据。

## 1. 运行信息（教师填写）

| 项 | 值 |
|---|---|
| 运行命令 | `boltz predict boltz_egfr_cetuximab.yaml --use_msa_server --recycling_steps 3 --diffusion_samples 5 --max_parallel_samples 1 --output_format pdb --out_dir boltz_out` |
| MSA 模式 | 在线 / 离线兜底（`--msa 2 --msa-depth 24`），圈出实际用的 |
| 运行日期与环境 | （如 app-boltz, NVIDIA A10, 用时约 x min） |
| 输出位置 | `materials/fallback/expected/boltz_out/` |

## 2. 基线结果（教师从 `boltz_out/` 运行摘要填入）

| 指标 | 数值 | 一句话解读（供学生对照） |
|---|---|---|
| iptm | 教师填写 | >0.8 说明预测自信，但不等于活性 |
| pAE 图 | 截图引用（见课程三讲义中回折判读部分） | 低 pAE 区域 = 预测稳定的区段 |
| 置信区域 | 教师填写 | 哪个结构域/哪些链最可靠 |
| 与 1YY9 叠合 RMSD（EGFR 胞外域公共段） | 教师填写 | 量级含义 |

叠合图：`fallback/expected/boltz_out/align_vs_1YY9_global.png`（全局）、`fallback/expected/boltz_out/align_vs_1YY9_interface.png`（界面放大）。

## 3. 参照结果使用边界

- 1YY9 的真实 PDB 链为 A=EGFR、C=Fab light、D=Fab heavy；Boltz 输入 ID 可以不同，但必须保存映射；
- ipTM/pAE/RMSD 只能用于结构复核和模型比较，不能写成 KD、阻断能力或药效；
- 新 mini-protein binder 的结果应记录在 `fallback/expected/design_out/` 和学生的 `refold_review.md`，不要塞入本表；
- 如果服务不可用，学生使用教师预生成结构时必须标明“教师预生成”。