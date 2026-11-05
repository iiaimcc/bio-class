# boltz_out：参照复合物的预测复核（已生成）

本目录是课程三"参照复合物"的 Boltz 复核结果，用来和实验结构 1YY9 对照。**产物已在本机构建完成**（2026-11-05），不是占位。

## 这次跑了什么

| 项目 | 值 |
|---|---|
| 工具 | Boltz-2（`boltz` 2.2.1，默认 checkpoint） |
| 输入 | `boltz_egfr_cetuximab_msa.yaml`：EGFR domain III（P00533 334–504）+ 西妥昔单抗 VL（1YY9 C 1–107）+ VH（1YY9 D 1–119） |
| MSA | ColabFold MMseqs2 服务器（`env` 模式，greedy 配对）；合并后 A 3090 / C 10728 / D 8718 条序列 |
| 运行 | `boltz predict … --recycling_steps 3 --diffusion_samples 1 --output_format pdb`，单卡 V100S |
| 日期 | 2026-11-05 |

## 与 1YY9 的对照结果

| 指标 | 值 | 复算方式 |
|---|---|---|
| EGFR domain III（CA）叠合 RMSD | **0.321 Å**（148 个 CA） | `align pred and chain A and name CA, ref and chain A and name CA` |
| 轻链 VL（CA，1–107） | 0.243 Å（89 CA） | 同上，链 C |
| 重链 VH（CA，1–107 区间） | 0.229 Å（94 CA） | 同上，链 D |
| EGFR 侧 4.5 Å 界面残基 | 实验 22 个 / 预测 21 个，**重合 21 个**，预测未多出，实验独有 1 个（A411） | `within 4.5 of` 命令见 `compare_vs_1yy9.py` |
| 界面原子数（EGFR 侧 4.5 Å 内） | 实验 77 / 预测 70 | 同上 |

置信度（`prediction/confidence_*.json`）：`ptm` 0.978、**`iptm` 0.972**、`complex_plddt` 0.981；链对 ipTM：A–C 0.966、A–D 0.969。

## 文件

| 文件 | 内容 |
|---|---|
| `boltz_egfr_cetuximab_msa.yaml` | 实际运行的输入（三条真实序列 + 各自 MSA 路径；链映射 A=DIII、C=VL、D=VH） |
| `prediction/boltz_egfr_cetuximab_msa_model_0.pdb` | 预测复合物结构 |
| `prediction/confidence_*.json` | 官方置信度指标（不得改数） |
| `align_vs_1yy9_global.png` | 预测（彩色）与 1YY9（灰、半透明）叠合总览 |
| `align_vs_1yy9_interface.png` | 界面放大 |
| `compare_vs_1yy9.json` | 上表所有数值的机器可读版本 |
| `compare_vs_1yy9.py` / `render_views.pml` | 复算脚本（PyMOL） |

## 使用规则

1. 这是**参照复合物**的复核：它验证的是"用 Boltz 重建已知抗体-抗原界面"的可靠性，**不是**学生设计的 binder 的结果；
2. 学生自己的 binder 复合物预测结果放在自己的 `boltz_out/`（课程目录外的个人产出），不要与这份混淆；
3. `iptm`/`pLDDT`/RMSD 都是**计算指标**：0.97 的 ipTM 说明模型内部自洽，叠合 0.32 Å 说明这次重建与实验一致，二者都不能写成亲和力、KD 或阻断能力；
4. 与学生的设计对比时，注意本次预测**使用了实验序列**（VH/VL/DIII 都来自 1YY9 与 UniProt），而学生设计的是全新序列——难度不可比。
