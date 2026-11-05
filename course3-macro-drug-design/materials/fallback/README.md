# 课程三兜底材料（fallback）

课程三的 1YY9 结构与 P00533 序列由学生**在课上现场下载**；RFdiffusion / ProteinMPNN / 回折属于在线服务，课上现场跑小规模任务。本目录分两部分：

- `inputs/`：与线上同源的原始文件，以及由它们派生的设计输入（domain III 目标结构、表位残基清单）；
- `expected/`：任务卡 1 的参考输出（链事实、界面接触表、CDR 区间），以及 `design_out/`、`boltz_out/` 的**生成 runbook 与填写清单**。

> **`design_out/` 与 `boltz_out/` 已经用实物填好**（2026-11-05 本机生成，见各自目录的 README）：
> - `design_out/`：BoltzGen 0.1.2 生成的 4 个 backbone + 设计序列 + 复折叠复合物 + 全部 100 条设计的指标 + 界面复核（含 1 个成功案例与 1 个失败案例）；
> - `boltz_out/`：Boltz-2 对 1YY9 参照复合物的复核（iptm 0.972、与实验结构叠合 RMSD 0.321 Å）。
>
> 两者都标注了工具、版本、参数与日期。第 4 节的 runbook 仍有意义：如果你希望课上让学生用**在线 RFdiffusion/ProteinMPNN** 跑出自己的产物，它就是那套步骤；本目录的实物用于"服务不可用/排队过久"时的兜底，以及课前预演比对。

## 1. 输入文件

抓取日期：2026-11-02。

| 文件 | 来源 | 说明 | SHA256 |
|---|---|---|---|
| `inputs/1YY9.pdb` | https://files.rcsb.org/download/1YY9.pdb | cetuximab Fab–EGFR 胞外域复合物，X-ray 2.61 Å | `e574b33a287c4100ff35f8346a6eadba1486d4a90464f4344f3189115f6d5fa8` |
| `inputs/P00533.fasta` | https://rest.uniprot.org/uniprotkb/P00533.fasta | EGFR canonical 1210 aa | `684e2da1f3ef05d0ce90d0aefe9125c6f7e9cdcfeac27ae72ed35243b64a73e8` |
| `inputs/EGFR_DIII_1YY9.pdb` | 由 `1YY9.pdb` A 链 310–480 提取 | design 的靶标结构；PDB 编号，未重新编号 | `15b2fc071b681add535ffac53a42cb42af0e131f5770842630146619f0a2cbf9` |
| `inputs/EGFR_DIII_hotspots.txt` | 由 4.5 Å 界面接触计算 | 22 个表位残基，供 RFdiffusion hotspot 参数使用 | `775f9b07a7d2eba3d05c295948ab21a1025cabc9916bcb5c8e2829c3c97adfc6` |

`EGFR_DIII_1YY9.pdb` 只保留 A 链 310–480（171 个残基），坐标、链名、编号与 1YY9 一致，便于与原始复合物直接 `align`。它不是"清洗过的设计级结构"：水、糖基、缺失残基都没有补；需要更干净的目标时，按 runbook 自行处理并记录。

## 2. 使用规则

1. 学生先现场下载 1YY9 与 P00533，记录 URL、访问日期、`sha256`；
2. 只有下载失败时才用 `inputs/` 副本，并在产出中写明"使用兜底文件 + 失败原因"；
3. 在线设计/回折失败时使用教师预生成结果，必须标注"教师预生成"及工具、日期、参数；
4. 不允许把兜底输入或教师结果写成自己现场获取/运行。

## 3. 课程三事实基线

- 1YY9：2.61 Å；蛋白链 A=EGFR 胞外域（PDB 1–618，对应 UniProt 25–642）、C=Fab 轻链（SEQRES 213）、D=Fab 重链（SEQRES 221）；B/E/F/G 是 N-糖链，不是蛋白链；
- 建模范围：A 2–614（缺 1、615–624），C 1–211（缺 212–213），D 1–220（缺 221）；
- cetuximab 靶向 EGFR domain III（UniProt 334–504 = PDB 310–480）；
- VH 在链 D（1–119 为可变区，之后接 CH1），VL 在链 C（1–107 为可变区，之后接 CL）；
- 4.5 Å 界面：EGFR 侧 22 个残基、轻链 7 个、重链 14 个，共 305 个原子对，最近原子距离 2.60 Å。

## 4. 教师预生成 runbook（`design_out/` 与 `boltz_out/`）

以下步骤在预演周用教师账号执行，产物放入本目录；每一步都要记录工具、版本、日期、参数、随机种子和输入文件哈希。

### 4.1 RFdiffusion backbone → `expected/design_out/`

1. 目标结构：`inputs/EGFR_DIII_1YY9.pdb`；hotspot：`inputs/EGFR_DIII_hotspots.txt` 中的 22 个残基（先选 4–6 个连续表位残基做小规模试跑，再扩到全部）；
2. 长度 60–90 aa，`num_designs` 按课时控制在 2–4 个；固定随机种子并记录；
3. 保存：`bb_<seed>_<n>.pdb`、运行 notebook 链接或截图、参数 JSON；
4. 用 PyMOL 检查 backbone 是否穿入 EGFR、是否靠近 hotspot、是否有明显碰撞，保留至少 1 个失败候选作为教学对照。

### 4.2 ProteinMPNN → `expected/design_out/`

1. 固定 EGFR 链，只设计 binder 链；`num_seq_per_target` ≥ 3；
2. 保存 `sequence_candidates.csv`（列：candidate、sequence、score、length、fixed_chains、designed_chain、seed）；
3. 至少 2 条进入回折，1 条作为淘汰对照。

### 4.3 回折（ColabFold 或在线 Boltz）→ `expected/design_out/`

1. 输入：EGFR domain III + binder 序列；记录服务名、日期、模型数、recycle 数；
2. 保存结构、pLDDT/PAE、链名映射；
3. 用 `pymol_interface_commands.pml` 的接触选择计算界面，结果写入 `refold_summary.csv`。

### 4.4 cetuximab 参照复核（可选）→ `expected/boltz_out/`

1. 用 `boltz_antibody_template.yaml` 填入 `EGFR_cetuximab_序列.md` 中的三条序列；
2. 按模板命令运行，保存 iptm、pAE 图、与 1YY9 的叠合 RMSD；
3. 把真实数值填入 `boltz_预生成_EGFR_cetuximab.md`，不要保留"教师填写"占位。

### 4.5 本机已完成的部分（2026-11-05）

| runbook 步骤 | 状态 | 产物 |
|---|---|---|
| 4.1 backbone 生成 | ✅ 已完成（用 BoltzGen，而非在线 RFdiffusion） | `expected/design_out/backbones/design_long_{14,29,86,99}.cif` |
| 4.2 序列设计 | ✅ 已完成 | `expected/design_out/sequence_candidates.csv`（4 条，含长度、ipTM、RMSD、表位重叠） |
| 4.3 复折叠复核 | ✅ 已完成 | `expected/design_out/refolded/refold_long_*.cif` + `all_designs_metrics.csv`（100 条） |
| 4.4 cetuximab 参照复核 | ✅ 已完成 | `expected/boltz_out/`（iptm 0.972、叠合 RMSD 0.321 Å、界面 21/22 重合） |

> 与在线 RFdiffusion/ProteinMPNN 的对应关系、以及"必须对学生说明工具不同"这一点，写在 `expected/design_out/README.md` 第 1 节。

### 4.6 自检（对照 runbook）

- [x] `expected/design_out/` 有 ≥2 个 backbone、≥3 条序列、≥2 个回折结构，附参数与日期（现为 4 / 4 / 4）；
- [x] `expected/boltz_out/` 有 iptm 数值与叠合图；
- [x] 所有文件能追溯到具体工具版本与参数（见两个子目录 README 的运行记录）；
- [x] 学生任务卡 2–4 的红灯路径能直接取用这些文件，且报告模板要求标注来源。
