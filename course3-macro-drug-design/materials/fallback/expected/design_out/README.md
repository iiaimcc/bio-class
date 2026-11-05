# design_out：教师预生成的 binder 设计产物（已生成）

本目录是课程三任务卡 2–4 的**红灯兜底**：在线 RFdiffusion / ProteinMPNN / 回折排队或不可用时，学生用这里的同参数产物继续完成"参数理解 → 候选筛选 → 结果复核"。**产物已在本机构建完成**（2026-11-05），不是占位。

## 1. 用什么生成的（必须对学生说明）

课程主线是**在线** RFdiffusion（生成 backbone）+ ProteinMPNN（设计序列）+ ColabFold/Boltz（回折）。本目录的兜底产物用本机安装的 **BoltzGen 0.1.2**（`app-boltzgen` conda 环境）一次跑完等价流程：

| 课程主线步骤 | 本目录对应产物 | BoltzGen 步骤 |
|---|---|---|
| RFdiffusion 生成 backbone | `backbones/design_*.cif` | `design`（扩散模型按表位约束生成骨架） |
| ProteinMPNN 设计序列 | `sequence_candidates.csv` 的 `designed_sequence` | `inverse_folding` |
| ColabFold/Boltz 回折 | `refolded/refold_*.cif` | `folding`（Boltz-2 复折叠复合物） |
| 指标与筛选 | `all_designs_metrics*.csv` | `design_folding / affinity / analysis / filtering` |

**不能说成同一个工具**：产物记录里写"本机 BoltzGen 预生成"，学生自己跑在线服务时写"RFdiffusion/ProteinMPNN 在线"。两者是同一类流程（扩散骨架 → 逆折叠 → 复折叠），但不是同一个软件。

## 2. 运行记录（两次运行，两个长度档）

| 项 | 主档（学生规格） | 长档（对照） |
|---|---|---|
| 设计规格 | `spec_egfr_diii_binder_main.yaml`：长度 **60..90** | `spec_egfr_diii_binder_long.yaml`：长度 **80..110** |
| 靶标与约束 | `../inputs/EGFR_DIII_1YY9.pdb`，22 个表位残基作为 `binding` 约束 | 同左 |
| 运行量 | `--num_designs 150`（142 条完成） | `--num_designs 100`（100 条完成） |
| 工具与环境 | BoltzGen 0.1.2，`--protocol protein-anything`，默认 step/noise 调度；`asglibc239` 容器 + V100S 单卡 | 同左 |
| 用时 | 约 2 小时 | 约 1 小时 40 分 |
| 日期 | 2026-11-05 | 2026-11-05 |

> `binding:` 用的是**链内相对编号**（1–171），不是 PDB 编号（310–480）——这是学生最常踩的坑之一，课前要提醒。

## 3. 选出的八个候选（4 成功 + 4 淘汰）

主档（60–90 aa）：

| 设计 | 长度 | 设计→靶标 ipTM | 界面 pAE | 复折叠自洽 RMSD | 氢键 | 与 22 个表位残基重叠 | 角色 |
|---|---|---|---|---|---|---|---|
| `main_082` | 61 | **0.633** | **2.77** | 0.70 Å | 3 | 15 / 22 | 成功（综合最好） |
| `main_078` | 87 | 0.566 | 4.20 | **0.48 Å** | 4 | 18 / 22 | 成功（最自洽） |
| `main_020` | 70 | 0.541 | 4.20 | 0.44 Å | **6** | 16 / 22 | 成功（界面化学最好） |
| `main_101` | 69 | 0.252 | 12.45 | 1.20 Å | 0 | 5 / 22 | **淘汰**：能复折叠但没落在表位上（ipTM 低、无氢键、表位重叠差） |

长档（80–110 aa）：

| 设计 | 长度 | ipTM | pAE | 自洽 RMSD | 氢键 | 表位重叠 | 角色 |
|---|---|---|---|---|---|---|---|
| `long_14` | 85 | 0.597 | 4.03 | 0.66 Å | 5 | 17 / 22 | 成功 |
| `long_29` | 98 | 0.551 | 4.01 | 0.85 Å | 4 | 19 / 22 | 成功（筛选排名第 1） |
| `long_86` | 82 | 0.502 | 5.08 | 0.71 Å | **9** | 20 / 22 | 成功（界面最大） |
| `long_99` | 80 | 0.153 | 16.47 | **13.77 Å** | 0 | 13 / 22 | **淘汰**：设计骨架与复折叠结构对不上，界面结论不可用 |

两个长度档合计 242 条设计里，同时满足"复折叠自洽（RMSD < 2 Å）+ 界面可信（ipTM > 0.4）"的只有 15 + 5 条；ipTM 中位数约 0.2，复折叠 RMSD 中位数约 1.1 Å。**这正是课上要让学生自己看出来的分布**：有输出 ≠ 成功，成功是少数。

## 4. 文件

| 文件 | 内容 |
|---|---|
| `spec_egfr_diii_binder_main.yaml` / `…_long.yaml` | 两次运行的设计规格（可直接用于 `boltzgen run`，或作为 RFdiffusion hotspot 参数的对照） |
| `backbones/design_main_*.cif`、`design_long_*.cif` | 8 个设计骨架（复折叠前） |
| `refolded/refold_main_*.cif`、`refold_long_*.cif` | 8 个复折叠复合物（链 A = binder，链 B = EGFR DIII，编号 1–171 = PDB 310–480） |
| `sequence_candidates.csv` | 8 条设计序列 + 关键指标 + 表位重叠 |
| `all_designs_metrics_main.csv` | 主档全部 142 条设计指标 |
| `all_designs_metrics.csv` | 长档全部 100 条设计指标 |
| `interface_summary.json` / `interface_summary_main.json` | 界面残基与表位重叠的机器可读结果 |
| `analyze_designs.py` / `analyze_main.py` | 复算脚本（PyMOL） |
| `interface_long14.png` / `interface_main082.png` | 两个成功案例的界面图（灰=EGFR DIII，橙=binder，红=落在西妥昔单抗表位上的残基） |

## 5. 课堂用法

1. **参数理解**：把 `spec_*.yaml` 的表位残基与自己的界面约束表逐条对照（编号空间不同：这里是 1–171，PDB 是 310–480）；
2. **候选筛选**：给出八个候选的指标，让学生按"复折叠自洽 → 表位重叠 → 界面化学（氢键/SASA）"排序，说明为什么 `main_082` 综合最好、`main_101` 和 `long_99` 必须淘汰（两种淘汰原因不同：一个没落在表位、一个根本没复折叠）；
3. **结果复核**：叠合 `refolded/refold_main_082.cif` 与 `backbones/design_main_082.cif`，量 RMSD、复算界面残基，与 1YY9 的表位对比；
4. **长度档对照**：同一条表位约束，60–90 与 80–110 两档的成功率与界面大小有差别——请学生解释为什么。

## 6. 边界

- 这些是**计算产物**：ipTM / pAE / RMSD 都不能写成 KD、亲和力、阻断能力或药效；
- 成功案例只是"计算上自洽"，离真实结合还差表达、纯化、SPR/BLI、细胞功能与选择性实验；
- 设计针对 **EGFR 胞外 domain III 的结合剂**，与课程一发现的胞内激酶域激活突变不是一回事（见 `docs/文献依据.md` 第 4 节）。
