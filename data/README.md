# 课程数据说明（data/）

本目录存放三门课程共用的示例数据，全部由 `make_synthetic_expression.py` 生成。

## 文件清单

| 文件 | 内容 | 学生可见 |
|---|---|---|
| `expression_matrix.csv` | 表达计数矩阵：10000 基因 × 12 样本（整数 counts） | 是 |
| `sample_meta.csv` | 样本注释：S01–S06 Normal，S07–S12 Tumor；两批次 B1/B2 | 是 |
| `gene_annotation.csv` | 基因注释：symbol、模拟 Ensembl ID、染色体、基因类型（protein_coding/lncRNA）、是否真值差异基因 | 是 |
| `druggable_panel.csv` | 可成药基因面板（247 个常见药物靶点基因，含类别标注） | 是 |
| `ground_truth_deg.csv` | **教师专用**真值表：真正的差异基因与 log2FC（不发给学生） | 否 |
| `make_synthetic_expression.py` | 数据生成脚本（纯 Python 标准库，可重复生成） | 是（供复现教学） |

## 数据设计（教学意图）

- 模拟"肺腺癌（LUAD）vs 正常"研究：12 个样本、约 10000 基因；
- 约 3%（300 个）真值差异基因，其中 EGFR/MAPK 信号通路核心基因（EGFR、KRAS、BRAF、MAPK1、VEGFA、MYC、AKT1 等）为上调主力，少数肿瘤抑制基因（PTEN、CDKN1A 等）下调；
- 真值差异基因中 36 个与可成药面板重合，按本课阈值实际检出约 18 个候选；EGFR 作为三课联动的主靶点被固定为强上调（真值 log2FC = 3.5、高基线表达），不依赖随机抽样是否"抽中"；
- 库深度与批次效应如实模拟（不同样本总 reads 不同、B1/B2 倍数差异），供教学演示"为什么需要归一化"；
- 包含约 1200 个低表达 lncRNA 类基因（部分零计数），贴近真实转录组的"寂静基因"现象；
- 数据为**统计合成**，不包含任何真实患者信息；统计特征参照公开肿瘤转录组文献（TCGA/LUAD）构造。

## 复现

```bash
python3 make_synthetic_expression.py            # 输出到仓库 data/（脚本所在目录，默认值）
python3 make_synthetic_expression.py -o /tmp/d  # 输出到指定目录
```

脚本使用固定随机种子（`SEED = 20260901`），且所有随机抽样都走同一个 `random.Random(SEED)` 实例；在任意机器上重复运行，输出逐字节一致。本目录内的 CSV 就是该脚本的发布版本，可用于期末核验学生是否复现。

**发布版本校验值**（`python3 make_synthetic_expression.py` 后应为同一组值）：

| 文件 | SHA256 |
|---|---|
| `expression_matrix.csv` | `0d5d33c9a4dbcc19dd0fff0c3c1e449b53a2c0216647b0088ba22a6281ca42ba` |
| `sample_meta.csv` | `61131f21b10ffd85d8fa53194a6b103b00c8b0601866220b9fe02ee69213f500` |
| `gene_annotation.csv` | `b717fc8177af3db4d9ee81b82e745b0efc3a6fa90f2a4aead48a969afe441189` |
| `druggable_panel.csv` | `b8e19d1bc2ba67f7b11e93eb94b02b5b337166e31a540feaf2eb9deb23494d6b` |
| `ground_truth_deg.csv` | `48cd067ce75bdd881776285ea7fb78e3978f31c0c194258be60356dfc424567b` |

```bash
sha256sum expression_matrix.csv sample_meta.csv gene_annotation.csv druggable_panel.csv ground_truth_deg.csv
```

> 注意：脚本改动（例如换主题、调整通路基因）会改变输出和所有教学基准。修改后必须重跑三份 R 脚本，并同步更新《实操手册》《教师准备清单》《讲义大纲》《教案》《课堂讲义》中的数字。

## 当前发布版本的教学基准

| 项目 | 值 |
|---|---|
| 显著 DEG（padj<0.05, \|log2FC\|>1） | 116（上调 76 / 下调 40） |
| EGFR（完整差异表，按 padj） | rank 2，log2FC ≈ 3.51 |
| 可成药面板 | 247 个基因 |
| shortlist 候选 | 18 个；EGFR 按 log2FC 排 rank 1 |
| shortlist top5 | EGFR、BRAF、MAPK3、TSC2、PIK3R1 |

## 数据替换指南（换主题）

如需换疾病/靶点主题（例如乳腺癌 EZH2、炎症 JAK2），修改脚本中三处即可：

1. `EGFR_PATHWAY` / `DOWN_PATHWAY` 列表换成新主题的通路基因；
2. 把新主题的核心靶基因加入 `ONCO_SIGNAL` 或新增一个列表并注册到 `assign_differential_genes`；
3. 重新生成数据后，**务必重跑三份 R 脚本**，确认新核心靶点进入候选 top10（先看真值表再定学生预期）。

> 提示：`druggable_panel.csv` 基于常见药物靶点基因整理，换主题时优先从面板内选"核心靶点"（如单抗/生物药已有临床应用记录的靶基因），保证课程三大分子设计案例（抗体-抗原为例）仍可衔接。