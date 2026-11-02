# 课程数据说明（data/）

本目录存放三门课程共用的**真实公开数据**：TCGA-LUAD（肺腺癌）的一个小子集，来自 NCI Genomic Data Commons（GDC）。数据已去标识化、公开可再分发，不含可识别的患者信息。

## 文件清单

| 文件 | 内容 | 学生可见 |
|---|---|---|
| `expression_matrix.csv` | 表达计数矩阵：60660 基因 × 30 样本（GDC STAR-Counts，unstranded 整数 counts） | 是 |
| `sample_meta.csv` | 样本注释：`sample_id, case_id, condition, sample_type, tss`；Tumor 15 / Normal 15 | 是 |
| `gene_annotation.csv` | 基因注释：Ensembl gene_id、symbol、gene_type | 是 |
| `druggable_panel.csv` | 可成药基因面板（247 个常见药物靶点基因，含 `panel_class`、`known_drug_class`） | 是 |
| `mutation_status.csv` | 突变长表：每行一个非沉默变异（gene、variant_class、hgvsp、chr、pos） | 是 |
| `mutation_matrix.csv` | 突变矩阵：14 例有 WXS 的肿瘤 × 14 个驱动基因（0/1） | 是 |
| `sample_manifest.csv` | 30 个表达文件的 GDC file_id / md5 / 样本信息（下载脚本的输入清单） | 是 |
| `mutation_manifest.csv` | 14 个 MAF 文件的 GDC file_id / md5 | 是 |
| `fetch_tcga_luad.py` | 下载 + 校验 + 合并脚本（纯 Python 标准库，联网时重建全部文件） | 是 |

## 数据来源与选择规则

- 项目：TCGA-LUAD（The Cancer Genome Atlas, Lung Adenocarcinoma）；
- 数据类型：RNA-Seq，分析流程 STAR - Counts（`*_star_gene_counts.tsv`），样本类型 Primary Tumor 与 Solid Tissue Normal，各取 15 例，来自不同患者（非配对设计）；
- 突变数据：Masked Somatic Mutation（MAF，`gdc-2.0.0-aliquot-merged-masked`），按同一批肿瘤取全部 MAF 文件；
- 选择规则写在 `fetch_tcga_luad.py --help` 与脚本注释里；`sample_manifest.csv` / `mutation_manifest.csv` 把 file_id 锁死，保证任何人重建得到同一批样本；
- 30 例不是随机抽的：在 TCGA-LUAD 的 STAR-Counts 文件里排除 3 个体积明显异常（23 KB~428 KB，正常约 4 MB）的文件后锁定 30 例，覆盖 15 个不同 TSS、包含 EGFR 激活突变（19del 与 L858R 各 1 例）。

## 为什么用真实数据

- 课上教的每一步（库深度差异、低表达过滤、多重检验、富集背景集、突变 vs 表达）都是真实数据里真实存在的问题，不需要"人为制造麻烦"；
- 结论可以被外部核查：任何人拿同样的 GDC 文件重跑，应该得到同一批数；
- 隐私：TCGA 数据已去标识化并获准公开再分发（GDC Data Use Certification 下的公开层）。学生作品里引用具体样本编号（如 TCGA-05-XXXX）没有问题，但不能写成"我们医院的患者"。

## 复现

```bash
# 完整重建（30 个文件 / 约 128 MB，需要能访问 api.gdc.cancer.gov；本机实测首次约 13 分钟，
# 受出口带宽影响可能更久，机房上课请用共享缓存或兜底副本）
python3 fetch_tcga_luad.py --out . --cache ./gdc_cache

# 只校验已有文件是否与发布版一致
python3 fetch_tcga_luad.py --out . --cache ./gdc_cache --verify-only
```

脚本对每个文件校验 md5（来自 manifest），全部通过后才合并 CSV；已下载且 md5 正确的文件会被缓存复用。
`--verify-only` 只报告脚本重建的 5 个 CSV（expression / gene_annotation / sample_meta / mutation_status / mutation_matrix）；`druggable_panel.csv` 是人工整理的面板、两个 manifest 是输入清单，都不在重建范围内，用下面的 `sha256sum` 逐个人工核对。

**发布版本校验值**（`sha256sum` 应得同一组值）：

| 文件 | SHA256 |
|---|---|
| `expression_matrix.csv` | `1da1bdca79bba6787bc69b9dfb1eafc3431b4a95928697871069342c9109c576` |
| `sample_meta.csv` | `e9826e41c85b9cd7d68613c74c96452e4e18548690f0f2532d92991d155cb1fd` |
| `gene_annotation.csv` | `32c0c37fad1edffe8c0cf138a18ba6937dd52f7dad23e260021d3ef5ed987713` |
| `druggable_panel.csv` | `b8e19d1bc2ba67f7b11e93eb94b02b5b337166e31a540feaf2eb9deb23494d6b` |
| `mutation_status.csv` | `f7caf3facb1c150dff8f48bc7853f1f8bbf785abdd7d5f9d06d63dc90b8d31a2` |
| `mutation_matrix.csv` | `349412bd4b4f59e198ca8f7f4b352880b56862214314bf74fbacbe07527471b8` |
| `sample_manifest.csv` | `98234d369cca5812523f21711fa135860d49409c06b3423d09379fbded30d1b4` |
| `mutation_manifest.csv` | `cde2eeb3d544495c9d775d307d868457c3b392938cda7c04edf3ff610f38eefa` |

```bash
sha256sum expression_matrix.csv sample_meta.csv gene_annotation.csv druggable_panel.csv \
          mutation_status.csv mutation_matrix.csv sample_manifest.csv mutation_manifest.csv
```

> 注意：改样本数、换样本或改注释版本都会改变输出和所有教学基准（DEG 数、富集条目数、候选数）。修改后必须重跑三份 R 脚本，并同步更新《实操手册》《教师准备清单》《讲义大纲》《教案》《课堂讲义》和 [docs/课程包自检结果.md](../docs/课程包自检结果.md) 中的数字。

## 断网兜底

`course1-omics-target-discovery/materials/fallback/inputs/` 存放同一份数据的压缩副本与学生结果参考（含 SHA256），供机房网络不可用时解压使用。详见该目录的 README。

## 教学基准（发布版实测）

| 项目 | 值 |
|---|---|
| counts 维度 | 60660 × 30 |
| 低表达过滤后（≥6 样本 counts≥10） | 23765 |
| 显著 DEG（padj<0.05, \|log2FC\|>1） | 7283（上调 4355 / 下调 2928） |
| EGFR 表达 | log2FC ≈ -0.21；padj ≈ 0.56（**不显著**） |
| EGFR 突变 | 2/15 例激活突变：E746_A750del、L858R |
| 可成药面板 | 247 个基因 |
| 表达驱动候选（shortlist） | 78 个，top5：TERT、HTR3A、TOP2A、CDKN2A、SLC6A3（EGFR 不在其中） |
| 突变驱动候选 | 9 个基因：TP53 8、KRAS 6、STK11 3、RET 2、KEAP1 2、EGFR 2、SMARCA4/PIK3CA/ROS1 各 1 例 |
| GO BP 显著条目 | 1054（top：extracellular structure organization，159 基因；extracellular matrix organization 158） |
| Reactome 显著条目 | 82（top：Extracellular matrix organization，138 基因） |
| 数据缺口 | TCGA-05-4245-01A 无 MAF（14/15 例有突变数据）；MAF 不含扩增/融合 |

> 上面每个值都由 `docs/课程包自检结果.md` 记录的运行环境实测得到；改动数据或脚本后必须重算。

## 数据替换指南（换主题）

本数据是直接用 GDC 官方文件重建的，不依赖手工编造，换主题只需改选择规则：

1. 修改 `fetch_tcga_luad.py` 里的项目 ID（例如 `TCGA-HNSC`）与样本筛选条件，或直接换掉 `sample_manifest.csv` 里的 file_id；
2. 调整 `druggable_panel.csv`（可按 Open Targets / DGIdb 重新整理面板类别）；
3. 突变基因列表（脚本里的驱动基因白名单）按新主题调整；
4. 重新生成数据后**务必重跑三份 R 脚本**，并按新的实测值更新所有文档中的基准数字（清单见上一节）。

> 换主题会影响三课衔接：课程二/三的结构与设计案例围绕 EGFR，换主题时需要同时替换晶体结构（1M17/1YY9）与抗体案例。
