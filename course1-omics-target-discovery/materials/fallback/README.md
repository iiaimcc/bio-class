# 课程一兜底材料（fallback）

课上的数据由学生**现场从 GDC 下载**（`data/fetch_tcga_luad.py`，任务卡第 0 步）。本目录只接住两件事：机房网络不可用、或学生想核对教师基准。

## 目录

```
fallback/
├── README.md                 # 本文件：来源、校验值、使用规则
├── inputs/
│   └── tcga_luad_subset.tar.gz   # 与 GDC 同源的 8 个数据文件（含 manifest）+ target_context.csv + drug_precedent.csv
└── expected/                 # 教师运行三个脚本后的参考输出
    ├── deg_symbols.txt               # 7283 个 DEG 的基因名（供 Enrichr 上传）
    ├── go_bp_enrichment.csv          # GO BP 富集结果（1054 条）
    ├── reactome_enrichment.csv       # Reactome 富集结果（82 条）
    ├── target_shortlist.csv          # 表达驱动候选（78 个）
    ├── mutation_driven_candidates.csv# 突变驱动候选（9 个基因）
    ├── integrated_candidates.csv     # 三层漏斗逐层结果（9 行）
    ├── target_recommendation.csv     # 三层全过的综合推荐（EGFR，rank 1）
    ├── volcano_plot.png
    ├── go_bp_dotplot.png
    └── reactome_dotplot.png
```

`deg_table_all.csv`（23765 行）不放进仓库（3.7 MB），教师用 `results_teacher/` 自行分发。

## 输入文件来源与校验

抓取日期：2026-11-02，来自 NCI GDC 公开层（TCGA-LUAD），下载后未做编辑，只做合并与压缩。

| 文件 | SHA256（解压后） |
|---|---|
| `expression_matrix.csv` | `1da1bdca79bba6787bc69b9dfb1eafc3431b4a95928697871069342c9109c576` |
| `sample_meta.csv` | `e9826e41c85b9cd7d68613c74c96452e4e18548690f0f2532d92991d155cb1fd` |
| `gene_annotation.csv` | `32c0c37fad1edffe8c0cf138a18ba6937dd52f7dad23e260021d3ef5ed987713` |
| `druggable_panel.csv` | `b8e19d1bc2ba67f7b11e93eb94b02b5b337166e31a540feaf2eb9deb23494d6b` |
| `mutation_status.csv` | `f7caf3facb1c150dff8f48bc7853f1f8bbf785abdd7d5f9d06d63dc90b8d31a2` |
| `mutation_matrix.csv` | `349412bd4b4f59e198ca8f7f4b352880b56862214314bf74fbacbe07527471b8` |
| `sample_manifest.csv` | `98234d369cca5812523f21711fa135860d49409c06b3423d09379fbded30d1b4` |
| `mutation_manifest.csv` | `cde2eeb3d544495c9d775d307d868457c3b392938cda7c04edf3ff610f38eefa` |
| `target_context.csv` | `537960a9c51a73433265a56c7ea7525eea994dca578053c322e8efaf6db6c816` |
| `drug_precedent.csv` | `d27e23d0af48140e4582b5ddebe2cecb89fe3b8bd7fcaf38945c9eac3719254f` |

压缩包本身：`846a4dfa9b7dc98b08e9ffcc0cb404db6bcc16d62afdc7d26031db24c740642`（2.2 MB，含 10 个文件：8 个数据文件 + `target_context.csv` + `drug_precedent.csv`）。

与 `data/` 下的发布版逐字节一致；`python3 data/fetch_tcga_luad.py --verify-only` 可现场复核。

## 使用规则

1. 先让学生走第 0 步的下载脚本。只有在机房网络不可用、或课堂时间不够时才解压本包：
   ```bash
   tar -xzf materials/fallback/inputs/tcga_luad_subset.tar.gz -C ./data
   ```
2. 变异功能表不在兜底包里：`variant_calls.csv` 是学生任务，断网时用 `materials/templates/variant_calls_reference.csv`（教师参考答案，含来源）核对，或按模板自查后填入。
3. `expected/` 只用于两处：教师快速批改；学生断网时对照判读（例如把 `deg_symbols.txt` 拖进网页版 Enrichr）。教师复现基准的完整命令：

   ```bash
   Rscript scripts/01_diff_analysis.R ../data ./results_teacher
   Rscript scripts/02_enrichment.R ./results_teacher ./results_teacher
   Rscript scripts/03_target_shortlist.R ../data ./results_teacher ../materials/templates/variant_calls_reference.csv
   ```
4. `expected/` 里的数字与 `data/README.md` 的教学基准一致；一旦数据或脚本改动，必须重跑并同步更新，不能只改一处。
