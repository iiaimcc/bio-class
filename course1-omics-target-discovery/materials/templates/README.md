# 变异功能标注（学生任务）

## 起点：先让脚本把能自动判的判掉

```bash
# 在 course1-omics-target-discovery/ 目录下
python3 ../data/make_variant_template.py ../data ./results/variant_calls.csv
```

生成的表 **27 行 = 每条非沉默变异一行**（不是每个基因一行；KRAS 的 p.G12C 与 p.G12V 各有两位患者，所以 25 个唯一变异对应 27 行）。列结构：

| 列 | 谁填 | 说明 |
|---|---|---|
| `symbol` | 脚本 | 基因名 |
| `sample_id` | 脚本 | 患者样本号（同一条变异出现在不同患者时，靠这一列区分） |
| `hgvsp` | 脚本 | 蛋白水平变异写法（如 `p.L858R`） |
| `variant_class` | 脚本 | MAF 的变异类型（如 `Missense_Mutation`、`Nonsense_Mutation`） |
| `functional_consequence` | **脚本填 4 行，其余你填** | `gain_of_function` / `loss_of_function` / `uncertain` |
| `evidence_level` | 你填 | `high` / `moderate` / `none found` |
| `evidence_source` | 你填 | 能被人打开核对的来源：ClinVar VCV 号、CIViC id、PMID（查不到就写清"checked <日期>，无记录"） |

截短变异（无义 `Nonsense_Mutation`、移码 `Frame_Shift_*`）由脚本自动判为 `loss_of_function`：蛋白被打坏了，不可能"打开"。**剪接位点、起始密码子、非终止变异不自动判**（本队列没有这类变异）——它们的效应取决于转录本、阅读框或根本不是截短，遇到时必须自己查。

剩下的 23 行（错义变异与框内缺失）是你填。`variant_calls_template.csv` 与生成结果同构，留作离线对照；`variant_calls_reference.csv` 是教师参考答案（含来源）。

## 为什么这一步不自动做

**MAF 只告诉你"哪个氨基酸变了"，不告诉你"变了以后蛋白是开着还是关着"。** 这个判断要查外部证据，而且逐条判：

- 同一个基因上的不同变异可以有不同结论。TP53 有 1 条无义变异（自动判功能缺失）+ 4 条错义变异（查库判功能受损）+ 4 条错义变异（查无记录，标未明）——不能因为基因里有一条截短变异，就把整个基因的所有变异都算成功能缺失；
- 同一条变异出现在不同患者时，两行都要填，脚本用 `symbol + sample_id + hgvsp + variant_class` 组合键核对完整性。

## 怎么做

1. 对每一条 `needs_review`，用现成脚本查两个数据库（只查、不下结论）：

   ```bash
   python3 ../data/lookup_variant.py EGFR p.L858R
   python3 ../data/lookup_variant.py PIK3CA p.P104L
   ```

   脚本把短写法（`p.L858R`）换算成 ClinVar 检索用的 HGVS（`p.Leu858Arg`），打印 ClinVar 记录（分级、审核状态、关联疾病）与 CIViC 已接受证据（证据类型 / 方向 / 等级）。

2. 需要网页原文时：ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>，CIViC <https://civicdb.org/>，OncoKB <https://www.oncokb.org/>（部分内容需注册）。

3. 把结论写进 `functional_consequence` / `evidence_level` / `evidence_source` 三列：

   - `evidence_level`：多个独立来源一致＝`high`；只有间接或单一语境＝`moderate`；查不到记录＝`none found`；
   - `evidence_source`：必须能被人打开核对（CIViC id、ClinVar VCV、PMID 或 URL 任选，多个来源用 `;` 分隔）。

   参考行的写法（截自 `variant_calls_reference.csv` 的 `evidence_source` 列）：

   ```text
   CIViC 78 (SENSITIVITYRESPONSE level A); ClinVar VCV001701193 / VCV000012578
   ```

4. 把表存成 `results/variant_calls.csv`，运行 `03_target_shortlist.R`（第三个参数省略时读这个文件）。脚本会逐项检查并在报错里点名具体行：

   | 检查 | 报错示例 |
   |---|---|
   | 必填列齐全 | `variant_calls.csv must have the columns: …` |
   | `functional_consequence` 只能是三个允许值（不能留 `needs_review`） | `Unresolved calls in …: TP53 p.Q317*` |
   | `evidence_level` 只能是 `high` / `moderate` / `none found` | `evidence_level must be one of: high, moderate, none found - fix these rows: …` |
   | `evidence_source` 不能留空或占位符（`NA`、`-`、`tbd` 等） | `evidence_source must be filled for every row …` |
   | 每条变异记录齐全（组合键 `symbol+sample_id+hgvsp+variant_class`） | `variant_calls.csv has no row for: KRAS / TCGA-55-8615-01A / p.G12C` |
   | 没有重复行 | `variant_calls.csv has duplicate rows for the same record: …` |
   | 没有不属于本队列的行 | `variant_calls.csv has rows that are not in mutation_status.csv: …` |

   通过后按基因汇总（`n_gof` / `n_lof` / `n_uncertain`；两类都有时标 `mixed`）。

## 三条判读纪律

- **查不到 ≠ 阴性。** 数据库没有记录，只能写 `uncertain`，不能写成"野生型"或"无关"。
- **注意证据来自哪个语境。** 一个变异可能在别的疾病（例如生殖系过度生长综合征）里被判定为致病，而在肿瘤里没有体细胞证据——照实写清楚，别丢掉线索，也别把它当成肿瘤驱动证据。
- **不要把"有没有药"混进 `functional_consequence`。** 这里只回答变异对蛋白功能做了什么；药物与文献属于事后佐证。
