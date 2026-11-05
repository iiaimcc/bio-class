# 变异功能标注（学生任务）

`variant_calls_reference.csv` 是教师填好的参考答案（含来源）。学生不直接用模板，而是先运行生成脚本，让脚本把能自动判的填好：

```bash
# 在 course1-omics-target-discovery/ 目录下
python3 ../data/make_variant_template.py ../data ./results/variant_calls.csv
```

生成的表 27 行（每条非沉默变异一行）：4 条截短变异已自动判为 `loss_of_function`，其余 23 行是 `needs_review`。`variant_calls_template.csv` 与生成结果同构，留作离线时的对照。

## 为什么这一步不能自动做

漏斗第一层（有没有非沉默突变）与"截短变异 = 蛋白被打坏"由数据本身就能判定，脚本已经自动完成（TP53、STK11、KEAP1 因此直接归为功能缺失）。剩下 6 个基因携带的是错义变异与框内缺失——**MAF 只告诉你"哪里的氨基酸变了"，不告诉你"变了以后蛋白是开着还是关着"**。这个判断需要查外部证据，所以留给你做。

## 怎么做

1. 对每一条 `needs_review`，用现成脚本查两个数据库（只查、不下结论）：

   ```bash
   # 在课程目录（course1-omics-target-discovery/）下运行：
   python3 ../data/lookup_variant.py EGFR p.L858R
   python3 ../data/lookup_variant.py PIK3CA p.P104L
   ```

   脚本把短写法（MAF 里的 `p.L858R`）换算成 ClinVar 检索用的 HGVS（`p.Leu858Arg`），打印 ClinVar 记录（分级、审核状态、关联疾病）与 CIViC 已接受证据（证据类型 / 方向 / 等级）。

2. 需要看网页原文时：ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>，CIViC <https://civicdb.org/>，OncoKB <https://www.oncokb.org/>（部分内容需注册）。

3. 把结论写进表格，每行一个变异：

   | 列 | 允许的取值 | 说明 |
   |---|---|---|
   | `functional_call` | `gain_of_function` / `loss_of_function` / `uncertain` | 这个变异把蛋白打开、关掉，还是证据不足 |
   | `evidence_level` | `high` / `moderate` / `none found` | 多个独立来源一致＝high；只有间接或单一语境＝moderate；查不到记录＝none found |
   | `source` | URL 或数据库编号（CIViC id、ClinVar VCV、PMID） | 必须能被人打开核对 |
   | `note` | 一句话 | 证据来自哪个语境、有什么限制 |

4. 把填好的表存成 `results/variant_calls.csv`（`results` 是你跑脚本的输出目录），再运行 `03_target_shortlist.R`。

## 三条判读纪律

- **查不到 ≠ 阴性。** 数据库没有记录，只能写 `uncertain`，不能写成"野生型"或"无关"。
- **注意证据来自哪个语境。** 一个变异可能在别的疾病（例如生殖系过度生长综合征）里被判定为致病，而在肿瘤里没有体细胞证据——照实写清楚，别丢掉线索，也别把它当成肿瘤驱动证据。
- **不要把"有没有药"混进 `functional_call`。** 这里只回答变异对蛋白功能做了什么；药物与文献属于事后佐证。
