# 1YY9 抗体链与 CDR 标注表

> 1YY9 的真实链为 A=EGFR、C=Fab light、D=Fab heavy。本表由 1YY9 SEQRES 提取 V 区序列、用 ANARCI（IMGT 方案）编号后换算到 PDB `resi`，可复算。
> 不要把 IMGT 编号直接当成 PDB `resi` 编号：IMGT 编号含缺口与插入位，PDB 编号是从 1 开始的连续编号。

## 链映射

| PDB 链 | 分子 | SEQRES 长度 | 建模范围 | 设计/分析用途 |
|---|---|---:|---|---|
| A | EGFR 胞外域 | 624 | 2–614 | 抗原与表位 |
| C | cetuximab Fab light chain（VL 1–107 + CL） | 213 | 1–211 | 轻链互补位 |
| D | cetuximab Fab heavy chain（VH 1–119 + CH1） | 221 | 1–220 | 重链互补位 |

## CDR 区间（IMGT 方案）

| 区域 | IMGT 区间 | PDB 链 | PDB `resi` 区间 | 序列（PDB 编号顺序） |
|---|---|---|---|---|
| CDR-L1 | 27–38 | C | 27–32 | `QSIGTN` |
| CDR-L2 | 56–65 | C | 50–52 | `YAS` |
| CDR-L3 | 105–117 | C | 89–97 | `QQNNNWPTT` |
| CDR-H1 | 27–38 | D | 26–33 | `GFSLTNYG` |
| CDR-H2 | 56–65 | D | 51–57 | `IWSGGNT` |
| CDR-H3 | 105–117 | D | 96–108 | `ARALTYYDYEFAY` |

说明：

- IMGT 区间是标准定义（CDR1 = 27–38、CDR2 = 56–65、CDR3 = 105–117）；本抗体的 V 区较短，区间内多处 IMGT 位点为空缺，所以换算到 PDB 后是较短的连续段，属正常现象；
- CDR-H3 的 IMGT 定义包含 `CAR` 中的 `AR`，与 Kabat 的 `ALTYYDYEFAY` 相差两个残基；本课程统一用 IMGT，不混用 Kabat/Chothia；
- 三条重链 CDR 与三条轻链 CDR 都参与 4.5 Å 界面（见 `fallback/expected/1YY9_界面参考.md`），这是"cetuximab 用 CDR 识别 domain III"的直接结构证据。

## 复算方法

```bash
# 1. 从 SEQRES 还原链序列（或 RCSB 页面下载 FASTA）
# 2. V 区截取：VH = D 链 1–119，VL = C 链 1–107（CH1/CL 起始处截断）
ANARCI -i vh.fa --scheme imgt --restrict heavy
ANARCI -i vl.fa --scheme imgt --restrict light
# 3. 按"非缺口位置顺序"把 IMGT 编号映射回 PDB resi
```

本次标注使用 ANARCI（bioconda，2026-11-02）+ 1YY9 SEQRES；换成 abYsis/IMGT 网页时，边界应一致，若不同先核对链序列版本。

## 教师核验步骤

1. 用上表在 PyMOL 中逐段 `select`，确认每段非空：
   ```pymol
   select cdrL1, chain C and resi 27-32
   select cdrL2, chain C and resi 50-52
   select cdrL3, chain C and resi 89-97
   select cdrH1, chain D and resi 26-33
   select cdrH2, chain D and resi 51-57
   select cdrH3, chain D and resi 96-108
   count_atoms cdrL1 or cdrL2 or cdrL3 or cdrH1 or cdrH2 or cdrH3
   ```
2. 展示 CDR 与 EGFR 接触残基的空间关系（`pymol_interface_commands.pml`）；
3. 说明 CDR 是序列/编号定义，4.5 Å 接触区是结构几何定义，两者不能互相替代。
