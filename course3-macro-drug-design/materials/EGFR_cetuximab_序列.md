# EGFR 胞外域 + 西妥昔单抗（Cetuximab）序列与结构资料

> 本文件中的序列全部由公开数据库文件提取，可复算：EGFR 来自 UniProt P00533，抗体来自 RCSB 1YY9 SEQRES。
> 学生课上自行下载同一批文件后核对，不应把本文件当成"无需下载"的替代品。

## 1. 事实基线

- **复合物结构**：PDB `1YY9`（西妥昔单抗 Fab–EGFR 胞外域复合物，X-ray，**2.61 Å**）。
  真实蛋白链为：链 A = EGFR 胞外域，链 C = Fab 轻链，链 D = Fab 重链。**链 B/E/F/G 是 N-糖链**（NAG/NAG-BMA-MAN 分支），不是蛋白链；A、D 链的蛋白残基上还连着 NAG。
- **EGFR（UniProt P00533）**：全长 1210 aa；信号肽 1–24；胞外域 25–645（UniProt 注释），
  由 domain I–IV 组成；配体 EGF 结合位点位于 domain I 与 III 之间的界面；
- **1YY9 的覆盖范围**：PDB A 链 1–618 对应 UniProt **25–642**（DBREF），建模到 2–614；即该晶体结构没有覆盖 643–645 与 C 端。
- **domain III**：UniProt **334–504**；在 1YY9 中是 PDB **310–480**（偏移 +24）。cetuximab 结合在 domain III。
- **西妥昔单抗**：嵌合（chimeric）单抗（鼠源 V 区 + 人源 C 区），靶向 EGFR domain III，
  通过空间位阻阻断 EGF 结合；FDA 获批用于结直肠癌/头颈癌；
- **Fab vs 全长 IgG**：晶体中为 Fab（VH-CH1 / VL-CL），本课用 VH+VL 两条链做参照。

## 2. 链与边界（可复算）

| 链 | 身份 | SEQRES | 建模范围 | V 区边界 |
|---|---|---:|---|---|
| A | EGFR 胞外域 | 624 | 2–614 | — |
| C | Fab 轻链 | 213 | 1–211 | VL = 1–107；CL 从 108 开始（`RTVAAPSVFIFPPSDEQLK`） |
| D | Fab 重链 | 221 | 1–220 | VH = 1–119；CH1 从 120 开始（`ASTKGPSVFPLAP`） |

复算：

```bash
curl -L https://files.rcsb.org/download/1YY9.pdb -o 1YY9.pdb
grep "^DBREF" 1YY9.pdb
grep "^SEQRES" 1YY9.pdb | awk '{print $3, $4, NF}'
```

## 3. CDR 标注（IMGT）

完整区间、换算到 PDB `resi` 的结果和复算方法见 [CDR标注表.md](CDR标注表.md)。要点：

- 重链 VH 在 **链 D**（不是链 B）；
- CDR-H1/H2/H3 的 PDB 区间为 D 26–33、51–57、96–108；
- 轻链 CDR-L1/L2/L3 为 C 27–32、50–52、89–97。

## 4. 序列（可直接用于 Boltz/ColabFold 输入）

```text
>EGFR_DIII UniProt P00533 334-504 (171 aa; 1YY9 PDB 310-480)
RKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQ

>Cetuximab_VH 1YY9 chain D 1-119 (119 aa)
QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSA

>Cetuximab_VL 1YY9 chain C 1-107 (107 aa)
DILLTQSPVILSVSPGERVSFSCRASQSIGTNIHWYQQRTNGSPRLLIKYASESISGIPSRFSGSGSGTDFTLSINSVESEDIADYYCQQNNNWPTTFGAGTKLELK
```

> 抗体序列来自 1YY9 SEQRES（晶体中 C 链缺 212–213、D 链缺 221，均在恒定区，不影响 VH/VL）。
> 保留"残基编号→序列位置"对照：VH 第 n 位 = 链 D `resi n`；VL 第 n 位 = 链 C `resi n`。
>
> **EGFR domain III 的序列与晶体坐标有一处单点差异**：UniProt canonical 498 位是 S，1YY9 SEQRES（PDB 474）是 K。
> 序列输入用 UniProt canonical（上表）；与 1YY9 做结构比对时以 PDB 坐标为准，不要按序列去改晶体文件。

## 5. 常见错误

| 错误写法 | 正确事实 |
|---|---|
| 1YY9 分辨率 3.4 Å | 2.61 Å（`REMARK   2 RESOLUTION`） |
| VH 从链 B 提取 | VH 在链 D；链 B 是 N-糖链 |
| 胞外域 1–621 或 25–645 混用 | UniProt 注释 25–645；1YY9 实际覆盖 25–642（PDB 1–618） |
| 1M17 的 T790 = `resi 790` | 1M17 中为 `resi 766`（偏移 +24） |

## 6. 参考文献

- Li S. et al. *Structural basis for inhibition of the epidermal growth factor receptor by cetuximab*. Cancer Cell 2005;7(4):301-11.（1YY9 出处）
- UniProt P00533 注释（EGFR）
- RCSB PDB 1YY9 页面
