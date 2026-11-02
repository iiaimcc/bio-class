# 课程二兜底材料（fallback）

课上的序列与结构数据由学生**现场从数据库获取**（UniProt、RCSB PDB、AlphaFold DB），这是任务卡 1–2 的操作内容。本目录只接住网络/服务失败，不替代操作。

## 目录

```
fallback/
├── README.md                 # 本文件：来源、校验值、使用规则
├── inputs/                   # 与线上同源的原始文件（教师课前抓取）
│   ├── EGFR_P00533.fasta
│   ├── 1M17.pdb
│   └── AF-P00533-F1-model_v6.pdb
└── expected/
    └── 结构事实参考.md        # 各步骤应得到的参考输出
```

## 输入文件来源与校验

抓取日期：2026-11-02（下载后未做任何编辑）。

| 文件 | 来源 URL | 大小 | SHA256 |
|---|---|---:|---|
| `EGFR_P00533.fasta` | https://rest.uniprot.org/uniprotkb/P00533.fasta | 1328 B | `684e2da1f3ef05d0ce90d0aefe9125c6f7e9cdcfeac27ae72ed35243b64a73e8` |
| `1M17.pdb` | https://files.rcsb.org/download/1M17.pdb | 252720 B | `5b020cd186307d657de80611793554b074f6936a3443d33f81d9b1a58004962e` |
| `AF-P00533-F1-model_v6.pdb` | https://alphafold.ebi.ac.uk/files/AF-P00533-F1-model_v6.pdb | 772334 B | `8b381730ae157a75f2fbc78043aff7b9e605a30420dfe749a56cd61f12125cb9` |

> AlphaFold DB 的 v4 文件已不再提供下载，课程统一使用 **v6**。模型内容若更新，以 `https://alphafold.ebi.ac.uk/api/prediction/P00533` 返回的 `latestVersion` 为准；更新后必须同步本表与 `expected/结构事实参考.md`。

## 使用规则

1. 先按任务卡现场下载，并在产出中记录 URL、访问日期和校验值；
2. 只有在下载失败、超时或文件损坏时才使用本目录，并在报告中写明：**使用了兜底文件、失败原因、与现场下载的差异**；
3. 不允许把兜底文件写成"我现场下载"；
4. 教师每学期开课前重新抓取一次并更新上表，防止数据库版本漂移。
