# boltz_out（待教师生成）

本目录是课程三可选的 Boltz 参照复核结果的存放位置，当前为空占位。仅当教师实际运行
`boltz_antibody_template.yaml` 后，把输出与指标放入本目录，并填写
`../../boltz_预生成_EGFR_cetuximab.md`。

## 需要放入的文件

| 文件 | 内容 |
|---|---|
| `boltz_egfr_cetuximab.yaml` | 实际运行的输入（三条真实序列 + 链映射） |
| `prediction/`（或服务默认输出目录） | 预测结构 PDB/CIF、置信度 JSON |
| `confidence_summary.md` | iptm/ptm、pAE 图路径、置信区域说明 |
| `align_vs_1YY9_global.png` | 与 1YY9 叠合的全局图 |
| `align_vs_1YY9_interface.png` | 界面放大图 |

## 规范

- 指标必须来自运行输出，禁止手填或估算；
- 记录运行日期、GPU、`--recycling_steps`、`--diffusion_samples`、MSA 模式；
- 叠合只比较 EGFR 胞外域公共段，写清用了哪些残基；
- ipTM/pAE/RMSD 不得写成 KD 或阻断能力。
