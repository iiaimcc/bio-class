# 课程三在线设计服务与课前实测

> 在线服务是课堂替代路径，不代表每次都可用。教师必须在课前一天实测登录、上传、运行和下载；学生报告必须写实际使用的工具。

## 1. RFdiffusion / backbone 生成

- ColabDesign RFdiffusion notebook：
  https://colab.research.google.com/github/sokrypton/ColabDesign/blob/main/rf/examples/diffusion.ipynb
- 公开代码仓库：
  https://github.com/RosettaCommons/RFdiffusion
- 可选无代码网页服务（教师先确认账号、额度与隐私条款）：
  https://app.tamarind.bio/tools/rfdiffusion

课堂固定：目标 PDB、目标表位、binder 长度、生成数量、随机种子。不要让学生在课堂中自行改变所有参数。

## 2. ProteinMPNN / backbone 到序列

- 官方 ProteinMPNN quick demo：
  https://colab.research.google.com/github/dauparas/ProteinMPNN/blob/main/colab_notebooks/quickdemo.ipynb
- 官方带 AF2 回折的 demo：
  https://colab.research.google.com/github/dauparas/ProteinMPNN/blob/main/colab_notebooks/quickdemo_wAF2.ipynb

学生必须记录设计链、固定链、序列数量和输出文件。ProteinMPNN 分数不能写成 KD。

## 3. 回折与复合物检查

- 教师可选 ColabFold：课前确认多链输入、GPU 配额和结果下载；
- 教师可选在线 Boltz-2：
  https://app.tamarind.bio/boltz
- 教师可选 Rowan Boltz-2：
  https://www.rowansci.com/features/free-online-boltz-2

网页服务的结果、输入保留政策和账号要求可能变化。课堂中以教师课前实测可用的一个服务为主，其余只作备选。

## 4. 课堂失败兜底

| 现象 | 课堂处理 | 学生应提交什么 |
|---|---|---|
| Colab 无 GPU/排队 | 使用教师预生成 backbone（BoltzGen 替代流程） | 参数表、预生成 CIF、结构审查 |
| 目标 PDB 上传失败 | 使用教师准备的目标局部结构 | 输入文件名和失败截图 |
| ProteinMPNN 不可用 | 教师提供脱敏候选序列表 | 候选比较和回折计划 |
| 回折服务不可用 | 使用教师预生成回折结果 | 判读表，不伪造运行成功 |

禁止上传患者级、未公开或受限数据；本课使用公开 PDB/UniProt 序列和教师整理的教学结构。
