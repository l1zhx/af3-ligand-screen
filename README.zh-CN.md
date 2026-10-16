# 基于 AF3 的蛋白–小分子筛选工作流

本项目把 β-catenin 约 6 万分子的实际计算经历整理成可复用的科研工程项目，覆盖 **数据准备 → 单蛋白 MSA → Slurm 推理 → 结果审计 → 结构检查 → 与 docking 对照 → 候选导出**。

原始项目配体清单包含 59,635 个分子，已报告完成服务器推理。当前确认的真实后处理证据是一个配体的 5 个 sample 全部被正确读取；全库分析成功率、最终候选数量、实验活性还没有在公开仓库中得到确认。示例目录的数据由程序人工生成，不作为科学发现。

## 先运行演示

```bash
python -m pip install -r requirements.txt
python scripts/demo.py --out demo_output
python -m unittest discover -s tests -v
```

不用 GPU、Slurm 或 AF3 权重。演示包含完整结果、重复 attempt、部分结果、缺失结果和有碰撞标记的结果。预期排名 2 个配体、保留 10 个 sample、排除 1 个重复 attempt、导出 1 个结构。再次运行时换一个输出目录。

## 在集群执行

按 [HPC 完整操作指南](docs/HPC.md) 执行；下面是各阶段入口。

| 阶段 | 文件 | 资源 |
| --- | --- | --- |
| SDF 清洗和去重 | `scripts/prepare.py library` | CPU + RDKit |
| 单蛋白 MSA 输入 | `scripts/prepare.py msa-input` | CPU |
| MSA/templates 搜索 | `slurm/msa.slurm` | CPU + AF3 数据库 |
| 复合物 JSON | `scripts/build_jobs.py` | CPU |
| 作业计划与提交 | `scripts/workflow.py plan/submit` | 登录节点轻量操作 |
| 推理 | `slurm/infer.slurm` | 每个数组任务 1 张 GPU |
| 置信度审计、结构检查 | `slurm/postprocess.slurm` | CPU + NumPy |
| docking 对照 | `scripts/merge_af3_glide.py` | CPU + RDKit |
| 公开汇总 | `scripts/make_report.py` | CPU；绘图可选 |

`workflow.py submit` 是提交辅助程序，实际调用的是 **`sbatch ... slurm/infer.slurm`**。文件扩展名不决定调度方式；只有经 `sbatch` 提交后才是 Slurm 作业。

`config/cluster.local.sh` 保存你自己的路径，不进入 Git。在旧项目上做后处理，只需将 `WORK_DIR` 指向已有 `jobs/`、`results/`、`ligands.csv` 所在目录。不要为了使用新仓库重跑已经完成的 AF3 推理。

## 看哪些结果

先检查 `summary/audit.json` 和 `coverage.csv`，确认实际发现的 sample 数量、完整配体数和错误，再看 `ligand_metrics.csv`。根目录最优结构是 sample 的副本，不是第 6 个 sample。文件名同时支持 `summary_confidences.json` 和带任务名前缀的形式。

候选按已知无碰撞优先、A–B ipTM 从高到低、再按 `ranking_score` 排序。`geometry/` 进一步检查配体 pLDDT、双向 PAE、4 Å 接触、极短接触和蛋白对齐后的姿势一致性，导出待人工复核的结构。**排名、姿势一致性都不能直接证明结合或抑制活性。**

`merge_af3_glide.py` 接收既有 docking 聚类表，通过标准化母体结构匹配，而不依赖两个流程编号恰好相同。若 docking 表只含 Top 1 万，未匹配的 AF3 配体不应直接判为阴性。详见 [方法与限制](docs/METHODS.md)。

## AI4S

展示科研问题、数据追踪、GPU 作业调度、失败诊断、结构评价和可复现性。英文首页适合快速查看，中文指南用于实际复现。[作品集说明](docs/PORTFOLIO.md) 提供GitHub 上传步骤和 Profile README 片段。所有吞吐量、成功率和命中率应来自实际审计；不要把脚本设置的并行上限写成实测性能。
