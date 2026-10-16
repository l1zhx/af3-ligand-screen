# AI4S portfolio and publishing guide

## Project positioning

Suggested repository name: **af3-ligand-screen**.

Suggested GitHub About text:

> Reproducible AlphaFold 3 protein–ligand screening with Slurm, MSA reuse, confidence and geometry auditing, and docking comparison.

Suggested topics: `alphafold3`, `ai-for-science`, `computational-biology`, `protein-ligand`, `slurm`, `rdkit`, `virtual-screening`, `reproducible-research`.

For AI4S applications, this project demonstrates turning a pretrained scientific model into a traceable research workflow: molecular data engineering, GPU orchestration, result quality control, evaluation design and scientific interpretation. It does not establish model architecture development or AF3 training experience.

## Resume wording

中文，可按你实际负责的工作调整：

> 面向 β-catenin 蛋白–小分子筛选，组织约 6 万分子规模的 AlphaFold 3 集群计算，整合 RDKit InChIKey 去重、单蛋白 MSA 复用、Slurm 作业数组与失败追踪；构建 sample 级结果审计、界面置信度排序和结构几何检查流程，并设计与 Glide 聚类结果的化学实体匹配及候选交叉复核。

English:

> Developed and operated an AlphaFold 3 screening workflow for a roughly 60,000-ligand β-catenin campaign, integrating RDKit identity tracking, reusable receptor MSA inputs, Slurm job arrays, and sample-level output auditing. Implemented interface-confidence ranking and structural triage, with standardized-parent matching to existing Glide clusters for candidate review.

Use “implemented” only for components you understand and can explain, and adjust the statement to your actual contribution. State that coding assistants were used when appropriate in an interview or contribution description. Do not claim measured acceleration, a validated hit rate, discovered inhibitors, model training or a completed full-library comparative analysis without corresponding evidence.

## Interview walkthrough

1. Explain why the receptor MSA is fixed and reusable, while AF3 inference still costs GPU time.
2. Show how a ligand moves from source SDF record to InChIKey, manifest ID, unique input, Slurm attempt and selected sample.
3. Demonstrate the failure case where bare `summary_confidences.json` filenames caused zero recognized samples, and the regression test that now covers it.
4. Explain why submission ledgers and runtime locks address different parts of accidental duplication, and how partial failures are recovered.
5. Compare A–B ipTM, ranking_score, ligand pLDDT, PAE and contact geometry without equating them with affinity.
6. Explain why importing docking Top 10,000 creates a selected comparison cohort and why standardized parents can merge molecular states.
7. Run the CPU demo and point to exactly which tests and scientific claims it supports.

## Publish the repository

The supplied project directory is a standalone repository bundle. Production data, local settings, results, caches and model weights are excluded. From **inside this project directory**, create a Git repository and review staged filenames:

```bash
git init -b main
git add README.md README.zh-CN.md LICENSE requirements.txt requirements-chem.txt .gitignore .gitattributes .github config scripts slurm tests docs examples
git status --short
git diff --cached --stat
git commit -m "Add AF3 ligand screening workflow and reproducible demo"
```

Create an empty repository in your GitHub account. Replace `YOUR_USERNAME` below with your actual username, then push:

```bash
git remote add origin https://github.com/YOUR_USERNAME/af3-ligand-screen.git
git push -u origin main
```

Before publication, confirm your right to share the project construct and any added research results. The bundle contains synthetic structures only, with the project FASTA as the biological example. The code license does not grant rights to AF3 parameters or private compound libraries.

## Add it to your GitHub profile

Use GitHub's profile customization to pin the project. If you maintain a profile README repository named exactly as your GitHub username, add the text in [profile_snippet.md](profile_snippet.md) after replacing its username placeholder. The project README is designed to stand on its own; a separate website is not required.

The upload commands are provided for you to execute. This bundle has not been published or pushed to a GitHub account automatically.
