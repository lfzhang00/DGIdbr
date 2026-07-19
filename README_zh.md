<p align="right">
  <a href="README.md">English</a> | <b>中文</b>
</p>

# DGIdbr

[![R CMD check](https://github.com/lancelotzhang0124/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lancelotzhang0124/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.2.0-blue)](https://github.com/lancelotzhang0124/DGIdbr)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

**基于统计检验的药物优先排序工具。**

给定一组差异表达基因，DGIdbr 查询 [DGIdb](https://dgidb.org/) 药物-基因互作数据库，
并通过**超几何富集检验**（而非简单的命中计数）对药物进行排序。它回答的核心问题是：
*"这个药物与我的基因集的重叠是否显著高于随机期望？"*

## 为什么重要

| 方法 | 做法 | 缺陷 |
|------|------|------|
| **简单计数** | 按药物命中了多少个输入基因排序 | 泛靶点药物（如广谱化疗药）永远排第一 |
| **DGIdbr (v1.2+)** | 超几何检验 + FDR + 富集倍数 | 特异性强、机制合理的药物脱颖而出 |

**示例**：输入 9 个癌症相关基因，旧版计数法将 RIBAVIRIN 和 CISPLATIN 排在最前面
——它们都靶向数百个基因。富集法将 CHEMBL1214407 排在第一：该化合物仅有 4 个已知靶点，
但其中 2 个在输入基因集中——**648 倍富集**（FDR = 0.001）。

## 功能特性

- **超几何富集检验**，带 Benjamini–Hochberg FDR 校正
- **动态背景校准**——自动查询 DGIdb 当前可成药基因总数（~11,665），而非硬编码猜测值
- **富集倍数**（enrichment ratio）作为直观的效应量指标
- **分组模式**：从 case-control 差异表达表中提取 up/down 基因集
- **亚型模式**：自动识别所有亚型，为每个亚型构建 up/down 基因集
- **FDA 批准过滤**——仅保留已批准药物，或纳入研究性药物
- **完全向后兼容**——设置 `enrichment = FALSE` 即可回退到旧的计数模式

## 安装

```r
# install.packages("remotes")
remotes::install_github("lancelotzhang0124/DGIdbr")
```

## 快速开始

### 分组模式（case-control）

输入 CSV 需包含 `gene` 和 `direction`（`up` / `down`）列：

```r
library(DGIdbr)

DGIdbr(
  mode           = "group",
  base_tables    = "path/to/input",
  group_filename = "group.csv",
  base_out       = "path/to/output",
  approve        = TRUE,       # 仅保留 FDA 已批准药物
  enrichment     = TRUE         # 超几何检验 + FDR（默认开启）
)
```

### 亚型模式

输入 CSV 需包含 `gene`、`direction` 和 `subtype` 列：

```r
DGIdbr(
  mode             = "subtype",
  base_tables      = "path/to/input",
  subtype_filename = "subtype.csv",
  base_out         = "path/to/output",
  approve          = TRUE
)
```

### 向后兼容（仅计数模式）

```r
DGIdbr(mode = "group", ..., enrichment = FALSE)
```

## 结果解读

每次运行输出 `dgidb_hits.csv`，包含以下关键列：

| 列名 | 含义 |
|------|------|
| `drug` | 药物名称 |
| `gene_count` | 输入基因中有多少个与该药物互作 |
| `total_targets` | 该药物在 DGIdb 中的全部靶基因数 |
| `total_score` | DGIdb 互作分数总和 |
| `enrichment_ratio` | $(k/n) \div (m/N)$ — 相对于随机期望的富集倍数 |
| `p_value` | 超几何检验 p 值 |
| `fdr` | Benjamini–Hochberg 校正后的 p 值 |
| `significance` | `***` < 0.001, `**` < 0.01, `*` < 0.05, `.` < 0.1, `ns` |

结果按 FDR 升序排列。筛选高置信度命中：

```r
hits <- read.csv("path/to/output/dgidb_group/up/dgidb_hits.csv")
strict <- subset(hits, fdr < 0.05 & enrichment_ratio > 5)
```

## 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `mode` | `"group"` | `"group"` 或 `"subtype"` |
| `base_tables` | — | 输入 CSV 文件所在目录 |
| `base_out` | — | 输出目录（子文件夹自动创建） |
| `group_filename` | `"group.csv"` | 分组模式输入文件名 |
| `subtype_filename` | `"subtype.csv"` | 亚型模式输入文件名 |
| `approve` | `TRUE` | 仅保留 FDA 已批准药物 |
| `enrichment` | `TRUE` | 是否进行超几何富集分析（设为 `FALSE` 回退旧行为） |
| `background_N` | `NULL` | 背景基因数。`NULL` = 自动从 DGIdb 获取（~11,665）。可设数字覆盖。 |

## 输入文件格式

- CSV 格式，UTF-8 编码，含表头
- **分组模式**：列 `gene`（基因符号）、`direction`（`up` / `down`）
- **亚型模式**：列 `gene`（基因符号）、`direction`（`up` / `down`）、`subtype`（字符串）
- 预处理：去除空白、合并重复、使用官方 HGNC 基因符号

## 环境变量

设置 `DGIDB_URL` 可覆盖默认 GraphQL 端点：

```bash
export DGIDB_URL="https://custom-dgidb-instance.org/api/graphql"
```

## 注意事项

本工具用于**假设生成**，不可用于临床推荐。请注意：

1. **研究偏差**——研究越多的药物和基因在 DGIdb 中的记录越多。富集框架能惩罚泛靶点药物，
   但无法补充未被文献收录的互作关系。
2. **药物作用方向**——当前版本不区分抑制剂与激活剂。抑制一个下调的抑癌基因可能反而有害。
   请务必人工审核 Top 候选药物的机制合理性。
3. **互作分数的跨源可比性**——DGIdb 汇总了多个数据库的评分，各数据库量纲不同。富集统计仅
   使用互作的有无（二值），不受此影响，但 `total_score` 列应谨慎解读。

详见 `vignette("DGIdbr")` 获取完整讨论和工作示例。

## 引用

L. Zhang (2025). *DGIdbr: DGIdb gene set query helper.* R package version
1.2.0. https://github.com/lancelotzhang0124/DGIdbr

```bibtex
@manual{DGIdbr,
  author = {Zhang, L.},
  title  = {DGIdbr: DGIdb gene set query helper},
  year   = {2025},
  version = {1.2.0},
  url    = {https://github.com/lancelotzhang0124/DGIdbr}
}
```

## 维护说明

本项目最初由我的个人 GitHub 账号开发并维护：

https://github.com/lancelotzhang0124/DGIdbr

为更好地区分个人与学术相关工作，后续科研软件开发及项目维护已迁移至我的学术 GitHub 账号：

https://github.com/lfzhang00

未来的版本发布、功能更新及持续维护将通过该学术账号进行。


## 许可证

MIT。详见 [LICENSE](LICENSE)。
