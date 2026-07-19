<p align="right">
  <a href="README.md">English</a> | <b>中文</b>
</p>

# DGIdbr

[![R CMD check](https://github.com/lancelotzhang0124/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lancelotzhang0124/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.3.0-blue)](https://github.com/lancelotzhang0124/DGIdbr)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

**基于统计检验的药物优先排序工具。**

给定一组差异表达基因，DGIdbr 查询 [DGIdb](https://dgidb.org/) 药物-基因互作数据库，
并通过**超几何富集检验**（而非简单的命中计数）对药物进行排序。它回答的核心问题是：
*"这个药物与我的基因集的重叠是否显著高于随机期望？"*

**v1.3+** 新增基于 **ChEMBL 的药物作用方向评分**和**药物卡片查询**功能——现在你
可以直接搜索一个药名，获取它的靶点、作用机制和临床适应症。

## 为什么重要

| 方法 | 做法 | 缺陷 |
|------|------|------|
| **简单计数** | 按药物命中了多少个输入基因排序 | 泛靶点药物（如广谱化疗药）永远排第一 |
| **DGIdbr (v1.2+)** | 超几何检验 + FDR + 富集倍数 | 特异性强、机制合理的药物脱颖而出 |

**示例**：输入 9 个癌症相关基因，旧版计数法将 RIBAVIRIN 和 CISPLATIN 排在最前面
——它们都靶向数百个基因。富集法将 CHEMBL1214407 排在第一：该化合物仅有 4 个已知靶点，
但其中 2 个在输入基因集中——**648 倍富集**（FDR = 0.001）。

## 功能特性

### 核心功能 (v1.0+)
- **超几何富集检验**，带 Benjamini--Hochberg FDR 校正
- **动态背景校准**——自动查询 DGIdb 当前可成药基因总数（~11,665），而非硬编码猜测值
- **富集倍数**（enrichment ratio）作为直观的效应量指标
- **分组模式**：从 case-control 差异表达表中提取 up/down 基因集
- **亚型模式**：自动识别所有亚型，为每个亚型构建 up/down 基因集
- **FDA 批准过滤**——仅保留已批准药物，或纳入研究性药物
- **完全向后兼容**——设置 `enrichment = FALSE` 即可回退到旧的计数模式

### v1.3.0 新增 —— ChEMBL 集成
- **药物作用方向一致性评分**：如果输入 CSV 包含 `direction` 列（`up` / `down`），
  DGIdbr 自动查询 ChEMBL API 获取每个药物的作用机制（抑制剂、激活剂、拮抗剂等），
  并与基因表达方向进行一致性评分。输出 CSV 新增三列：
  `n_with_direction_data`、`n_direction_consistent`、`direction_ratio`。
- **药物卡片查询**：无需输入基因集，直接搜索药名即可获取其 ChEMBL 分子信息、
  作用机制、解析后的靶基因符号以及带临床阶段的适应症列表。

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

当 CSV 包含 `direction` 列时，**方向评分自动触发**。

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

### 药物卡片查询（v1.3.0 新增）

直接搜药名，获取靶点、机制和适应症：

```r
library(DGIdbr)

# 打印摘要卡片
drug_card("ASPIRIN")

# 使用返回的结构化数据
card <- drug_card("CISPLATIN")
card$target_genes             # 解析后的基因符号
card$mechanisms$action_type   # 例如 "INHIBITOR"
card$indications$disease      # 临床适应症
```

输出示例：

```
== Drug Card: ASPIRIN ==
=> ChEMBL ID: CHEMBL25  状态: Approved

=> 靶点与作用机制:
  * INHIBITOR (PTGS2) -- Cyclooxygenase inhibitor

=> 适应症（前 15 / 共 49）:
  * Fever  (Approved)
  * Myocardial Infarction  (Approved)
  * Pain  (Approved)
  * Stroke  (Approved)
  * Thrombosis  (Approved)
  * Atrial Fibrillation  (Phase 3)
  * Breast Neoplasms  (Phase 3)
  * ...
```

## 结果解读

### 主运行输出: `dgidb_hits.csv`

| 列名 | 含义 |
|------|------|
| `drug` | 药物名称 |
| `gene_count` | 输入基因中有多少个与该药物互作 |
| `total_targets` | 该药物在 DGIdb 中的全部靶基因数 |
| `total_score` | DGIdb 互作分数总和 |
| `enrichment_ratio` | (k/n) / (m/N) -- 相对于随机期望的富集倍数 |
| `p_value` | 超几何检验 p 值 |
| `fdr` | Benjamini--Hochberg 校正后的 p 值 |
| `significance` | `***` < 0.001, `**` < 0.01, `*` < 0.05, `.` < 0.1, `ns` |
| `n_with_direction_data` | 该药物有 ChEMBL 机制注释的靶基因数 |
| `n_direction_consistent` | 药物作用方向与基因表达方向一致的靶基因数 |
| `direction_ratio` | 方向一致率（consistent / 有注释的总数）|

结果按 FDR 升序排列。筛选高置信度命中：

```r
hits <- read.csv("path/to/output/dgidb_group/up/dgidb_hits.csv")
strict <- subset(hits, fdr < 0.05 & enrichment_ratio > 5)
```

### 药物卡片输出: `drug_card()`

返回一个命名列表：
- `$drug` -- 药物名称
- `$molecule_chembl_id` -- ChEMBL 分子 ID
- `$max_phase` -- 最高临床阶段
- `$mechanisms` -- 机制数据框（target_chembl_id, action_type, mechanism_of_action）
- `$target_genes` -- 命名向量（target_chembl_id -> 基因符号）
- `$indications` -- 适应症数据框（含阶段标签）

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

**方向评分自动启用**——无需额外参数。当 CSV 包含 `direction` 列时自动查询 ChEMBL
并添加方向列到输出。如果 ChEMBL 不可用（无网络、限流或服务暂时中断），这些列会
以 `NA` 填充，核心分析不受影响。

## 输入文件格式

- CSV 格式，UTF-8 编码，含表头
- **分组模式**：列 `gene`（基因符号）、`direction`（`up` / `down`）
- **亚型模式**：列 `gene`（基因符号）、`direction`（`up` / `down`）、`subtype`（字符串）
- 预处理：去除空白、合并重复、使用官方 HGNC 基因符号

## 环境变量

| 变量 | 用途 | 默认值 |
|------|------|--------|
| `DGIDB_URL` | 覆盖 DGIdb GraphQL 端点 | `https://dgidb.org/api/graphql` |
| `NO_PROXY` | 为 ChEMBL API 绕过代理（如使用 HTTP 代理） | — |

如果网络需要通过代理但 ChEMBL API 连接失败，试试：

```bash
export NO_PROXY="www.ebi.ac.uk,ebi.ac.uk"
```

## 注意事项

本工具用于**假设生成**，不可用于临床推荐。请注意：

1. **研究偏差**——研究越多的药物和基因在 DGIdb 中的记录越多。富集框架能惩罚泛靶点药物，
   但无法补充未被文献收录的互作关系。

2. **药物作用方向**——v1.3+ 支持通过 ChEMBL API 自动查询药物作用机制
   （抑制剂/激活剂/拮抗剂等），并与基因表达方向（上调/下调）进行一致性评分。
   结果中的 `direction_ratio` 列量化了匹配程度。需要注意：
   - ChEMBL 覆盖范围有限，无机制数据的药物会以 `NA` 填充
   - 抑制/激活分类是启发式的，对于 Top 候选药物，仍需人工审核其
     在特定组织环境和基因类型（癌基因 vs 抑癌基因）下的机制合理性
   - ChEMBL 不可用时，方向列以 `NA` 填充，不影响核心分析流程

3. **互作分数的跨源可比性**——DGIdb 汇总了多个数据库的评分，各数据库量纲不同。
   富集统计仅使用互作的有无（二值），不受此影响，但 `total_score` 列应谨慎解读。

详见 `vignette("DGIdbr")` 获取完整讨论和工作示例。

## 引用

L. Zhang (2025). *DGIdbr: DGIdb gene set query helper.* R package version
1.3.0. https://github.com/lancelotzhang0124/DGIdbr

```bibtex
@manual{DGIdbr,
  author = {Zhang, L.},
  title  = {DGIdbr: DGIdb gene set query helper},
  year   = {2025},
  version = {1.3.0},
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
