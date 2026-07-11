<p align="right">
  <a href="#english"><b>English</b></a> |
  <a href="#chinese--中文">中文</a>
</p>

# DGIdbr

[![R CMD check](https://github.com/lancelotzhang0124/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lancelotzhang0124/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.2.0-blue)](https://github.com/lancelotzhang0124/DGIdbr)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

---

<div id="english"></div>

**Statistically grounded drug prioritisation from gene sets.**

Given a list of differentially expressed genes, DGIdbr queries the
[DGIdb](https://dgidb.org/) drug–gene interaction database and ranks drugs by
**hypergeometric enrichment** — not by naïve hit counting. It answers: *"Is
this drug's overlap with my gene set greater than expected by chance?"*

## Why this matters

| Approach | What it does | Problem |
|----------|-------------|---------|
| **Naïve counting** | Ranks drugs by how many of your genes they target | Promiscuous drugs (e.g. broad chemotherapies) always win |
| **DGIdbr (v1.2+)** | Hypergeometric test + FDR + enrichment ratio | Specific, mechanistically plausible drugs surface to the top |

**Example**: With 9 cancer genes as input, the old counting method ranks
RIBAVIRIN and CISPLATIN first — both hit hundreds of genes. The enrichment
method instead ranks CHEMBL1214407 first: only 4 known targets, but 2 of them
are in your gene set — a **648× enrichment** over random expectation (FDR = 0.001).

## Features

- **Hypergeometric enrichment test** with Benjamini–Hochberg FDR correction
- **Dynamic background calibration** — automatically queries DGIdb for the
  true count of druggable genes (~11,665), not a hardcoded guess
- **Enrichment ratio** as an intuitive effect-size metric
- **Group mode**: up/down gene sets from case-control differential expression
- **Subtype mode**: auto-detects all subtypes and builds up/down sets per subtype
- **FDA approval filter** — keep only approved drugs or include investigational
- **Full backward compatibility** — set `enrichment = FALSE` for the old
  counting-based behaviour

## Installation

```r
# install.packages("remotes")
remotes::install_github("lancelotzhang0124/DGIdbr")
```

## Quick start

### Group (case–control) mode

Input CSV needs columns `gene` and `direction` (`up` / `down`):

```r
library(DGIdbr)

DGIdbr(
  mode           = "group",
  base_tables    = "path/to/input",
  group_filename = "group.csv",
  base_out       = "path/to/output",
  approve        = TRUE,       # FDA-approved only
  enrichment     = TRUE         # hypergeometric test + FDR (default)
)
```

### Subtype mode

Input CSV needs columns `gene`, `direction`, and `subtype`:

```r
DGIdbr(
  mode             = "subtype",
  base_tables      = "path/to/input",
  subtype_filename = "subtype.csv",
  base_out         = "path/to/output",
  approve          = TRUE
)
```

### Backward-compatible (counting only)

```r
DGIdbr(mode = "group", ..., enrichment = FALSE)
```

## Interpreting the output

Each run writes `dgidb_hits.csv` with these key columns:

| Column | Meaning |
|--------|---------|
| `drug` | Drug name |
| `gene_count` | How many input genes interact with this drug |
| `total_targets` | Total genes this drug targets (from DGIdb) |
| `total_score` | Sum of DGIdb interaction scores |
| `enrichment_ratio` | $(k/n) \div (m/N)$ — fold enrichment over random expectation |
| `p_value` | Hypergeometric p-value |
| `fdr` | Benjamini–Hochberg corrected p-value |
| `significance` | `***` < 0.001, `**` < 0.01, `*` < 0.05, `.` < 0.1, `ns` |

Results are sorted by FDR ascending. Filter for high-confidence hits:

```r
hits <- read.csv("path/to/output/dgidb_group/up/dgidb_hits.csv")
strict <- subset(hits, fdr < 0.05 & enrichment_ratio > 5)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mode` | `"group"` | `"group"` or `"subtype"` |
| `base_tables` | — | Directory containing input CSV files |
| `base_out` | — | Output directory (subfolders auto-created) |
| `group_filename` | `"group.csv"` | Group-mode input file name |
| `subtype_filename` | `"subtype.csv"` | Subtype-mode input file name |
| `approve` | `TRUE` | Keep only FDA-approved drugs |
| `enrichment` | `TRUE` | Run hypergeometric enrichment (set `FALSE` for old behaviour) |
| `background_N` | `NULL` | Background gene count. `NULL` = auto-detect from DGIdb (~11,665). Set a number to override. |

## Input file format

- CSV with header, UTF-8 encoded
- **Group mode**: columns `gene` (symbol), `direction` (`up` / `down`)
- **Subtype mode**: columns `gene` (symbol), `direction` (`up` / `down`), `subtype` (string)
- Clean genes beforehand: remove blanks, collapse duplicates, use official HGNC symbols

## Environment

Set `DGIDB_URL` to override the default GraphQL endpoint:

```bash
export DGIDB_URL="https://custom-dgidb-instance.org/api/graphql"
```

## Caveats

This tool is for **hypothesis generation**, not clinical recommendation.
Please be aware of:

1. **Study bias** — well-studied drugs and genes have more documented
   interactions in DGIdb. The enrichment framework penalises promiscuity but
   cannot recover interactions that haven't been curated yet.
2. **Drug mechanism** — the current version does not distinguish inhibitor
   from activator. A drug that inhibits a down-regulated tumour suppressor
   could be harmful. Always review top candidates manually.
3. **Interaction score comparability** — DGIdb aggregates scores from
   multiple databases with different scales. The enrichment statistics use
   binary presence/absence and are unaffected, but the `total_score` column
   should be interpreted cautiously.

See `vignette("DGIdbr")` for a full discussion with worked examples.

## Citation

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

## License

MIT. See [LICENSE](LICENSE).

---

<div id="chinese--中文"></div>

<p align="right">
  <a href="#english"><b>English</b></a> |
  <a href="#chinese--中文">中文</a>
</p>

# DGIdbr <small>中文</small>

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

### 分组模式

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

### 向后兼容（仅计数）

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

## 许可证

MIT。详见 [LICENSE](LICENSE)。
