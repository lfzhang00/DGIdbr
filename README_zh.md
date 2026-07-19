<p align="right">
  <a href="README.md">English</a> | <b>中文</b>
</p>

# DGIdbr

[![R CMD check](https://github.com/lfzhang00/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lfzhang00/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.3.1-blue)](https://github.com/lfzhang00/DGIdbr)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

查询 DGIdb 药物-基因互作数据库，通过超几何富集检验对药物排序。
也可直接搜索药物名称，从 ChEMBL 获取其靶点、作用机制和临床适应症。

## 安装

```r
remotes::install_github("lfzhang00/DGIdbr")
```

> 注意：原仓库 `lancelotzhang0124/DGIdbr` 已归档，不再维护。

## 函数说明

### `DGIdbr()` — 基因集驱动的药物排序

输入差异表达基因 CSV，查询 [DGIdb](https://dgidb.org/) 并用富集检验排序（而非简单计数）。
CSV 包含 `direction` 列（`up`/`down`）时**自动**触发 ChEMBL 方向一致性评分。

```r
# 分组模式：CSV 需包含 gene, direction 列
DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv",
       base_out = ".", approve = TRUE, enrichment = TRUE)

# 亚型模式：CSV 需包含 gene, direction, subtype 列
DGIdbr(mode = "subtype", base_tables = ".", subtype_filename = "subtype.csv",
       base_out = ".")

# 关闭富集，回退到计数模式
DGIdbr(mode = "group", ..., enrichment = FALSE)
```

#### 输出文件 `dgidb_hits.csv`

| 列名 | 说明 |
|------|------|
| `drug`, `gene_count`, `total_targets`, `total_score` | 基本统计 |
| `enrichment_ratio`, `p_value`, `fdr`, `significance` | 富集检验 |
| `n_with_direction_data`, `n_direction_consistent`, `direction_ratio` | ChEMBL 方向评分（CSV 有 direction 列时自动添加）|

#### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `mode` | `"group"` | `"group"` 或 `"subtype"` |
| `base_tables` | — | 输入 CSV 目录 |
| `base_out` | — | 输出目录 |
| `group_filename` | `"group.csv"` | 分组模式输入文件 |
| `subtype_filename` | `"subtype.csv"` | 亚型模式输入文件 |
| `approve` | `TRUE` | 仅保留 FDA 已批准药物 |
| `enrichment` | `TRUE` | 超几何富集检验 + FDR |
| `background_N` | `NULL` | 背景基因数（`NULL` 自动获取）|

### `drug_card()` — 药物查询

搜索药名，获取靶点、作用机制和适应症。

```r
card <- drug_card("ASPIRIN")
card$target_genes           # 基因符号（如 PTGS2）
card$mechanisms$action_type # 如 "INHIBITOR"
card$indications$disease    # 适应症名称

drug_card("ASPIRIN", phase = "approved")   # 仅已批准
drug_card("ASPIRIN", phase = "trial")       # 仅试验中
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `drug_name` | — | 药物名称（如 `"ASPIRIN"`、`"CISPLATIN"`） |
| `phase` | `"all"` | 适应症筛选：`"all"`、`"approved"`、`"trial"` |

## 环境变量

| 变量 | 用途 | 默认值 |
|------|------|--------|
| `DGIDB_URL` | 覆盖 DGIdb GraphQL 端点 | `https://dgidb.org/api/graphql` |
| `NO_PROXY` | 为 ChEMBL 绕过代理 | — |

## 引用

L. Zhang (2025). *DGIdbr: DGIdb gene set query helper.* R package version 1.3.1. https://github.com/lfzhang00/DGIdbr

```bibtex
@manual{DGIdbr,
  author = {Zhang, L.},
  title  = {DGIdbr: DGIdb gene set query helper},
  year   = {2025},
  version = {1.3.1},
  url    = {https://github.com/lfzhang00/DGIdbr}
}
```

## 许可证

MIT。详见 [LICENSE](LICENSE)。
