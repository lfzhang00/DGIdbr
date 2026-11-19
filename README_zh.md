<p align="right">
  <a href="README.md">English</a> | <b>中文</b>
</p>

# DGIdbr

[![R CMD check](https://github.com/lfzhang00/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lfzhang00/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.3.2-blue)](https://github.com/lfzhang00/DGIdbr/releases/tag/1.3.2)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

查询 DGIdb 药物-基因互作数据库，通过超几何富集检验对药物排序。
也可直接搜索药物名称，从 ChEMBL 获取其靶点、作用机制和临床适应症。

## 安装

```r
remotes::install_github("lfzhang00/DGIdbr", ref = "1.3.2")
library(DGIdbr)
```

> 注意：原仓库 `lancelotzhang0124/DGIdbr` 已归档，不再维护。

## 如何判断 1.3.2 的查询结果

```r
run_gene_set("SIRT5", "SIRT5", out_dir = "sirt5", approve = FALSE, enrichment = FALSE)
read.csv("sirt5/dgidb_raw.csv")[c("drug", "interaction_types", "approved")]
```

- 主查询遇到 DNS、HTTP、GraphQL 或 JSON 解析错误时会向 R 调用方报错，不会将失败写成零命中。批处理会继续其他基因集，最后汇总失败项。建议每次使用新输出目录：失败不会删除旧结果。
- 查询成功但没有互作时，两张 CSV 都只有表头。`dgidb_raw.csv` 保留全部互作和数据库报告的类型（如 `inhibitor`）；`approve = TRUE` 只筛选 `dgidb_hits.csv`，仅保留 DGIdb 明确标记 `approved = TRUE` 的药物。因此，原始表非空但汇总表为空，也可能是有效筛选结果。
- DGIdb 每批最多查询 100 个名称，临时错误最多重试两次。每次请求默认超时 120 秒，可用 `options(DGIdbr.timeout = 180)` 调整。ChEMBL 机制查询按返回的 `next` 链接翻页，未取完的批次不会缓存为完整结果。
- 启用方向评分时，两张 CSV 在 ChEMBL 评分前写出，完成后更新。设 `direction = FALSE` 可跳过评分，同时保留上下调分组。

## 函数说明

### `DGIdbr()` — 基因集驱动的药物排序

输入差异表达基因 CSV，查询 [DGIdb](https://dgidb.org/) 并用富集检验排序（而非简单计数）。
CSV 包含 `direction` 列（`up`/`down`）时**自动**触发 ChEMBL 方向一致性评分。

```r
# 分组模式：gene 为必填列，direction 可选
DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv",
       base_out = ".", approve = TRUE, enrichment = TRUE)

# 亚型模式：CSV 需包含 gene, direction, subtype 列
DGIdbr(mode = "subtype", base_tables = ".", subtype_filename = "subtype.csv",
       base_out = ".")

# 关闭富集，回退到计数模式
DGIdbr(mode = "group", ..., enrichment = FALSE)
```

分组 CSV 只有 `gene` 列时，结果写入 `dgidb_group/all`。提供 `direction` 后按
`up` 和 `down` 分组。自动背景以 DGIdb 中匹配到的输入基因作为样本，指定
`background_N` 时使用全部去重输入基因。自动背景查询失败时保留命中记录，富集统计量记为缺失。某药对同一基因的作用方向冲突时记为缺失。


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
| `approve` | `TRUE` | 汇总表仅保留 DGIdb 明确标记为已获批的药物 |
| `direction` | `TRUE` | 有上下调标签时计算方向一致性 |
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

L. Zhang (2026). *DGIdbr: DGIdb gene set query helper.* R package version 1.3.2. https://github.com/lfzhang00/DGIdbr

```bibtex
@manual{DGIdbr,
  author = {Zhang, L.},
  title  = {DGIdbr: DGIdb gene set query helper},
  year   = {2026},
  version = {1.3.2},
  url    = {https://github.com/lfzhang00/DGIdbr}
}
```

## 许可证

MIT。详见 [LICENSE.md](LICENSE.md)。
