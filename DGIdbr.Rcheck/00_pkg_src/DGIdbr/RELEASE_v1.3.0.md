# DGIdbr v1.3.0 -- ChEMBL 集成：药物方向评分 + 药物卡片

## 概览

v1.3.0 将 ChEMBL API 集成为包的一等公民，无需用户额外配置即可自动获取
药物作用机制和适应症数据。新增两个核心能力：

1. **方向一致性评分**（自动）-- 输入 CSV 只要有 `direction` 列，自动
  查询 ChEMBL 并计算方向匹配度
2. **药物卡片查询**（交互）-- 搜药名就能看靶点、机制、适应症

## 新增内容

### 新功能

- `drug_card(drug_name)` -- 公开函数。搜索药物，返回靶基因、作用机制和
  临床适应症。输出可读摘要 + 结构化列表。

- **自动方向评分** -- `run_gene_set()` 新增 `gene_direction` 参数，
  `build_gene_sets()` 自动从 CSV 提取方向信息。当 CSV 含 `direction`
  列时，输出 `dgidb_hits.csv` 自动新增三列：
  - `n_with_direction_data` -- 有 ChEMBL 注释的靶基因数
  - `n_direction_consistent` -- 方向一致的靶基因数
  - `direction_ratio` -- 方向一致率

### 底层新增（内部函数）

| 函数 | 用途 |
|------|------|
| `chembl_map_drug_to_molecule()` | 药物名 → ChEMBL 分子 ID（批量、缓存） |
| `chembl_fetch_mechanisms()` | 获取药物作用机制（action_type + 靶点） |
| `chembl_target_to_gene()` | ChEMBL 靶点 ID → HGNC 基因符号 |
| `chembl_fetch_indications()` | 获取药物适应症（含临床阶段） |
| `classify_action_direction()` | 将 action_type 分类为 suppress/enhance |
| `compute_direction_consistency()` | 核心：计算方向一致性评分 |

### 文件变更

| 文件 | 状态 | 说明 |
|------|------|------|
| `R/direction.R` | **新增** | ChEMBL 客户端 + 方向评分模块 (~590行) |
| `R/DGIdbr.R` | 修改 | `build_gene_sets()` 添加方向元数据；`run_gene_set()` 集成方向分析 |
| `NAMESPACE` | 修改 | 新增 `export(drug_card)` |
| `DESCRIPTION` | 修改 | 版本号 1.2.0 → 1.3.0 |
| `README.md` / `README_zh.md` | 修改 | 完整更新，新增 ChEMBL 文档和环境变量说明 |
| `man/*.Rd` | 修改 | roxygen2 重新生成全部手册页 |

## 使用方法

```r
library(DGIdbr)

# 用法 1：基因集 → 药物排序（方向评分自动触发）
DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv")

# 用法 2：搜药
card <- drug_card("ASPIRIN")
card$target_genes       # PTGS2
card$indications        # 49 条适应症
```

## 行为变更

- 输入 CSV 含 `direction` 列时，**方向分析自动启用**（无侵入式，ChEMBL
  不可用时优雅降级）
- `dgidb_hits.csv` 新增三列方向数据（若 ChEMBL 不可用则为 NA）
- 输出文件 "dgidb_raw.csv" 新增 `direction_match` 列

## 向后兼容性

- **完全向后兼容**：没有 `direction` 列的旧版 CSV 行为不变
- `enrichment = FALSE` 计数模式仍可使用
- `DGIdbr()` 所有已有参数签名未变

## 依赖

- 新增运行时依赖：ChEMBL REST API（`https://www.ebi.ac.uk/chembl/api/data`）
- R 包依赖无变化（已有：httr, jsonlite, utils）
- 如果网络使用 HTTP 代理且 ChEMBL 访问失败，设置环境变量：
  ```bash
  export NO_PROXY="www.ebi.ac.uk,ebi.ac.uk"
  ```

## 已知限制

- ChEMBL 的药物机制覆盖不完整（如 METFORMIN 无机制数据），此类药物
  的方向列会以 NA 填充
- 某些靶点是非蛋白实体（如 CISPLATIN 靶向 DNA），无法解析为基因符号
- R CMD check 已知警告（..Rcheck 目录、缺少 inst/doc）与 v1.2.0 一致
