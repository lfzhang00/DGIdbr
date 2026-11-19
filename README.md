<p align="right">
  <b>English</b> | <a href="README_zh.md">中文</a>
</p>

# DGIdbr

[![R CMD check](https://github.com/lfzhang00/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lfzhang00/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.3.2-blue)](https://github.com/lfzhang00/DGIdbr/releases/tag/1.3.2)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

Query DGIdb drug-gene interactions and rank drugs by hypergeometric enrichment.
Also search any drug via ChEMBL for its targets, mechanism of action, and
clinical indications.

## Installation

```r
remotes::install_github("lfzhang00/DGIdbr", ref = "1.3.2")
library(DGIdbr)
```

> Note: the legacy repository `lancelotzhang0124/DGIdbr` is archived and no longer maintained.

## Reading results in 1.3.2

```r
run_gene_set("SIRT5", "SIRT5", out_dir = "sirt5", approve = FALSE, enrichment = FALSE)
read.csv("sirt5/dgidb_raw.csv")[c("drug", "interaction_types", "approved")]
```

- A failed primary query (DNS, HTTP, GraphQL, or invalid JSON) raises an R error, not an empty result. Batch runs identify failed sets after processing the others. Use a fresh output directory: a failed query does not remove files from earlier runs.
- A successful query with no interactions writes two header-only CSVs. `dgidb_raw.csv` retains all interactions and their reported types (e.g., `inhibitor`); `approve = TRUE` filters only `dgidb_hits.csv`, keeping explicit DGIdb `approved = TRUE` values. An empty hits table with a nonempty raw table can therefore be a valid filter result.
- DGIdb reads use batches of up to 100 names and retry transient failures up to twice. Set `options(DGIdbr.timeout = 180)` to override the default 120 seconds per attempt. ChEMBL mechanism queries follow the returned `next` link; incomplete batches are not cached as complete.
- With direction scoring enabled, both CSVs are written before ChEMBL scoring and updated afterward. Set `direction = FALSE` to skip scoring without changing the up/down groups.

## Functions

### `DGIdbr()` — Gene-set driven drug prioritisation

Given a CSV of differentially expressed genes, queries [DGIdb](https://dgidb.org/)
and ranks drugs by enrichment (not hit count). If the CSV has a `direction` column
(`up`/`down`), ChEMBL mechanism-of-action data is **automatically** added to the
output.

```r
# Group mode: gene required, direction optional
DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv",
       base_out = ".", approve = TRUE, enrichment = TRUE)

# Subtype mode: CSV needs columns gene, direction, subtype
DGIdbr(mode = "subtype", base_tables = ".", subtype_filename = "subtype.csv",
       base_out = ".")

# Disable enrichment to fall back to counting
DGIdbr(mode = "group", ..., enrichment = FALSE)
```

A group CSV with only `gene` produces one result in `dgidb_group/all`. With a
`direction` column, results remain in `up` and `down`. The automatic enrichment
background counts DGIdb gene records and uses matched input symbols as the
sample. Unavailable background counts leave enrichment statistics missing. An explicit `background_N` uses all distinct input symbols. Conflicting
ChEMBL action directions for a drug–gene pair are reported as missing.

#### Output: `dgidb_hits.csv`

| Column | Description |
|--------|-------------|
| `drug`, `gene_count`, `total_targets`, `total_score` | Basic stats |
| `enrichment_ratio`, `p_value`, `fdr`, `significance` | Enrichment test |
| `n_with_direction_data`, `n_direction_consistent`, `direction_ratio` | ChEMBL direction scoring (auto-added when CSV has `direction`) |

#### Parameters

| Param | Default | Description |
|-------|---------|-------------|
| `mode` | `"group"` | `"group"` or `"subtype"` |
| `base_tables` | — | Directory with input CSVs |
| `base_out` | — | Output directory |
| `group_filename` | `"group.csv"` | Group-mode input file |
| `subtype_filename` | `"subtype.csv"` | Subtype-mode input file |
| `approve` | `TRUE` | Restrict the hits table to drugs marked approved by DGIdb |
| `direction` | `TRUE` | Score direction when up/down labels are available |
| `enrichment` | `TRUE` | Hypergeometric test with FDR |
| `background_N` | `NULL` | Background gene count (auto-detect if `NULL`) |

### `drug_card()` — Drug lookup via ChEMBL

Search any drug name for its targets, mechanism, and indications.

```r
card <- drug_card("ASPIRIN")
card$target_genes           # gene symbols (e.g. PTGS2)
card$mechanisms$action_type # e.g. "INHIBITOR"
card$indications$disease    # condition names

drug_card("ASPIRIN", phase = "approved")   # approved only
drug_card("ASPIRIN", phase = "trial")       # investigational only
```

| Param | Default | Description |
|-------|---------|-------------|
| `drug_name` | — | Drug name (e.g. `"ASPIRIN"`, `"CISPLATIN"`) |
| `phase` | `"all"` | Filter indications: `"all"`, `"approved"`, or `"trial"` |

## Environment

| Variable | Purpose | Default |
|----------|---------|---------|
| `DGIDB_URL` | Override DGIdb GraphQL endpoint | `https://dgidb.org/api/graphql` |
| `NO_PROXY` | Bypass proxy for ChEMBL | — |

## Citation

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

## License

MIT. See [LICENSE.md](LICENSE.md).
