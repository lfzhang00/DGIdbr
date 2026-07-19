<p align="right">
  <b>English</b> | <a href="README_zh.md">中文</a>
</p>

# DGIdbr

[![R CMD check](https://github.com/lfzhang00/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lfzhang00/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.3.1-blue)](https://github.com/lfzhang00/DGIdbr)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

Query DGIdb drug-gene interactions and rank drugs by hypergeometric enrichment.
Also search any drug via ChEMBL for its targets, mechanism of action, and
clinical indications.

## Installation

```r
remotes::install_github("lfzhang00/DGIdbr")
```

> Note: the legacy repository `lancelotzhang0124/DGIdbr` is archived and no longer maintained.

## Functions

### `DGIdbr()` — Gene-set driven drug prioritisation

Given a CSV of differentially expressed genes, queries [DGIdb](https://dgidb.org/)
and ranks drugs by enrichment (not hit count). If the CSV has a `direction` column
(`up`/`down`), ChEMBL mechanism-of-action data is **automatically** added to the
output.

```r
# Group mode: CSV needs columns gene, direction
DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv",
       base_out = ".", approve = TRUE, enrichment = TRUE)

# Subtype mode: CSV needs columns gene, direction, subtype
DGIdbr(mode = "subtype", base_tables = ".", subtype_filename = "subtype.csv",
       base_out = ".")

# Disable enrichment to fall back to counting
DGIdbr(mode = "group", ..., enrichment = FALSE)
```

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
| `approve` | `TRUE` | Only FDA-approved drugs |
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

## License

MIT. See [LICENSE](LICENSE).
