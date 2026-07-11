<p align="right">
  <b>English</b> | <a href="README_zh.md">中文</a>
</p>

# DGIdbr

[![R CMD check](https://github.com/lancelotzhang0124/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lancelotzhang0124/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.2.0-blue)](https://github.com/lancelotzhang0124/DGIdbr)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

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
