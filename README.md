<p align="right">
  <b>English</b> | <a href="README_zh.md">中文</a>
</p>

# DGIdbr

[![R CMD check](https://github.com/lfzhang00/DGIdbr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lfzhang00/DGIdbr/actions)
[![Version](https://img.shields.io/badge/version-1.3.0-blue)](https://github.com/lfzhang00/DGIdbr)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

**Statistically grounded drug prioritisation from gene sets.**

Given a list of differentially expressed genes, DGIdbr queries the
[DGIdb](https://dgidb.org/) drug--gene interaction database and ranks drugs by
**hypergeometric enrichment** -- not by naive hit counting. It answers: *"Is
this drug's overlap with my gene set greater than expected by chance?"*

**v1.3+** extends the package with **ChEMBL-powered drug mechanism direction
scoring** and a **drug card lookup** -- you can now search a drug name and get
its targets, mechanism of action, and clinical indications.

## Why this matters

| Approach | What it does | Problem |
|----------|-------------|---------|
| **Naive counting** | Ranks drugs by how many of your genes they target | Promiscuous drugs (e.g. broad chemotherapies) always win |
| **DGIdbr (v1.2+)** | Hypergeometric test + FDR + enrichment ratio | Specific, mechanistically plausible drugs surface to the top |

**Example**: With 9 cancer genes as input, the old counting method ranks
RIBAVIRIN and CISPLATIN first -- both hit hundreds of genes. The enrichment
method instead ranks CHEMBL1214407 first: only 4 known targets, but 2 of them
are in your gene set -- a **648x enrichment** over random expectation (FDR = 0.001).

## Features

### Core (v1.0+)
- **Hypergeometric enrichment test** with Benjamini--Hochberg FDR correction
- **Dynamic background calibration** -- automatically queries DGIdb for the
  true count of druggable genes (~11,665), not a hardcoded guess
- **Enrichment ratio** as an intuitive effect-size metric
- **Group mode**: up/down gene sets from case-control differential expression
- **Subtype mode**: auto-detects all subtypes and builds up/down sets per subtype
- **FDA approval filter** -- keep only approved drugs or include investigational
- **Full backward compatibility** -- set `enrichment = FALSE` for the old
  counting-based behaviour

### New in v1.3.0 -- ChEMBL integration
- **Direction consistency scoring**: If your input CSV has a `direction` column
  (`up` / `down`), DGIdbr automatically queries the ChEMBL API for each drug's
  mechanism of action (inhibitor, activator, antagonist, etc.) and scores the
  consistency between drug action direction and gene expression direction.
  Added columns: `n_with_direction_data`, `n_direction_consistent`,
  `direction_ratio` in the output CSV.
- **Drug card lookup**: Search any drug name and get its ChEMBL molecule info,
  mechanism of action, resolved target genes, and clinical indications with
  development phases -- without needing an input gene set.

## Installation

> **Note**: Development has moved to [github.com/lfzhang00/DGIdbr](https://github.com/lfzhang00/DGIdbr).
> The legacy repository `lancelotzhang0124/DGIdbr` is **archived and no longer maintained**.

```r
# install.packages("remotes")
remotes::install_github("lfzhang00/DGIdbr")
```

## Quick start

### Group (case--control) mode

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

Direction scoring is **automatic** when the CSV contains a `direction` column.

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

### Drug card lookup (new in v1.3.0)

Search a drug by name and get its targets, mechanism, and indications:

```r
library(DGIdbr)

# Print a summary card
drug_card("ASPIRIN")

# Use the returned structured data
card <- drug_card("CISPLATIN")
card$target_genes             # resolved gene symbols
card$mechanisms$action_type   # e.g. "INHIBITOR"
card$indications$disease      # clinical indications
```

Example output:

```
== Drug Card: ASPIRIN ==
=> ChEMBL ID: CHEMBL25  Status: Approved

=> Targets & Mechanism of Action:
  * INHIBITOR (PTGS2) -- Cyclooxygenase inhibitor

=> Indications (top 15 of 49):
  * Fever  (Approved)
  * Myocardial Infarction  (Approved)
  * Pain  (Approved)
  * Stroke  (Approved)
  * Thrombosis  (Approved)
  * Atrial Fibrillation  (Phase 3)
  * Breast Neoplasms  (Phase 3)
  * ...
```

## Interpreting the output

### Main run output: `dgidb_hits.csv`

| Column | Meaning |
|--------|---------|
| `drug` | Drug name |
| `gene_count` | How many input genes interact with this drug |
| `total_targets` | Total genes this drug targets (from DGIdb) |
| `total_score` | Sum of DGIdb interaction scores |
| `enrichment_ratio` | (k/n) / (m/N) -- fold enrichment over random expectation |
| `p_value` | Hypergeometric p-value |
| `fdr` | Benjamini--Hochberg corrected p-value |
| `significance` | `***` < 0.001, `**` < 0.01, `*` < 0.05, `.` < 0.1, `ns` |
| `n_with_direction_data` | Number of (drug, target) pairs with ChEMBL mechanism data |
| `n_direction_consistent` | Pairs where drug action direction matches gene expression direction |
| `direction_ratio` | Direction consistency rate (consistent / annotated total) |

Results are sorted by FDR ascending. Filter for high-confidence hits:

```r
hits <- read.csv("path/to/output/dgidb_group/up/dgidb_hits.csv")
strict <- subset(hits, fdr < 0.05 & enrichment_ratio > 5)
```

### Drug card output: `drug_card()`

Returns a named list:
- `$drug` -- drug name
- `$molecule_chembl_id` -- ChEMBL molecule ID
- `$max_phase` -- highest clinical phase reached
- `$mechanisms` -- data.frame of mechanisms (target_chembl_id, action_type, mechanism_of_action)
- `$target_genes` -- named character vector (target_chembl_id -> gene symbol)
- `$indications` -- data.frame of indications with phase labels

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mode` | `"group"` | `"group"` or `"subtype"` |
| `base_tables` | -- | Directory containing input CSV files |
| `base_out` | -- | Output directory (subfolders auto-created) |
| `group_filename` | `"group.csv"` | Group-mode input file name |
| `subtype_filename` | `"subtype.csv"` | Subtype-mode input file name |
| `approve` | `TRUE` | Keep only FDA-approved drugs |
| `enrichment` | `TRUE` | Run hypergeometric enrichment (set `FALSE` for old behaviour) |
| `background_N` | `NULL` | Background gene count. `NULL` = auto-detect from DGIdb (~11,665). Set a number to override. |

**Direction scoring is automatic** -- no extra parameter needed. When the input
CSV contains a `direction` column, the ChEMBL API is queried and direction
columns are added to the output. If ChEMBL is unreachable (no internet,
rate-limited, or temporary API outage), these columns are silently filled with
`NA` and the core analysis continues without interruption.

## Input file format

- CSV with header, UTF-8 encoded
- **Group mode**: columns `gene` (symbol), `direction` (`up` / `down`)
- **Subtype mode**: columns `gene` (symbol), `direction` (`up` / `down`), `subtype` (string)
- Clean genes beforehand: remove blanks, collapse duplicates, use official HGNC symbols

## Environment variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `DGIDB_URL` | Override the DGIdb GraphQL endpoint | `https://dgidb.org/api/graphql` |
| `NO_PROXY` | Bypass proxy for ChEMBL API (if behind HTTP proxy) | -- |

If your network requires a proxy but ChEMBL API fails, try:

```bash
export NO_PROXY="www.ebi.ac.uk,ebi.ac.uk"
```

## Caveats

This tool is for **hypothesis generation**, not clinical recommendation.
Please be aware of:

1. **Study bias** -- well-studied drugs and genes have more documented
   interactions in DGIdb. The enrichment framework penalises promiscuity but
   cannot recover interactions that haven't been curated yet.

2. **Drug mechanism direction** -- v1.3+ queries ChEMBL automatically for
   drug mechanism of action (inhibitor / activator / antagonist etc.) and
   scores consistency with gene expression direction (up/down). The
   `direction_ratio` column quantifies this alignment. Note that:
   - ChEMBL coverage varies; drugs without mechanism data show `NA`
   - The classification (suppress vs. enhance) is a heuristic;
     context-dependent effects (e.g. tissue-specific signalling, oncogene vs.
     tumour suppressor) should always be reviewed manually for top candidates
   - If ChEMBL is unavailable, direction columns are filled with `NA` without
     affecting the core analysis

3. **Interaction score comparability** -- DGIdb aggregates scores from
   multiple databases with different scales. The enrichment statistics use
   binary presence/absence and are unaffected, but the `total_score` column
   should be interpreted cautiously.

See `vignette("DGIdbr")` for a full discussion with worked examples.

## Citation

L. Zhang (2025). *DGIdbr: DGIdb gene set query helper.* R package version
1.3.0. https://github.com/lfzhang00/DGIdbr

```bibtex
@manual{DGIdbr,
  author = {Zhang, L.},
  title  = {DGIdbr: DGIdb gene set query helper},
  year   = {2025},
  version = {1.3.0},
  url    = {https://github.com/lfzhang00/DGIdbr}
}
```

## Maintenance notice

This package is maintained at [github.com/lfzhang00/DGIdbr](https://github.com/lfzhang00/DGIdbr).
The legacy repository `lancelotzhang0124/DGIdbr` is **archived and no longer maintained**.
All future releases will only be published under the academic profile.

## Version history

| Version | New functions | Key changes |
|---------|---------------|-------------|
| **1.3.0** | `drug_card()` | ChEMBL integration: direction consistency scoring, drug card lookup, automatic direction columns (`n_with_direction_data`, `n_direction_consistent`, `direction_ratio`). Internal: `chembl_map_drug_to_molecule()`, `chembl_fetch_mechanisms()`, `chembl_target_to_gene()`, `chembl_fetch_indications()`, `classify_action_direction()`, `compute_direction_consistency()`. |
| **1.2.0** | — | Hypergeometric enrichment test with FDR correction, enrichment ratio, dynamic background calibration. Internal: `fetch_drug_target_counts()`, `fetch_druggable_gene_count()`, `compute_drug_enrichment()`. |
| **1.1.0** | — | Subtype mode, FDA approval filter, grouped output directory structure. |
| **1.0.0** | `DGIdbr()`, `run_gene_set()`, `build_gene_sets()` | Initial release: DGIdb GraphQL queries, drug-gene interaction aggregation, group mode (up/down). |

## License
## License

MIT. See [LICENSE](LICENSE).
