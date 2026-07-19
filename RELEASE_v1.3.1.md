# DGIdbr v1.3.1 — Drug Card Phase Filter + ChEMBL Integration

## Highlights

**v1.3.1** builds on the ChEMBL integration with a new phase filter for
indications, plus all the foundational work from v1.3.0.

### v1.3.1 — Drug card phase filter

`drug_card()` gains a `phase` parameter to filter indications by development
stage:

```r
drug_card("ASPIRIN")                              # all 49 indications
drug_card("ASPIRIN", phase = "approved")           # Phase 4 / approved only
drug_card("ASPIRIN", phase = "trial")              # Phase 1-3 investigational
```

### v1.3.0 — ChEMBL Integration (first release)

Two new capabilities with zero user configuration:

1. **Direction consistency scoring** (automatic) — when your input CSV has a
   `direction` column, DGIdbr automatically queries ChEMBL for each drug's
   mechanism of action and scores how well the drug's action direction
   (inhibit vs. activate) aligns with the gene's expression direction
   (up vs. down).
2. **Drug card lookup** (interactive) — search any drug name to instantly get
   its targets, mechanism of action, and clinical indications with
   development phases.

## What's new

### New / updated functions

| Function | Description |
|----------|-------------|
| `drug_card(drug_name, phase)` | **Updated.** `phase` parameter (`"all"`, `"approved"`, `"trial"`) filters indications by development stage. Default `"all"`. |

### Automatic direction scoring

- `run_gene_set()` gains a new `gene_direction` parameter (internal, populated
  automatically from the input CSV)
- `build_gene_sets()` now extracts per-gene direction metadata from the input
- Output CSV `dgidb_hits.csv` now includes three new columns:
  - `n_with_direction_data` — how many of this drug's targets have ChEMBL mechanism annotations
  - `n_direction_consistent` — how many targets have action direction consistent with expression direction
  - `direction_ratio` — consistent / annotated (proportion of matches)

### New internal functions

| Function | Purpose |
|----------|---------|
| `chembl_map_drug_to_molecule()` | Drug name → ChEMBL molecule ID (batched, cached) |
| `chembl_fetch_mechanisms()` | Fetch mechanism of action records (action_type + target) |
| `chembl_target_to_gene()` | Resolve ChEMBL target IDs to HGNC gene symbols |
| `chembl_fetch_indications()` | Fetch clinical indications with development phases |
| `classify_action_direction()` | Map action types (INHIBITOR, AGONIST, …) to suppress/enhance |
| `compute_direction_consistency()` | Core scoring function for direction consistency |

### Bug fixes (v1.3.0)

- Fixed drug name lookup: replaced broken `/drug.json?name__in=` with
  `/molecule.json?pref_name__in=` (the drug endpoint's `name` field is always
  `None`)
- Fixed pagination metadata key: ChEMBL uses `page_meta` (not `page_metadata`)
  with a `next` field
- Fixed target synonym field names: ChEMBL returns `component_synonym` and
  `syn_type` (not `synonym` and `synonym_type`)
- Removed duplicate scoring block left over from a partial edit

### Documentation

- `README.md` / `README_zh.md` fully updated: new features, environment
  variables, direction scoring explanation, version history table
- Vignette updated with Drug Card Lookup section and `phase` parameter
- All man pages regenerated via roxygen2

## Migration notes

- **Fully backward compatible** — input CSVs without a `direction` column
  behave exactly as before
- `enrichment = FALSE` (counting mode) still works
- `DGIdbr()` function signature is unchanged — no new parameters to learn
- Direction scoring is **automatic** when a `direction` column is present;
  if ChEMBL is unreachable, direction columns are silently filled with `NA`
- `drug_card(phase = "all")` is the default — all existing `drug_card()` calls
  continue to work unchanged

## Installation

```r
# install.packages("remotes")
remotes::install_github("lfzhang00/DGIdbr")
```

## Quick example

```r
library(DGIdbr)

# Gene-set driven drug prioritisation (direction scoring automatic)
DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv")

# Drug card lookup with phase filter
card <- drug_card("ASPIRIN", phase = "approved")
card$target_genes           # "PTGS2"
card$mechanisms$action_type # "INHIBITOR"
card$indications$disease    # 5 approved indications
```

## Dependencies

- ChEMBL REST API (`https://www.ebi.ac.uk/chembl/api/data`)
- No new R package dependencies (httr, jsonlite, utils — unchanged)
- Behind an HTTP proxy? Set `NO_PROXY=www.ebi.ac.uk,ebi.ac.uk`

## Known limitations

- ChEMBL mechanism coverage is incomplete; drugs without mechanism data get
  `NA` in direction columns (e.g. METFORMIN)
- Some targets are non-protein entities (e.g. cisplatin targets DNA) and
  cannot be resolved to a gene symbol
