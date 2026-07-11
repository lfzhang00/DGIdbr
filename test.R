# Local usage example for DGIdbr.
#
# Features:
#   - enrichment = TRUE  --> hypergeometric test + FDR correction (default)
#   - enrichment = FALSE --> original behaviour (gene_count ranking only)
#   - background_N       --> set background gene count (default 20000)

# 1) Install from local source (uncomment if needed)
remotes::install_local(".", force = TRUE)

# 2) Load the package
library(DGIdbr)

# 3) Group mode, enrichment ON (default) — keep all drugs
DGIdbr(
  mode = "group",
  base_tables = ".",
  group_filename = "group.csv",
  base_out = "dgidb_group_out",
  approve = FALSE,
  enrichment = TRUE,
  background_N = NULL   # NULL = auto-detect from DGIdb (~11665)
)

# 4) Subtype mode, enrichment ON, only FDA approved drugs
DGIdbr(
  mode = "subtype",
  base_tables = ".",
  subtype_filename = "subtype.csv",
  base_out = "dgidb_subtype_out",
  approve = TRUE,
  enrichment = TRUE,
  background_N = NULL   # NULL = auto-detect from DGIdb (~11665)
)

# 5) Backward-compatible: enrichment OFF (original behaviour)
# DGIdbr(
#   mode = "group",
#   base_tables = ".",
#   group_filename = "group.csv",
#   base_out = "dgidb_group_out",
#   approve = FALSE,
#   enrichment = FALSE
# )
