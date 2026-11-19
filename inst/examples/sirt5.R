# Run with Rscript after installing DGIdbr. Output goes to a new directory.
out <- tempfile("sirt5_", tmpdir = getwd())
DGIdbr::run_gene_set("SIRT5", "SIRT5", out, approve = FALSE, enrichment = FALSE)
raw <- read.csv(file.path(out, "dgidb_raw.csv"))
hits <- read.csv(file.path(out, "dgidb_hits.csv"))
stopifnot(all(raw$gene == "SIRT5"), setequal(raw$drug_concept_id, hits$drug_concept_id))
print(raw[c("drug", "drug_concept_id", "interaction_types", "sources")], row.names = FALSE)
