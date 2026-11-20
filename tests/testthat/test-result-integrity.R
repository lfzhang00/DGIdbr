test_that("a missing input file fails the batch before any output", {
  root <- tempfile()
  for (mode in c("group", "subtype")) {
    expect_error(DGIdbr(mode, base_tables = root, base_out = root), "file not found")
  }
  expect_false(dir.exists(root))
})

test_that("distinct subtype labels cannot use the same output directory", {
  root <- tempfile()
  dir.create(root)
  for (labels in list(c("A", "a"), c("A B", "A_B"))) {
    input <- data.frame(gene = c("SIRT5", "TFAM"), direction = "up", subtype = labels)
    write.csv(input, file.path(root, "subtype.csv"), row.names = FALSE)
    expect_error(suppressWarnings(build_gene_sets("subtype", root, root)), "same output")
  }
  expect_false(dir.exists(file.path(root, "dgidb_subtype")))
})

test_that("unavailable target counts keep a readable unscored hits table", {
  item <- list(drug = list(name = "A", conceptId = "a", approved = TRUE), interactionScore = 1)
  local_mocked_bindings(
    dgidb_named_nodes = function(...) list(list(name = "SIRT5", interactions = list(item))),
    fetch_drug_target_counts = function(...) NULL,
    fetch_druggable_gene_count = function(...) stop("unexpected background query"))
  out <- tempfile()
  expect_warning(run_gene_set("sample", "SIRT5", out), "target counts unavailable")
  hits <- read.csv(file.path(out, "dgidb_hits.csv"))
  expect_identical(hits$drug, "A")
  fields <- c("total_targets", "enrichment_ratio", "p_value", "fdr", "significance")
  expect_true(all(fields %in% names(hits)))
  expect_true(all(is.na(hits[fields])))
})

test_that("invalid background counts are not rounded or cached", {
  count <- 3.5
  local_mocked_bindings(.cache_bg = new.env(parent = emptyenv()),
    run_gql = function(...) list(status = 200L, data = list(genes = list(totalCount = count))))
  for (count in list(3.5, -1, "20", c(12, 13))) {
    expect_warning(result <- fetch_druggable_gene_count(), "background unavailable")
    expect_identical(result, NA_integer_)
  }
  count <- 120L
  expect_identical(fetch_druggable_gene_count(), 120L)
})
