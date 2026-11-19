test_that("a timed-out GraphQL read and a 503 are retried with both timeouts set", {
  calls <- 0L
  options_seen <- list()
  local_mocked_bindings(POST = function(url, body, encode, ...) {
    calls <<- calls + 1L
    options_seen <<- lapply(list(...), function(x) x$options)
    if (calls == 1L) {
      stop(structure(list(message = "timed out"), class = c("curl_error_operation_timedout",
                                                             "curl_error", "error", "condition")))
    }
    if (calls == 2L) return(response("unavailable", 503L))
    response('{"data":{"genes":{"nodes":[]}}}')
  }, .package = "httr")
  expect_type(run_gql("{ genes { nodes { name } } }")$data, "list")
  expect_equal(calls, 3L)
  expect_equal(options_seen[[1]]$timeout_ms, 120000)
  expect_equal(options_seen[[2]]$connecttimeout_ms, 120000)
})

test_that("permanent HTTP and GraphQL errors do not trigger retries", {
  reply <- response("bad query", 400L)
  calls <- 0L
  local_mocked_bindings(POST = function(...) { calls <<- calls + 1L; reply }, .package = "httr")
  expect_error(run_gql("query"), "HTTP 400")
  reply <- response('{"errors":[{"message":"Unknown field"}]}')
  expect_error(run_gql("query"), "Unknown field")
  expect_equal(calls, 2L)
})

test_that("direction interruption leaves readable enrichment and successful scoring updates both files", {
  out <- tempfile()
  params <- list(label = "sample", genes = "SIRT5", out_dir = out,
                 background_N = 100, gene_direction = c(SIRT5 = "up"))
  item <- list(drug = list(name = "A", conceptId = "a", approved = TRUE),
               interactionScore = 3, interactionTypes = list(list(type = "inhibitor")),
               sources = list(list(sourceDbName = "test")))
  done <- FALSE
  hits_path <- file.path(out, "dgidb_hits.csv")
  raw_path <- file.path(out, "dgidb_raw.csv")
  local_mocked_bindings(
    dgidb_named_nodes = function(...) list(list(name = "SIRT5", interactions = list(item))),
    fetch_drug_target_counts = function(...) data.frame(drug = "A", drug_concept_id = "a", total_targets = 2),
    compute_direction_consistency = function(agg, df, gene_direction) {
      expect_equal(read.csv(hits_path)$fdr, 0.02)
      expect_identical(read.csv(raw_path)$interaction_types, "inhibitor")
      if (!done) signalCondition(structure(list(message = "cancelled"), class = c("interrupt", "condition")))
      agg$direction_ratio <- 1
      df$direction_match <- "consistent"
      list(agg = agg, raw = df)
    })
  interrupted <- tryCatch({ do.call(run_gene_set, params); FALSE }, interrupt = function(e) TRUE)
  expect_true(interrupted)
  expect_false("direction_ratio" %in% names(read.csv(hits_path)))
  done <- TRUE
  do.call(run_gene_set, params)
  expect_equal(read.csv(hits_path)$direction_ratio, 1)
  expect_identical(read.csv(raw_path)$direction_match, "consistent")
})

test_that("direction FALSE retains separate groups and never calls ChEMBL", {
  dir <- tempfile()
  dir.create(dir)
  write.csv(data.frame(gene = c("SIRT5", "TFAM"), direction = c("up", "down")),
            file.path(dir, "group.csv"), row.names = FALSE)
  item <- list(drug = list(name = "A", conceptId = "a", approved = TRUE), interactionScore = 1)
  local_mocked_bindings(
    dgidb_named_nodes = function(field, values, ...) list(list(name = values, interactions = list(item))),
    compute_direction_consistency = function(...) stop("unexpected ChEMBL request"))
  DGIdbr(base_tables = dir, base_out = dir, direction = FALSE, enrichment = FALSE)
  for (group in c("up", "down")) {
    raw <- read.csv(file.path(dir, "dgidb_group", group, "dgidb_raw.csv"))
    expect_identical(raw$gene, if (group == "up") "SIRT5" else "TFAM")
  }
})
