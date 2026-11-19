response <- function(text, status = 200L) {
  structure(list(status_code = status, content = charToRaw(text),
                 headers = list("content-type" = "application/json"), url = "https://test.invalid"),
            class = "response")
}

test_that("API failures reach the caller and empty hits replace old results", {
  reply <- response("temporarily unavailable", 503L)
  local_mocked_bindings(POST = function(...) reply, .package = "httr")
  out <- tempfile()
  run <- function() run_gene_set("sample", "SIRT5", out, enrichment = FALSE)
  expect_error(run(), "HTTP 503")
  reply <- response("<html/>")
  expect_error(run(), "Invalid DGIdb JSON")
  reply <- response('{"data":null,"errors":[{"message":"unknown field"}]}')
  expect_error(run(), "unknown field")
  reply <- response('{"data":{"genes":{}}}')
  expect_error(run(), "invalid genes.nodes")
  input <- tempfile()
  dir.create(input)
  write.csv(data.frame(gene = "SIRT5", direction = "up"), file.path(input, "group.csv"))
  expect_error(DGIdbr(base_tables = input, base_out = out), "group_up")
  reply <- response('{"data":{"genes":{"nodes":[]}}}')
  dir.create(out)
  files <- file.path(out, c("dgidb_hits.csv", "dgidb_raw.csv"))
  for (f in files) writeLines("stale", f)
  run()
  for (f in files) expect_equal(nrow(read.csv(f)), 0L)
})

test_that("known approval and nested types survive API parsing", {
  drug <- function(id, approval, types) {
    list(drug = list(name = id, conceptId = id, approved = approval),
         interactionTypes = lapply(types, function(x) list(type = x)),
         interactionScore = 1, sources = list(list(sourceDbName = "DTC")))
  }
  interactions <- list(drug("A", TRUE, c("inhibitor", "binder")),
                       drug("B", FALSE, "inhibitor"), drug("C", NULL, character()))
  node <- list(name = "SIRT5", interactions = interactions)
  json <- jsonlite::toJSON(list(data = list(genes = list(nodes = list(node)))),
                           auto_unbox = TRUE, null = "null")
  local_mocked_bindings(POST = function(url, body, ...) {
    expect_match(body$query, "conceptId approved")
    response(as.character(json))
  }, .package = "httr")
  for (approved_only in c(TRUE, FALSE)) {
    out <- tempfile()
    run_gene_set("sample", "SIRT5", out, approve = approved_only, enrichment = FALSE)
    raw <- read.csv(file.path(out, "dgidb_raw.csv"))
    hits <- read.csv(file.path(out, "dgidb_hits.csv"))
    expect_identical(raw$interaction_types[1:2], c("inhibitor;binder", "inhibitor"))
    expect_setequal(hits$drug, if (approved_only) "A" else c("A", "B", "C"))
  }
})

test_that("mechanism pagination retries a failed batch and follows next", {
  urls <- character()
  fail <- TRUE
  local_mocked_bindings(GET = function(url, ...) {
    urls <<- c(urls, url)
    second <- grepl("offset=1", url, fixed = TRUE)
    if (second && fail) return(response("unavailable", 503L))
    row <- list(molecule_chembl_id = "TEST_PAGE", action_type = "INHIBITOR",
                target_chembl_id = if (second) "T2" else "T1")
    more <- if (!second) "/chembl/api/data/mechanism.json?offset=1"
    doc <- list(mechanisms = list(row), page_meta = list("next" = more))
    response(as.character(jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null")))
  }, .package = "httr")
  expect_warning(first <- chembl_fetch_mechanisms("TEST_PAGE"), "503")
  expect_equal(nrow(first), 0L)
  fail <- FALSE
  second <- chembl_fetch_mechanisms("TEST_PAGE")
  expect_identical(second$target_chembl_id, c("T1", "T2"))
  expect_equal(urls[c(2, 4)], rep("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?offset=1", 2))
})
