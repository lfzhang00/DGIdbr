test_that("CSV inputs normalize symbols and do not duplicate undirected gene sets", {
  d <- tempfile()
  dir.create(d)
  f <- file.path(d, "group.csv")
  write.csv(data.frame(gene = c(" sirt5", "SIRT5", "TFAM ")), f)
  sets <- build_gene_sets("group", d, d)
  expect_length(sets, 1L)
  expect_identical(sets[[1]]$genes, c("SIRT5", "TFAM"))
  write.csv(data.frame(gene = c(" sirt5 ", "tfam"), direction = c(" Up ", "DOWN")), f)
  sets <- build_gene_sets("group", d, d)
  expect_identical(sets[[1]]$gene_direction, c(SIRT5 = "up"))
  expect_identical(sets[[2]]$genes, "TFAM")
  write.csv(data.frame(gene=" tfam ", direction=" Down ", subtype=" B "), file.path(d,"subtype.csv"))
  expect_identical(build_gene_sets("subtype", d, d)[[1]]$gene_direction, c(TFAM = "down"))
})

test_that("impossible counts cannot yield significant enrichment", {
  agg <- data.frame(gene_count = c(1, 3, 1), total_targets = c(4, 2, 0))
  expect_warning(out <- compute_drug_enrichment(agg, 3, 100), "counts")
  expect_equal(out$p_value[1], stats::phyper(0, 4, 96, 3, lower.tail = FALSE))
  expect_true(all(is.na(out[2:3, c("enrichment_ratio", "p_value", "fdr")])))
  expect_error(compute_drug_enrichment(agg, 3, 2), "background")
})

test_that("known ChEMBL identifiers bypass name lookup", {
  local_mocked_bindings(GET = function(...) stop("unexpected name lookup"), .package = "httr")
  known <- c("CHEMBL:CHEMBL477936", "chembl25")
  expect_equal(unname(unlist(chembl_map_drug_to_molecule(known))), c("CHEMBL477936","CHEMBL25"))
})

test_that("conflicting mechanisms remain unresolved in either response order", {
  records <- data.frame(molecule_chembl_id = c("M1", "M1", "M2"), target_chembl_id = "T",
                        action_type = c("AGONIST", "INHIBITOR", "INHIBITOR"))
  local_mocked_bindings(chembl_map_drug_to_molecule = function(...) list(A="M1", B="M2"),
                        chembl_fetch_mechanisms = function(...) records,
                        chembl_target_to_gene = function(...) c(T="SIRT5"))
  agg <- data.frame(drug=c("A", "B"))
  raw <- data.frame(drug=c("A", "B"), gene="SIRT5")
  for (i in 1:2) {
    out <- compute_direction_consistency(agg, raw, c(SIRT5="up"))
    expect_identical(out$raw$direction_match, c(NA_character_, "consistent"))
    expect_identical(out$agg$drug, agg$drug)
    records <- records[rev(seq_len(nrow(records))), ]
  }
})

test_that("synonyms work without preferred-name hits and transient errors can be retried", {
  attempts <- 0L
  local_mocked_bindings(chembl_get_json = function(url, query) {
    if (!is.null(query$pref_name__in)) {
      attempts <<- attempts + 1L
      return(if (attempts > 1L) list(molecules=list()))
    }
    ids <- if (query$molecule_synonyms__molecule_synonym__iexact == "TEST_ALIAS") "CHEMBL123" else c("X","Y")
    list(molecules=data.frame(molecule_chembl_id=ids))
  })
  expect_null(chembl_map_drug_to_molecule("TEST_ALIAS")[[1]])
  resolved <- chembl_map_drug_to_molecule(c("TEST_ALIAS","TEST_AMBIGUOUS"))
  expect_identical(resolved$TEST_ALIAS, "CHEMBL123")
  expect_null(resolved$TEST_AMBIGUOUS)
})

test_that("human single-component targets can be retried without mislabelling complexes", {
  comps <- data.frame(accession = "ABC")
  comps$target_component_synonyms <- list(data.frame(component_synonym="TFAM",syn_type="GENE_SYMBOL"))
  first <- TRUE
  local_mocked_bindings(chembl_get_json = function(url, ...) {
    if (first) return(NULL)
    list(organism=if (grepl("MOUSE",url)) "Mus musculus" else "Homo sapiens",
         target_components=if (grepl("COMPLEX",url)) rbind(comps, comps) else comps)
  })
  expect_true(is.na(chembl_target_to_gene("RETRY")))
  first <- FALSE
  expect_identical(unname(chembl_target_to_gene(c("RETRY","COMPLEX","MOUSE"))),
                   c("TFAM",NA_character_,NA_character_))
})

test_that("multi-gene enrichment uses the matched population and preserves drug rank", {
  dir <- tempfile()
  dir.create(dir)
  write.csv(data.frame(gene=c("TFAM", " sirt5 ", "Unmatched")), file.path(dir,"group.csv"))
  endpoint_key <- paste0("https://", basename(dir), ".invalid")
  previous_endpoint <- Sys.getenv("DGIDB_URL", unset=NA_character_)
  on.exit(if (is.na(previous_endpoint)) Sys.unsetenv("DGIDB_URL") else Sys.setenv(DGIDB_URL=previous_endpoint))
  Sys.setenv(DGIDB_URL=endpoint_key)
  item <- list(drug=list(name="TEST_DRUG", conceptId="test:1", approved=TRUE),
               interactionScore=1, sources=list(), interactionTypes=list())
  gene <- function(name) list(name=name, interactions=list(item))
  local_mocked_bindings(POST = function(url, body, ...) {
    q <- body$query
    if (grepl("totalCount", q, fixed=TRUE)) {
      data <- list(genes=list(totalCount=100))
    } else if (grepl("drugs(names:", q, fixed=TRUE)) {
      targets <- lapply(c("SIRT5","TFAM","X"), function(x) list(gene=list(name=x)))
      data <- list(drugs=list(nodes=list(list(name="TEST_DRUG", conceptId="test:1", interactions=targets))))
    } else {
      last <- grepl('after: "next"', q, fixed=TRUE)
      data <- list(genes=list(nodes=list(gene(if (last) "TFAM" else "SIRT5")),
                              pageInfo=list(hasNextPage=!last, endCursor="next")))
    }
    response(as.character(jsonlite::toJSON(list(data=data), auto_unbox=TRUE)))
  }, .package="httr")
  DGIdbr(base_tables=dir, base_out=dir)
  hits <- read.csv(file.path(dir,"dgidb_group","all","dgidb_hits.csv"))
  expect_equal(hits$gene_count, 2L)
  expect_equal(hits$total_targets, 3L)
  expect_equal(hits$p_value, choose(3,2) / choose(100,2))
})

test_that("a repeated DGIdb cursor is rejected", {
  local_mocked_bindings(run_gql=function(...) list(data=list(
    genes=list(nodes=list(),pageInfo=list(hasNextPage=TRUE,endCursor="old")))))
  expect_error(dgidb_named_nodes("genes","TEST","name",20),"cursor")
})

test_that("background failure preserves hits and can recover in the same session", {
  calls <- 0L
  local_mocked_bindings(run_gql=function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("unavailable")
    list(status=200L, data=list(genes=list(totalCount=101L)))
  }, .cache_bg=new.env(parent=emptyenv()))
  agg <- data.frame(drug="test", gene_count=1L, total_targets=4L)
  expect_warning(first <- compute_drug_enrichment(agg, 2L), "background unavailable")
  expect_identical(first$drug, "test")
  expect_true(all(is.na(first[c("p_value","fdr","enrichment_ratio")])))
  second <- compute_drug_enrichment(agg, 2L)
  expect_equal(second$p_value, stats::phyper(0,4,97,2,lower.tail=FALSE))
})
