# DGIdbr statistical enrichment helpers
#
# These functions add hypergeometric enrichment testing with FDR correction
# to transform the drug-gene lookup into a statistically grounded analysis.

# ---- In-memory cache for drug target counts ----
# Avoid re-querying the same drug across multiple gene sets in one session.
.cache_drug_targets <- new.env(parent = emptyenv())

# ---- Cached druggable gene count (background N) ----
# Stored in an environment so it can be mutated after namespace lock.
# Fetched once from DGIdb per session.
.cache_bg <- new.env(parent = emptyenv())

#' Fetch the total number of DGIdb gene records
#'
#' Queries the DGIdb genes connection for its \code{totalCount}. This is the
#' background size for all registered genes, including those without
#' a documented drug interaction.  Result is cached in memory for the session.
#'
#' @return Integer; number of gene records in DGIdb.
#' @keywords internal
fetch_druggable_gene_count <- function() {
  service <- Sys.getenv("DGIDB_URL", endpoint)
  if (!is.null(.cache_bg[[service]])) {
    return(.cache_bg[[service]])
  }

  q <- "{ genes(first: 1) { totalCount } }"
  resp <- tryCatch(run_gql(q), error = function(e) e)

  count <- if (!inherits(resp, "error") && identical(as.integer(resp$status), 200L) &&
               is.null(resp$errors)) resp$data$genes$totalCount
  if (is.numeric(count) && length(count) == 1L && is.finite(count) &&
      count > 0 && count <= .Machine$integer.max && count == floor(count)) {
    .cache_bg[[service]] <- as.integer(count)
    cat("DGIdb background gene count:", .cache_bg[[service]], "\n")
  } else {
    warning("DGIdb background unavailable; enrichment statistics set to NA.", call. = FALSE)
    return(NA_integer_)
  }

  .cache_bg[[service]]
}

#' Fetch total gene target counts for a list of drugs from DGIdb
#'
#' Queries DGIdb for each drug's full interaction list and counts the number
#' of unique genes that the drug is known to interact with. Results are cached
#' in memory so the same drug is only queried once per session.
#'
#' @param drug_names Character vector of drug names
#' @param chunk Integer; number of drugs per API request (default 150)
#' @return data.frame with columns: drug, drug_concept_id, total_targets
#' @keywords internal
fetch_drug_target_counts <- function(drug_names, chunk = 150) {
  service <- Sys.getenv("DGIDB_URL", endpoint)
  if (is.null(.cache_drug_targets[[service]])) {
    .cache_drug_targets[[service]] <- new.env(parent = emptyenv())
  }
  cache <- .cache_drug_targets[[service]]
  drug_names <- unique(drug_names)
  cached <- intersect(drug_names, ls(cache, all.names = TRUE))
  results <- list()
  for (nm in cached) results <- c(results, cache[[nm]])
  to_query <- setdiff(drug_names, cached)
  for (batch in split(to_query, ceiling(seq_along(to_query) / chunk))) {
    nd <- tryCatch(dgidb_named_nodes("drugs", batch, "name conceptId interactions { gene { name } }"),
                   error = function(e) {
                     warning("Cannot fetch target counts: ", conditionMessage(e), call. = FALSE)
                     list()
                   })
    for (d in nd) {
      inters <- d$interactions
      gene_names <- character(0)

      if (is.data.frame(inters)) {
        # interactions is a data.frame — extract gene$name column
        gene_names <- if ("gene" %in% names(inters)) {
          gcol <- inters$gene
          if (is.data.frame(gcol) && "name" %in% names(gcol)) {
            gcol$name
          } else if (is.list(gcol)) {
            sapply(gcol, function(x) x$name %||% NA_character_)
          } else {
            character(0)
          }
        } else {
          character(0)
        }
      } else if (is.list(inters) && length(inters) > 0) {
        gene_names <- sapply(inters, function(it) {
          if (is.list(it) && !is.null(it$gene)) {
            it$gene$name %||% NA_character_
          } else {
            NA_character_
          }
        })
      }

      gene_names <- unique(gene_names[!is.na(gene_names) & nzchar(gene_names)])
      n_targets <- length(gene_names)

      entry <- data.frame(
        drug            = d$name,
        drug_concept_id = d$conceptId,
        total_targets   = n_targets,
        stringsAsFactors = FALSE
      )

      # Store every node returned for this name: one name can carry several
      # concept IDs (salt forms, source-specific IDs) and the merge downstream
      # is by (drug, drug_concept_id), so keeping only one would silently drop
      # the other pairs in later gene sets of the same session.
      cache[[d$name]] <- c(cache[[d$name]], list(entry))

      results[[length(results) + 1]] <- entry
    }
  }

  if (length(results) == 0) return(NULL)
  out <- do.call(rbind, results)
  out <- out[!duplicated(out[, c("drug", "drug_concept_id")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Compute drug enrichment statistics
#'
#' Adds hypergeometric-test p-values, Benjamini-Hochberg FDR, enrichment
#' ratio, and significance stars to the aggregated drug table.
#'
#' @param agg data.frame from the drug aggregation step; must contain
#'   `gene_count` and `total_targets` columns.
#' @param n_genes_input Integer; number of genes in the input gene set.
#' @param background_N Integer or NULL; total background gene count. If NULL
#'   (default), automatically queries DGIdb for the number of DGIdb gene records.
#'   Set an explicit integer to override.
#' @return The same data.frame with additional columns:
#'   enrichment_ratio, p_value, fdr, significance
#' @keywords internal
compute_drug_enrichment <- function(agg, n_genes_input, background_N = NULL) {
  if (!"total_targets" %in% names(agg)) {
    warning("No 'total_targets' column in aggregated data; skipping enrichment calculation")
    return(agg)
  }

  if (nrow(agg) == 0) return(agg)

  # Resolve background N
  if (is.null(background_N)) {
    background_N <- fetch_druggable_gene_count()
    if (is.na(background_N)) {
      agg$enrichment_ratio <- agg$p_value <- agg$fdr <- NA_real_
      agg$significance <- NA_character_
      return(agg)
    }
  }

  k <- agg$gene_count          # genes from input set that interact with this drug
  m <- agg$total_targets       # total genes this drug targets in DGIdb

  if (anyNA(m) || anyNA(k)) {
    warning(sprintf("%d of %d drugs have no target count; their enrichment statistics are NA.",
                    sum(is.na(m) | is.na(k)), length(m)), call. = FALSE)
  }

  for (size in list(background_N, n_genes_input)) {
    if (length(size) != 1L || !is.numeric(size) || !is.finite(size) ||
        size <= 0 || size != floor(size)) {
      stop("Sample and background sizes must be positive whole numbers.", call. = FALSE)
    }
  }
  if (n_genes_input > background_N) stop("Sample exceeds background size.", call. = FALSE)
  valid <- is.finite(k) & is.finite(m) & k == floor(k) & m == floor(m) &
    k > 0 & k <= pmin(n_genes_input, m) & m <= background_N &
    n_genes_input - k <= background_N - m
  if (any(!valid & !is.na(m))) {
    warning("Incompatible overlap or target counts; affected statistics set to NA.", call. = FALSE)
  }
  agg$enrichment_ratio <- NA_real_
  agg$p_value <- NA_real_
  agg$enrichment_ratio[valid] <- (k[valid] / n_genes_input) / (m[valid] / background_N)
  agg$p_value[valid] <- stats::phyper(k[valid] - 1, m[valid], background_N - m[valid],
                                    n_genes_input, lower.tail = FALSE)

  # Benjamini-Hochberg FDR correction
  agg$fdr <- stats::p.adjust(agg$p_value, method = "BH")

  # Significance stars
  agg$significance <- ifelse(agg$fdr < 0.001, "***",
                      ifelse(agg$fdr < 0.01,  "**",
                      ifelse(agg$fdr < 0.05,  "*",
                      ifelse(agg$fdr < 0.1,   ".",
                      "ns"))))

  agg
}

# Helper: null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x
