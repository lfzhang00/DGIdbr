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

#' Fetch the total number of druggable genes in DGIdb
#'
#' Queries the DGIdb genes connection for its \code{totalCount}. This is the
#' correct background N for the hypergeometric test — genes that have at least
#' one documented drug interaction.  Result is cached in memory for the session.
#'
#' @return Integer; number of unique genes with drug interactions in DGIdb.
#' @keywords internal
fetch_druggable_gene_count <- function() {
  if (!is.null(.cache_bg$count)) {
    return(.cache_bg$count)
  }

  q <- "{ genes(first: 1) { totalCount } }"
  resp <- tryCatch(run_gql(q, timeout_sec = 15), error = function(e) e)

  if (!inherits(resp, "error") &&
      resp$status == 200 &&
      is.null(resp$errors) &&
      !is.null(resp$data$genes$totalCount)) {
    .cache_bg$count <- as.integer(resp$data$genes$totalCount)
    cat("DGIdb background gene count:", .cache_bg$count, "\n")
  } else {
    # Fallback: reasonable estimate based on DGIdb v5 publications
    .cache_bg$count <- 12000L
    warning("Could not query DGIdb gene count; using fallback N = ",
            .cache_bg$count)
  }

  .cache_bg$count
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
  # Check cache first
  cached <- mget(drug_names, envir = .cache_drug_targets, ifnotfound = NA)
  cached_hits <- !is.na(cached)
  if (all(cached_hits)) {
    return(do.call(rbind, cached))
  }

  # Only query drugs not in cache
  to_query <- drug_names[!cached_hits]
  results <- list()

  # Add cached results
  for (nm in drug_names[cached_hits]) {
    results[[length(results) + 1]] <- cached[[nm]]
  }

  n <- length(to_query)
  if (n == 0) {
    return(do.call(rbind, results))
  }

  for (i in seq(1, n, by = chunk)) {
    sub <- to_query[i:min(i + chunk - 1, n)]
    sub_arg <- paste(sprintf('"%s"', sub), collapse = ", ")

    q <- sprintf(
      '{
        drugs(names: [%s]) {
          nodes {
            name
            conceptId
            interactions {
              gene { name }
            }
          }
        }
      }',
      sub_arg
    )

    resp <- tryCatch(run_gql(q, timeout_sec = 60), error = function(e) e)

    if (inherits(resp, "error")) {
      warning("fetch_drug_target_counts: request error for chunk starting at index ", i, ": ",
              conditionMessage(resp))
      next
    }
    if (resp$status != 200) {
      warning("fetch_drug_target_counts: HTTP ", resp$status, " for chunk starting at index ", i)
      next
    }
    if (!is.null(resp$errors)) {
      warning("fetch_drug_target_counts: GraphQL errors for chunk starting at index ", i)
      next
    }

    nd <- resp$data$drugs$nodes
    if (is.null(nd)) next

    # Normalise nodes to a list-of-lists
    if (is.data.frame(nd)) {
      nd <- lapply(seq_len(nrow(nd)), function(j) {
        list(
          name       = nd$name[[j]],
          conceptId  = nd$conceptId[[j]],
          interactions = nd$interactions[[j]]
        )
      })
    }

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

      # Store in cache
      .cache_drug_targets[[d$name]] <- entry

      results[[length(results) + 1]] <- entry
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
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
#'   (default), automatically queries DGIdb for the number of druggable genes.
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
  }

  k <- agg$gene_count          # genes from input set that interact with this drug
  m <- agg$total_targets       # total genes this drug targets in DGIdb

  # Enrichment ratio: (k/n) / (m/N)
  #  > 1 = drug's targets are over-represented in the input gene set
  agg$enrichment_ratio <- (k / n_genes_input) / (m / background_N)

  # Hypergeometric test: P(X >= k)
  #   phyper(k-1, m, N-m, n, lower.tail = FALSE)
  agg$p_value <- mapply(
    function(k_i, m_i) {
      if (is.na(m_i) || m_i <= 0 || k_i <= 0)
        return(NA_real_)
      stats::phyper(k_i - 1, m_i, background_N - m_i, n_genes_input,
                    lower.tail = FALSE)
    },
    k_i = k, m_i = m
  )

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
