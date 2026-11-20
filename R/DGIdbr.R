#' DGIdb GraphQL endpoint
#'
#' Uses env var `DGIDB_URL` if set, otherwise defaults to public endpoint.
endpoint <- "https://dgidb.org/api/graphql"

#' Run a GraphQL POST request and parse JSON
#' @keywords internal
run_gql <- function(query, variables = NULL, timeout_sec = getOption("DGIdbr.timeout", 120)) {
  if (!is.numeric(timeout_sec) || length(timeout_sec) != 1L ||
      !is.finite(timeout_sec) || timeout_sec <= 0) {
    stop("DGIdbr.timeout must be a positive number of seconds.", call. = FALSE)
  }
  transient_errors <- c("curl_error_operation_timedout", "curl_error_couldnt_connect",
                        "curl_error_couldnt_resolve_host", "curl_error_couldnt_resolve_proxy",
                        "curl_error_recv_error", "curl_error_send_error", "curl_error_partial_file",
                        "curl_error_got_nothing", "curl_error_http2", "curl_error_http2_stream",
                        "curl_error_ssl_connect_error")
  for (attempt in seq_len(3L)) {
    resp <- tryCatch(httr::POST(
      url = Sys.getenv("DGIDB_URL", endpoint),
      body = list(query = query, variables = variables), encode = "json",
      httr::timeout(timeout_sec), httr::config(connecttimeout_ms = ceiling(1000 * timeout_sec))
    ), curl_error = function(e) e)
    delay <- 2^(attempt - 1L)
    if (inherits(resp, "error")) {
      if (attempt == 3L || !inherits(resp, transient_errors)) stop(resp)
      reason <- conditionMessage(resp)
    } else {
      status <- httr::status_code(resp)
      if (attempt == 3L || !status %in% c(408L, 429L, 500L, 502L, 503L, 504L)) break
      reason <- paste("HTTP", status)
      header <- httr::headers(resp)[["retry-after"]]
      if (!is.null(header)) {
        seconds <- suppressWarnings(as.numeric(header))
        if (is.na(seconds)) seconds <- as.numeric(httr::parse_http_date(header)) - as.numeric(Sys.time())
        if (is.finite(seconds)) delay <- max(delay, seconds)
      }
    }
    message("DGIdb request failed (", reason, "). Retry ", attempt, "/2 in ", delay, " s.")
    Sys.sleep(delay)
  }

  text_resp <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (status != 200L) stop("DGIdb HTTP ", status, call. = FALSE)
  parsed <- tryCatch(jsonlite::fromJSON(text_resp, simplifyVector = TRUE),
                     error = function(e) stop("Invalid DGIdb JSON: ", conditionMessage(e),
                                              call. = FALSE))
  if (!is.list(parsed) || !any(c("data", "errors") %in% names(parsed))) {
    stop("Invalid DGIdb response: expected data or errors", call. = FALSE)
  }
  if (length(parsed$errors)) {
    stop("DGIdb GraphQL error: ", paste(parsed$errors$message, collapse = "; "),
         call. = FALSE)
  }

  list(
    status = status,
    data   = parsed$data,
    errors = parsed$errors
  )
}

# Collect name-filtered connections without exceeding the requested page size.
dgidb_named_nodes <- function(field, values, selection, timeout_sec = getOption("DGIdbr.timeout", 120)) {
  result <- list()
  for (batch in split(values, ceiling(seq_along(values) / 100L))) {
    after <- "null"
    seen <- character()
    repeat {
      query <- sprintf("{ %s(names: %s, first: 100, after: %s) {
                          nodes { %s } pageInfo { hasNextPage endCursor } } }",
                       field, jsonlite::toJSON(batch), after, selection)
      connection <- run_gql(query, timeout_sec = timeout_sec)$data[[field]]
      nodes <- connection$nodes
      if (!is.list(nodes)) stop("Missing or invalid ", field, ".nodes", call. = FALSE)
      if (is.data.frame(nodes)) {
        nodes <- lapply(seq_len(nrow(nodes)), function(i) lapply(nodes, `[[`, i))
      }
      result <- c(result, nodes)
      if (!isTRUE(connection$pageInfo$hasNextPage)) break
      cursor <- connection$pageInfo$endCursor
      if (length(cursor) != 1L || is.na(cursor) || !nzchar(cursor) || cursor %in% seen) {
        stop("Invalid pagination cursor for DGIdb ", field, call. = FALSE)
      }
      seen <- c(seen, cursor)
      after <- jsonlite::toJSON(cursor, auto_unbox = TRUE)
    }
  }
  result
}

#' Run DGIdb query for one gene set and write CSV outputs
#'
#' @param direction Run ChEMBL direction scoring when direction data are supplied (default TRUE). Set FALSE to skip ChEMBL while preserving gene-set grouping.
#' @param label Label for the gene set (used in logs/output paths)
#' @param genes Character vector of gene symbols
#' @param out_dir Output directory; created if missing
#' @param approve Logical; restrict the hits table to drugs marked approved by DGIdb. The raw table is unfiltered.
#' @param enrichment Logical; if TRUE compute hypergeometric enrichment
#'   p-values and FDR. Requires additional API calls per drug. Default TRUE.
#' @param background_N Integer or NULL; number of background genes for
#'   hypergeometric test. If NULL (default), automatically queries DGIdb
#'   for the total gene record count. Set a number to override.
#' @param gene_direction Optional named character vector mapping gene symbols
#'   to \code{"up"} or \code{"down"}. When provided, ChEMBL is queried for
#'   drug mechanism of action and the consistency between drug action direction
#'   (inhibitor/suppressor vs. activator/enhancer) and gene expression direction
#'   is scored. Added columns: \code{n_with_direction_data},
#'   \code{n_direction_consistent}, \code{direction_ratio} in the output.
#'   Default \code{NULL} (skip direction analysis).
#' @return Invisibly returns NULL; writes `dgidb_hits.csv` and `dgidb_raw.csv`
#' @export
run_gene_set <- function(label, genes, out_dir, approve = TRUE,
                         enrichment = TRUE, background_N = NULL,
                         gene_direction = NULL, direction = TRUE) {
  if (!is.logical(approve) || length(approve) != 1L || is.na(approve)) {
    stop("approve must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(direction) || length(direction) != 1L || is.na(direction)) {
    stop("direction must be TRUE or FALSE.", call. = FALSE)
  }
  genes <- toupper(trimws(genes))
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0) {
    cat(label, ": no genes, skipped.\n\n")
    return(invisible())
  }
  if (direction && length(gene_direction)) {
    direction_names <- toupper(trimws(names(gene_direction)))
    gene_direction <- stats::setNames(tolower(trimws(gene_direction)), direction_names)
  }
  cat("== Running set:", label, "Genes:", length(genes), "==\n")
  nodes <- dgidb_named_nodes("genes", genes, "
      name
      interactions {
        drug { name conceptId approved }
        interactionTypes { type }
        interactionScore
        sources { sourceDbName }
      }")
  cat(label, ": matched gene records:", length(nodes), "\n")

  rows <- list()
  if (!is.null(nodes) && length(nodes) > 0) {
    for (g in nodes) {
      inters <- g$interactions
      if (!is.list(inters)) stop(label, ": invalid gene interactions", call. = FALSE)
      if (is.data.frame(inters)) {
        inters <- lapply(seq_len(nrow(inters)), function(i) {
          list(drug = inters$drug[i, , drop = FALSE],
               interactionScore = inters$interactionScore[[i]],
               interactionTypes = inters$interactionTypes[[i]],
               sources = inters$sources[[i]])
        })
      }
      for (it in inters) {
        types <- unlist(it$interactionTypes, use.names = FALSE)
        sources <- unlist(it$sources, use.names = FALSE)
        sources <- sort(unique(sources[!is.na(sources) & nzchar(sources)]))
        if (isTRUE(approve) && !"approved" %in% names(it$drug)) {
          stop(label, ": missing drug approval field", call. = FALSE)
        }
        approval <- unlist(it$drug$approved, use.names = FALSE)
        if (length(approval) && (!is.logical(approval) || length(approval) != 1L)) {
          stop(label, ": invalid drug approval field", call. = FALSE)
        }
        if (isTRUE(approve) && length(approval) && is.na(approval)) {
          stop(label, ": missing drug approval value", call. = FALSE)
        }
        rows[[length(rows) + 1]] <- data.frame(
          gene = g$name,
          drug = it$drug$name,
          drug_concept_id = it$drug$conceptId,
          interaction_score = if (length(it$interactionScore)) it$interactionScore else NA_real_,
          interaction_types = paste(types[!is.na(types)], collapse = ";"),
          sources = if (length(sources)) paste(sources, collapse = ";") else NA_character_,
          approved = if (length(approval)) approval else NA,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) {
    cat(label, ": No gene-drug interactions returned.\n")
    df <- data.frame(gene = character(), drug = character(), drug_concept_id = character(),
                     interaction_score = numeric(), interaction_types = character(),
                     sources = character(), approved = logical(), key = character())
    agg <- data.frame(drug = character(), drug_concept_id = character(), key = character(),
                      total_score = numeric(), max_score = numeric(), gene_count = integer(),
                      genes = character(), sources = character(), approved = logical())
  } else {
    df <- do.call(rbind, rows)
    df$key <- paste(df$drug, df$drug_concept_id, sep = "|")
    agg <- stats::aggregate(
      interaction_score ~ drug + drug_concept_id + key,
      data = df,
      na.action = stats::na.pass,
      FUN = function(x) {
        scored <- x[!is.na(x)]
        if (!length(scored)) return(c(total = NA_real_, max = NA_real_))
        c(total = sum(scored), max = max(scored))
      }
    )
    scores_mat <- agg$interaction_score
    agg$total_score <- scores_mat[, "total"]
    agg$max_score <- scores_mat[, "max"]
    agg$interaction_score <- NULL

    genes_by_key <- tapply(df$gene, df$key, function(x) sort(unique(x)))
    agg$gene_count <- sapply(agg$key, function(k) length(genes_by_key[[k]]))
    agg$genes <- sapply(agg$key, function(k) paste(genes_by_key[[k]], collapse = ";"))
    collapse_sources <- function(vals) {
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) == 0) return(NA_character_)
      vals <- unlist(strsplit(vals, ";"))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) == 0) return(NA_character_)
      paste(sort(unique(vals)), collapse = ";")
    }
    unique_keys <- unique(df$key)
    sources_by_key <- sapply(unique_keys, function(k) collapse_sources(df$sources[df$key == k]), USE.NAMES = TRUE)
    agg$sources <- sources_by_key[agg$key]
    agg$approved <- df$approved[match(agg$key, df$key)]
  }

  if (isTRUE(approve)) {
    approved_rows <- !is.na(agg$approved) & agg$approved
    agg <- agg[approved_rows, ]
    cat(label, ": Filtered to approved drugs only. Remaining:", nrow(agg), "\n")
  }

  # ---- Enrichment analysis ----
  if (isTRUE(enrichment) && nrow(agg) > 0) {
    cat("\n", label, ": Computing enrichment statistics...\n")
    drugs_to_query <- unique(agg$drug)
    tgt <- fetch_drug_target_counts(drugs_to_query)
    if (!is.null(tgt)) {
      agg <- merge(agg, tgt, by = c("drug", "drug_concept_id"), all.x = TRUE, sort = FALSE)
      sample_n <- if (is.null(background_N)) length(unique(vapply(nodes, `[[`, character(1), "name"))) else length(genes)
      agg <- compute_drug_enrichment(agg, n_genes_input = sample_n, background_N = background_N)
      # Sort by FDR ascending, then enrichment_ratio descending
      agg <- agg[order(agg$fdr, -agg$enrichment_ratio, na.last = TRUE), ]
      cat("Enrichment done. Drugs with FDR < 0.05:", sum(agg$fdr < 0.05, na.rm = TRUE), "\n")
    } else {
      warning("DGIdb target counts unavailable; enrichment statistics set to NA.", call. = FALSE)
      agg$total_targets <- NA_integer_
      agg$enrichment_ratio <- agg$p_value <- agg$fdr <- NA_real_
      agg$significance <- NA_character_
      agg <- agg[order(-agg$gene_count, -agg$total_score, -agg$max_score), ]
    }
  } else {
    agg <- agg[order(-agg$gene_count, -agg$total_score, -agg$max_score), ]
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, "dgidb_hits.csv")
  out_raw <- file.path(out_dir, "dgidb_raw.csv")
  write_result_tables(agg, df, out_file, out_raw)

  if (direction && length(gene_direction) > 0 && nrow(agg) > 0) {
    dir_result <- compute_direction_consistency(agg, df, gene_direction)
    agg <- dir_result$agg
    df <- dir_result$raw
    write_result_tables(agg, df, out_file, out_raw)
  }

  cat("==", label, "aggregated drug hits ==\n")
  cat("Total candidate drugs:", nrow(agg), "\n")
  if (isTRUE(enrichment) && "fdr" %in% names(agg)) {
    cat("Top 10 by FDR / enrichment:\n")
    cols <- intersect(c("drug", "gene_count", "total_targets",
                        "enrichment_ratio", "p_value", "fdr", "significance"),
                      names(agg))
    print(utils::head(agg[, cols], 10))
  } else {
    cat("Top 10 by gene_count/score:\n")
    print(utils::head(agg[, c("drug", "gene_count", "total_score", "max_score")], 10))
  }
  cat("\nFiles written to:\n -", out_file, "\n -", out_raw, "\n\n")
}

# Serialize both tables before replacing either checkpoint. Deferring interrupts
# across the two renames keeps a user interrupt from publishing half an update.
write_result_tables <- function(hits, raw, hits_path, raw_path) {
  paths <- c(hits_path, raw_path)
  staged <- tempfile("dgidbr-", tmpdir = dirname(hits_path), fileext = c("-hits.csv", "-raw.csv"))
  on.exit(unlink(staged), add = TRUE)
  utils::write.csv(hits, staged[1], row.names = FALSE)
  utils::write.csv(raw, staged[2], row.names = FALSE)
  suspendInterrupts({
    for (i in seq_along(paths)) {
      if (!file.rename(staged[i], paths[i])) stop("Could not publish ", paths[i], call. = FALSE)
    }
  })
  invisible(NULL)
}

#' Warn about direction values that cannot be used
#' @keywords internal
warn_unrecognised_direction <- function(direction, where) {
  values <- tolower(trimws(direction))
  known <- values %in% c("up", "down")
  blank <- is.na(values) | !nzchar(values)
  if (any(blank)) {
    warning(sprintf("%s: %d row(s) have an empty direction and are not used.",
                    where, sum(blank)), call. = FALSE)
  }
  unknown <- unique(values[!known & !blank])
  if (length(unknown)) {
    warning(sprintf("%s: direction value(s) %s are not 'up' or 'down'; %d row(s) are not used.",
                    where, paste0("'", sort(unknown), "'", collapse = ", "),
                    sum(!known & !blank)), call. = FALSE)
  }
  invisible(NULL)
}

#' Make a value usable as one path component
#' @keywords internal
safe_path_component <- function(value, what) {
  raw <- tolower(trimws(value))
  safe <- gsub("[^a-z0-9._-]+", "_", raw)
  safe <- gsub("^\\.+", "_", safe)
  if (!nzchar(safe)) safe <- "unknown"
  if (!identical(safe, raw)) {
    warning(sprintf("%s '%s' contains characters that are unsafe in a directory name; using '%s'.",
                    what, value, safe), call. = FALSE)
  }
  safe
}

#' Read genes from CSV, optionally filtered by direction
#' @keywords internal
read_genes_with_dir <- function(file, direction_value = NULL) {
  df <- if (is.data.frame(file)) file else utils::read.csv(file, stringsAsFactors = FALSE)
  if (!"gene" %in% names(df)) stop("Missing gene column")
  if (!is.null(direction_value) && "direction" %in% names(df)) {
    df <- df[tolower(trimws(df$direction)) %in% direction_value, , drop = FALSE]
  }
  genes <- toupper(trimws(df$gene))
  unique(genes[!is.na(genes) & nzchar(genes)])
}

#' Build gene sets for a chosen mode
#'
#' Each gene set is a list with components \code{label} (character),
#' \code{genes} (character vector), \code{gene_direction} (named character
#' vector \code{"up"}/\code{"down"}, or \code{NULL}), and \code{out}
#' (output directory path).
#'
#' @param mode "group" or "subtype"
#' @param base_tables Directory containing input tables
#' @param base_out Output base directory
#' @param group_filename Case-control file name (requires gene; direction is optional)
#' @param subtype_filename Subtype file name (needs gene/direction/subtype)
#' @return List of gene set lists (label, genes, gene_direction, out)
#' @export
build_gene_sets <- function(mode = c("group", "subtype"),
                            base_tables,
                            base_out,
                            group_filename = "group.csv",
                            subtype_filename = "subtype.csv") {
  mode <- match.arg(mode)
  gene_sets <- list()

  if (mode == "group") {
    group_file <- file.path(base_tables, group_filename)
    if (file.exists(group_file)) {
      df_all <- utils::read.csv(group_file, stringsAsFactors = FALSE)
      if ("direction" %in% names(df_all)) {
        warn_unrecognised_direction(df_all$direction, basename(group_file))
        directions <- c("up", "down")
      } else {
        directions <- "all"
      }
      gene_sets <- lapply(directions, function(dir) {
        genes <- read_genes_with_dir(df_all, dir)
        list(label = paste0("group_", dir), genes = genes,
             gene_direction = if (dir != "all") stats::setNames(rep(dir, length(genes)), genes),
             out = file.path(base_out, "dgidb_group", dir))
      })
    } else {
      stop("Group file not found: ", group_file, call. = FALSE)
    }
    return(gene_sets)
  }

  # subtype mode
  subtype_file <- file.path(base_tables, subtype_filename)
  if (!file.exists(subtype_file)) {
    stop("Subtype file not found: ", subtype_file, call. = FALSE)
  }

  df_sub <- utils::read.csv(subtype_file, stringsAsFactors = FALSE)
  needed_cols <- c("gene", "direction", "subtype")
  if (!all(needed_cols %in% names(df_sub))) {
    stop("Subtype file must contain columns: ", paste(needed_cols, collapse = ", "))
  }
  df_sub$direction <- tolower(trimws(df_sub$direction))
  df_sub$subtype <- trimws(as.character(df_sub$subtype))
  warn_unrecognised_direction(df_sub$direction, basename(subtype_file))
  dropped <- is.na(df_sub$subtype) | !nzchar(df_sub$subtype)
  if (any(dropped)) {
    warning(sprintf("%s: %d row(s) have an empty subtype and are not used.",
                    basename(subtype_file), sum(dropped)), call. = FALSE)
  }

  subtypes <- sort(unique(df_sub$subtype[!is.na(df_sub$subtype) & nzchar(df_sub$subtype)]))
  if (length(subtypes) == 0) {
    cat("Warning: no subtype values detected in", subtype_file, "\n")
    return(gene_sets)
  }

  subtype_dirs <- vapply(subtypes, safe_path_component, character(1), what = "subtype")
  names(subtype_dirs) <- subtypes
  shared <- duplicated(subtype_dirs) | duplicated(subtype_dirs, fromLast = TRUE)
  if (any(shared)) {
    stop("Distinct subtype labels map to the same output directory: ",
         paste(shQuote(subtypes[shared]), collapse = ", "), ". Rename the labels.", call. = FALSE)
  }
  for (st in subtypes) {
    sub_df <- df_sub[df_sub$subtype == st, ]
    for (dirv in c("up", "down")) {
      genes <- read_genes_with_dir(sub_df, dirv)
      if (length(genes) == 0) next
      gd_vec <- stats::setNames(rep(dirv, length(genes)), genes)
      gene_sets[[length(gene_sets) + 1]] <- list(
        label = sprintf("subtype_%s_%s", subtype_dirs[[st]], dirv),
        genes = genes,
        gene_direction = gd_vec,
        out   = file.path(base_out, "dgidb_subtype", subtype_dirs[[st]], dirv)
      )
    }
  }

  gene_sets
}

#' Run a list of gene sets
#' @keywords internal
run_gene_sets <- function(gene_sets, approve = TRUE, enrichment = TRUE, background_N = NULL,
                          direction = TRUE) {
  failures <- character()
  for (gs in gene_sets) {
    tryCatch(
      run_gene_set(gs$label, gs$genes, gs$out,
                   approve = approve, enrichment = enrichment,
                   background_N = background_N, gene_direction = gs$gene_direction,
                   direction = direction),
      error = function(e) {
        detail <- paste(gs$label, conditionMessage(e), sep = ": ")
        message(detail)
        failures <<- c(failures, detail)
      }
    )
  }
  if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)
}

#' Main entry point for DGIdbr package
#'
#' @param mode "group" or "subtype"
#' @param base_tables Input directory for tables
#' @param base_out Output base directory
#' @param direction Run ChEMBL direction scoring when direction data are supplied (default TRUE). Set FALSE to skip ChEMBL while preserving gene-set grouping.
#' @param group_filename Group file name (default group.csv)
#' @param subtype_filename Subtype file name (default subtype.csv)
#' @param approve Logical; restrict the hits table to drugs marked approved by DGIdb. The raw table is unfiltered.
#' @param enrichment Logical; if TRUE compute hypergeometric enrichment
#'   p-values and FDR for each drug. Default TRUE.
#' @param background_N Integer; number of background genes for
#'   hypergeometric test (default: current DGIdb gene count).
#' @details
#' \strong{Direction consistency scoring:} When the input CSV contains a
#' \code{direction} column (values \code{"up"} / \code{"down"}), the package
#' automatically queries ChEMBL for drug mechanism of action and scores the
#' consistency between drug action direction (inhibitor/suppressor vs.
#' activator/enhancer) and gene expression direction. See the output columns
#' \code{n_with_direction_data}, \code{n_direction_consistent}, and
#' \code{direction_ratio} in \code{dgidb_hits.csv}.  Requires internet access
#' to the ChEMBL API.  If ChEMBL is unavailable or no mechanism data is found,
#' these columns are filled with \code{NA} (the analysis continues without
#' error).
#' @return Invisibly returns NULL
#' @export
#' @name DGIdbr
#' @examples
#' \dontrun{
#' DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv", base_out = ".")
#' }
DGIdbr <- function(mode = "group",
                   base_tables = file.path("results", "nsNMF", "diff", "Tables"),
                   base_out = file.path("results", "nsNMF", "drug", "Tables"),
                   group_filename = "group.csv",
                   subtype_filename = "subtype.csv",
                   approve = TRUE,
                   enrichment = TRUE,
                   background_N = NULL, direction = TRUE) {

  gs <- build_gene_sets(
    mode = mode,
    base_tables = base_tables,
    base_out = base_out,
    group_filename = group_filename,
    subtype_filename = subtype_filename
  )
  if (length(gs) == 0) {
    cat("No gene sets found in", base_tables, "\n")
    return(invisible(NULL))
  }
  run_gene_sets(gs, approve = approve,
                enrichment = enrichment, background_N = background_N, direction = direction)
  cat("DGIdb run finished.\n")
  invisible(NULL)
}
