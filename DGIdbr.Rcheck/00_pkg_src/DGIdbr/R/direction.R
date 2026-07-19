# DGIdbr -- Drug mechanism-of-action direction analysis
#
# Queries the ChEMBL API to determine drug action direction
# (inhibitor/suppressor vs. activator/enhancer per target gene),
# then scores consistency with input gene expression direction.
#
# =========================================================================
# ChEMBL API reference: https://chembl.gitbook.io/chembl-interface-documentation
# Molecule endpoint:    /chembl/api/data/molecule (use pref_name__in)
# Mechanism endpoint:   /chembl/api/data/mechanism
# Target endpoint:      /chembl/api/data/target/{id}.json
# IMPORTANT: The /drug.json endpoint has broken name filtering -- always use
#            /molecule.json with pref_name__in for name lookups!
# =========================================================================

CHEMBL_API_BASE <- "https://www.ebi.ac.uk/chembl/api/data"

# ---- Session caches ----
# Avoid re-querying the same drug or target across calls.
.cache_mol_id     <- new.env(parent = emptyenv())   # drug name -> molecule_chembl_id
.cache_mechanisms <- new.env(parent = emptyenv())   # molecule_chembl_id -> mechanism data.frame
.cache_target_gene <- new.env(parent = emptyenv())  # target_chembl_id -> gene symbol

# =========================================================================
# Step 1 -- Drug name -> ChEMBL molecule ID
# =========================================================================

#' Map drug names to ChEMBL molecule IDs via the molecule endpoint
#'
#' Queries \verb{GET /molecule.json?pref_name__in=...} for exact name match,
#' then falls back to \verb{molecule_synonyms} search.  The drug endpoint
#' (\verb{/drug.json}) should NOT be used -- its \code{name} field is often
#' \code{None} and the filter returns incorrect results.  Results are cached
#' in the session environment \code{.cache_mol_id}.
#'
#' @param drug_names Character vector of drug names.
#' @param chunk      Batch size for the \code{__in} filter (default 20).
#' @return Named list: drug name -> molecule_chembl_id (character), or
#'   \code{NULL} if no match found.
#' @keywords internal
chembl_map_drug_to_molecule <- function(drug_names, chunk = 20) {
  drug_names <- unique(drug_names[!is.na(drug_names) & nzchar(drug_names)])
  if (length(drug_names) == 0) return(stats::setNames(list(), character(0)))

  # --- cache read ---
  out <- stats::setNames(vector("list", length(drug_names)), drug_names)
  to_query <- character(0)
  for (nm in drug_names) {
    if (exists(nm, envir = .cache_mol_id, inherits = FALSE)) {
      out[[nm]] <- .cache_mol_id[[nm]]
    } else {
      to_query <- c(to_query, nm)
    }
  }
  if (length(to_query) == 0) return(out)

  # --- batch query via molecule endpoint ---
  # Use pref_name__in for exact name matching; molecule endpoint has correct
  # data for pref_name and molecule_synonyms, unlike the drug endpoint.
  for (i in seq(1, length(to_query), by = chunk)) {
    batch <- to_query[i:min(i + chunk - 1, length(to_query))]
    q <- paste(vapply(batch, utils::URLencode, character(1)), collapse = ",")
    url  <- paste0(CHEMBL_API_BASE, "/molecule.json?pref_name__in=", q,
                   "&limit=100")

    doc <- tryCatch(
      jsonlite::fromJSON(
        httr::content(httr::GET(url, httr::timeout(20)),
                       "text", encoding = "UTF-8"),
        simplifyVector = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(doc) || is.null(doc$molecules) ||
        !is.data.frame(doc$molecules) || nrow(doc$molecules) == 0) {
      for (nm in batch) { .cache_mol_id[[nm]] <- NULL; out[[nm]] <- NULL }
      next
    }

    molecules <- doc$molecules
    # Build lookup: lowercase pref_name -> molecule_chembl_id
    name_map <- stats::setNames(molecules$molecule_chembl_id,
                                tolower(molecules$pref_name))

    for (nm in batch) {
      key <- tolower(nm)
      if (key %in% names(name_map)) {
        val <- name_map[[key]]
        .cache_mol_id[[nm]] <- val
        out[[nm]] <- val
        next
      }
      # Fallback: search molecule_synonyms
      found <- FALSE
      for (j in seq_len(nrow(molecules))) {
        syns <- molecules$molecule_synonyms[[j]]
        if (is.data.frame(syns) && "synonym" %in% names(syns) &&
            any(key %in% tolower(syns$synonym))) {
          .cache_mol_id[[nm]] <- molecules$molecule_chembl_id[j]
          out[[nm]] <- molecules$molecule_chembl_id[j]
          found <- TRUE
          break
        }
      }
      if (!found) { .cache_mol_id[[nm]] <- NULL; out[[nm]] <- NULL }
    }
  }

  out
}

# =========================================================================
# Step 2 -- Fetch mechanisms for ChEMBL molecule IDs
# =========================================================================

#' Fetch mechanisms of action for ChEMBL molecule IDs
#'
#' Queries \verb{GET /mechanism.json?molecule_chembl_id__in=...} and collects
#' the (target, action_type) pairs.  Handles pagination transparently.
#' Results are cached per molecule ID.
#'
#' @param mol_ids Character vector of ChEMBL molecule IDs (e.g. CHEMBL11359).
#' @param chunk   Number of molecule IDs per API request (default 50).
#' @return data.frame with columns:
#'   \code{molecule_chembl_id}, \code{target_chembl_id},
#'   \code{action_type}, \code{mechanism_of_action}.
#'   Empty data.frame (0 rows) if none found or API unavailable.
#' @keywords internal
chembl_fetch_mechanisms <- function(mol_ids, chunk = 50) {
  mol_ids <- unique(mol_ids[!is.na(mol_ids) & nzchar(mol_ids)])
  if (length(mol_ids) == 0) return(data.frame())

  # --- cache read ---
  cached_list <- list()
  to_query <- character(0)
  for (mid in mol_ids) {
    if (exists(mid, envir = .cache_mechanisms, inherits = FALSE)) {
      cached_list[[mid]] <- .cache_mechanisms[[mid]]
    } else {
      to_query <- c(to_query, mid)
    }
  }
  if (length(to_query) == 0) {
    return(do.call(rbind, cached_list))
  }

  # --- batch query ---
  for (i in seq(1, length(to_query), by = chunk)) {
    batch  <- to_query[i:min(i + chunk - 1, length(to_query))]
    q      <- paste(batch, collapse = ",")
    accum  <- list()
    page   <- 1
    max_pages <- 5  # safety valve

    repeat {
      url <- paste0(CHEMBL_API_BASE, "/mechanism.json",
                    "?molecule_chembl_id__in=", q,
                    "&limit=1000&page=", page)

      doc <- tryCatch(
        jsonlite::fromJSON(
          httr::content(httr::GET(url, httr::timeout(30)),
                         "text", encoding = "UTF-8"),
          simplifyVector = TRUE
        ),
        error = function(e) NULL
      )

      if (is.null(doc) || is.null(doc$mechanisms)) break
      mechs <- doc$mechanisms
      if (!is.data.frame(mechs) || nrow(mechs) == 0) break

      accum[[length(accum) + 1]] <- mechs

      # Pagination: check if there's a "next" page in page_meta
      # (ChEMBL uses `page_meta`, NOT `page_metadata`)
      meta <- doc$page_meta
      if (!is.list(meta) || is.null(meta[["next"]])) {
        # No "next" URL -> last page (or single-page result)
        break
      }
      if (page >= max_pages) break
      page <- page + 1
    }

    if (length(accum) == 0) {
      for (mid in batch) .cache_mechanisms[[mid]] <- data.frame()
      next
    }

    all_mechs <- do.call(rbind, accum)

    # Subset to relevant columns (NOTE: mechanism endpoint does NOT
    # return target_name -- only target_chembl_id is available)
    needed <- c("molecule_chembl_id", "target_chembl_id",
                "action_type", "mechanism_of_action")
    present <- intersect(needed, names(all_mechs))
    if (length(present) < 2) {
      for (mid in batch) .cache_mechanisms[[mid]] <- data.frame()
      next
    }
    all_mechs <- all_mechs[, present, drop = FALSE]

    # Cache per molecule
    for (mid in batch) {
      sub <- all_mechs[all_mechs$molecule_chembl_id == mid, , drop = FALSE]
      .cache_mechanisms[[mid]] <- sub
      cached_list[[mid]] <- sub
    }
  }

  do.call(rbind, cached_list)
}

# =========================================================================
# Step 3 -- Resolve target_chembl_id to human gene symbol
# =========================================================================

#' Resolve ChEMML target IDs to gene symbols
#'
#' For one or more ChEMBL target IDs, queries the target endpoint and
#' extracts the \code{GENE_SYMBOL} synonym from its components.
#' Results are cached.
#'
#' @param target_ids Character vector of ChEMBL target IDs (e.g. CHEMBL255).
#' @return Named character vector: target_id -> gene symbol (or \code{NA}).
#' @keywords internal
chembl_target_to_gene <- function(target_ids) {
  target_ids <- unique(target_ids[!is.na(target_ids) & nzchar(target_ids)])
  if (length(target_ids) == 0) return(stats::setNames(character(0), character(0)))

  out <- stats::setNames(rep(NA_character_, length(target_ids)), target_ids)
  to_query <- character(0)
  for (tid in target_ids) {
    if (exists(tid, envir = .cache_target_gene, inherits = FALSE)) {
      out[tid] <- .cache_target_gene[[tid]]
    } else {
      to_query <- c(to_query, tid)
    }
  }
  if (length(to_query) == 0) return(out)

  for (tid in to_query) {
    url <- paste0(CHEMBL_API_BASE, "/target/", tid, ".json")
    doc <- tryCatch(
      jsonlite::fromJSON(
        httr::content(httr::GET(url, httr::timeout(15)),
                       "text", encoding = "UTF-8"),
        simplifyVector = TRUE
      ),
      error = function(e) NULL
    )

    symbol <- NA_character_
    if (!is.null(doc)) {
      comps <- doc$target_components
      if (is.data.frame(comps) && nrow(comps) > 0 &&
          "target_component_synonyms" %in% names(comps)) {
        for (k in seq_len(nrow(comps))) {
          syns <- comps$target_component_synonyms[[k]]
          if (is.data.frame(syns) && "component_synonym" %in% names(syns) &&
              "syn_type" %in% names(syns)) {
            gs <- syns$component_synonym[syns$syn_type == "GENE_SYMBOL"]
            if (length(gs) > 0 && nzchar(gs[1])) {
              symbol <- gs[1]
              break
            }
          }
        }
      }
    }

    .cache_target_gene[[tid]] <- symbol
    out[tid] <- symbol
  }

  out
}

# =========================================================================
# Step 4 -- Classify action type into direction category
# =========================================================================

#' Classify ChEMBL action types into suppress/enhance/unknown
#'
#' Maps raw action types (e.g. INHIBITOR, AGONIST) to a high-level
#' direction category for consistency scoring.
#'
#' @param x Character vector of action types (or \code{NULL}/\code{NA}).
#' @return Character vector of the same length, with values
#'   \code{"suppress"}, \code{"enhance"}, or \code{"unknown"}.
#' @keywords internal
classify_action_direction <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))

  # Normalise
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA")] <- NA_character_

  # Suppress / inhibit categories
  suppress <- c(
    "INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR",
    "NEGATIVE_ALLOSTERIC_MODULATOR", "INVERSE AGONIST", "INVERSE_AGONIST",
    "PARTIAL ANTAGONIST", "PARTIAL_ANTAGONIST", "ALLOSTERIC ANTAGONIST",
    "ALLOSTERIC_ANTAGONIST", "GATING INHIBITOR", "GATING_INHIBITOR",
    "SUPPRESSOR", "SUPPRESSOR", "DOWNREGULATOR", "DOWN-REGULATOR",
    "INHIBITORY ALLOSTERIC MODULATOR", "INHIBITORY_ALLOSTERIC_MODULATOR",
    "DECREASED EXPRESSION", "DECREASED_EXPRESSION",
    "NEGATIVE MODULATOR", "NEGATIVE_MODULATOR",
    "ANTISENSE INHIBITOR", "ANTISENSE_INHIBITOR"
  )

  # Enhance / activate categories
  enhance <- c(
    "ACTIVATOR", "AGONIST", "PARTIAL AGONIST", "PARTIAL_AGONIST",
    "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE_ALLOSTERIC_MODULATOR",
    "ALLOSTERIC AGONIST", "ALLOSTERIC_AGONIST",
    "STIMULATOR", "STIMULATOR",
    "POSITIVE MODULATOR", "POSITIVE_MODULATOR",
    "ACTIVATING ALLOSTERIC MODULATOR", "ACTIVATING_ALLOSTERIC_MODULATOR",
    "UPREGULATOR", "UP-REGULATOR",
    "INCREASED EXPRESSION", "INCREASED_EXPRESSION",
    "RELEASER", "CHAPERONE",
    "COACTIVATOR", "CO-ACTIVATOR",
    "POTENTIATOR"
  )

  result <- rep("unknown", length(x))
  result[is.na(x)] <- "unknown"
  result[x %in% suppress] <- "suppress"
  result[x %in% enhance]  <- "enhance"
  result
}

# =========================================================================
# Step 5 -- Compute direction consistency scores
# =========================================================================

#' Compute drug-target direction consistency scores
#'
#' For each drug in an aggregated DGIdb result, queries ChEMBL for the
#' mechanism of action, resolves targets to gene symbols, and compares
#' the drug's action direction (suppress/enhance) with the gene's
#' expression direction (up/down).  The following rule is applied:
#'
#' \itemize{
#'   \item UP gene + suppress action = \strong{consistent}
#'         (drug reduces overactivity)
#'   \item DOWN gene + enhance action = \strong{consistent}
#'         (drug boosts underactivity)
#'   \item UP gene + enhance action = \strong{inconsistent}
#'         (drug could further activate an already overactive gene)
#'   \item DOWN gene + suppress action = \strong{inconsistent}
#'         (drug could further suppress an already underactive gene)
#'   \item unknown action type = \strong{unknown}
#' }
#'
#' @param agg_df       data.frame; the aggregated drug table from
#'                     \code{run_gene_set} (must contain columns
#'                     \code{drug} and \code{genes}).
#' @param raw_df       data.frame; the raw gene-drug interaction table
#'                     (must contain \code{gene}, \code{drug}).
#' @param gene_direction Named character vector mapping gene symbols to
#'                     \code{"up"} or \code{"down"}.  Typically passed
#'                     from the input CSV.
#' @return A list with two elements:
#'   \describe{
#'     \item{agg}{The \code{agg_df} with added columns:
#'       \code{n_with_direction_data}, \code{n_direction_consistent},
#'       \code{direction_ratio}.}
#'     \item{raw}{The \code{raw_df} with added column
#'       \code{direction_match} (\code{"consistent"},
#'       \code{"inconsistent"}, or \code{NA}).}
#'   }
#'   If ChEMBL is unreachable or no direction data is available, the
#'   added columns are filled with \code{NA} and a message is printed.
#' @keywords internal
compute_direction_consistency <- function(agg_df, raw_df, gene_direction) {
  # Safety: skip if no gene direction info
  if (is.null(gene_direction) || length(gene_direction) == 0) {
    return(.add_empty_direction_cols(agg_df, raw_df))
  }

  # Get unique drug names
  drugs <- unique(agg_df$drug)
  drugs <- drugs[!is.na(drugs) & nzchar(drugs)]
  if (length(drugs) == 0) {
    return(.add_empty_direction_cols(agg_df, raw_df))
  }

  cat("Checking drug mechanism direction (ChEMBL)...\n")

  # ---- Step A: drug name -> molecule ID ----
  mol_map <- chembl_map_drug_to_molecule(drugs)
  found <- !vapply(mol_map, is.null, logical(1))
  if (!any(found)) {
    cat("  No ChEMBL molecule IDs found for these drugs; skipping direction.\n")
    return(.add_empty_direction_cols(agg_df, raw_df))
  }

  # ---- Step B: fetch mechanisms ----
  mol_ids <- unique(as.character(mol_map[found]))
  mechs <- chembl_fetch_mechanisms(mol_ids)
  if (!is.data.frame(mechs) || nrow(mechs) == 0) {
    cat("  No mechanism data returned from ChEMBL; skipping direction.\n")
    return(.add_empty_direction_cols(agg_df, raw_df))
  }

  # Ensure required columns exist (NOTE: mechanism endpoint does NOT return
  # target_name -- only target_chembl_id is available)
  for (col in c("molecule_chembl_id", "target_chembl_id", "action_type")) {
    if (!col %in% names(mechs)) {
      cat("  ChEMBL response missing column '", col, "'; skipping direction.\n", sep = "")
      return(.add_empty_direction_cols(agg_df, raw_df))
    }
  }

  # ---- Step C: classify action types ----
  mechs$direction_cat <- classify_action_direction(mechs$action_type)

  # ---- Step D: resolve targets to gene symbols ----
  # Collect all unique target_chembl_ids and resolve to gene symbols via API
  all_target_ids <- unique(
    mechs$target_chembl_id[!is.na(mechs$target_chembl_id) & nzchar(mechs$target_chembl_id)]
  )
  target_gene_map <- chembl_target_to_gene(all_target_ids)

  # Build reverse map: molecule -> list of drug names
  mol_to_drugs <- list()
  for (nm in names(mol_map)) {
    mid <- mol_map[[nm]]
    if (!is.null(mid)) {
      mol_to_drugs[[mid]] <- c(mol_to_drugs[[mid]], nm)
    }
  }

  # ---- Step E: build drug -> gene -> action lookup ----
  # drug_target_action[[drug_name]][[gene_symbol]] = "suppress" / "enhance"
  drug_target_action <- list()

  for (k in seq_len(nrow(mechs))) {
    mid  <- mechs$molecule_chembl_id[k]
    tid  <- mechs$target_chembl_id[k]
    dcat <- mechs$direction_cat[k]

    if (dcat == "unknown" || is.na(tid) || !nzchar(tid)) next

    # Resolve target_chembl_id -> gene symbol
    gene <- target_gene_map[tid]
    if (is.na(gene) || !nzchar(gene)) next
    if (!(gene %in% names(gene_direction))) next  # not in input genes

    # Find drug names for this molecule
    drug_names <- mol_to_drugs[[mid]]
    if (is.null(drug_names)) next

    for (dn in drug_names) {
      if (is.null(drug_target_action[[dn]])) drug_target_action[[dn]] <- list()
      drug_target_action[[dn]][[gene]] <- dcat
    }
  }

  # ---- Step F: score each (drug, gene) pair ----
  raw_df$direction_match <- NA_character_

  for (k in seq_len(nrow(raw_df))) {
    drug_name <- raw_df$drug[k]
    gene_name <- raw_df$gene[k]

    if (is.na(drug_name) || is.na(gene_name)) next
    if (is.null(drug_target_action[[drug_name]])) next

    found_action <- drug_target_action[[drug_name]][[gene_name]]
    if (is.null(found_action) || is.na(found_action)) next

    gene_dir <- gene_direction[gene_name]
    if (is.na(gene_dir)) next

    if ((gene_dir == "up" && found_action == "suppress") ||
        (gene_dir == "down" && found_action == "enhance")) {
      raw_df$direction_match[k] <- "consistent"
    } else {
      raw_df$direction_match[k] <- "inconsistent"
    }
  }

  # ---- Step G: aggregate direction stats to drug level ----
  if (!"direction_match" %in% names(raw_df) ||
      all(is.na(raw_df$direction_match))) {
    cat("  No direction matches could be resolved.\n")
    return(.add_empty_direction_cols(agg_df, raw_df))
  }

  direction_stats <- stats::aggregate(
    direction_match ~ drug,
    data = raw_df[!is.na(raw_df$direction_match), , drop = FALSE],
    FUN = function(x) {
      c(
        n_with_data   = sum(!is.na(x)),
        n_consistent  = sum(x == "consistent", na.rm = TRUE)
      )
    },
    simplify = TRUE
  )

  if (nrow(direction_stats) == 0) {
    return(.add_empty_direction_cols(agg_df, raw_df))
  }

  mat <- direction_stats$direction_match
  direction_stats$n_with_direction_data  <- mat[, "n_with_data"]
  direction_stats$n_direction_consistent  <- mat[, "n_consistent"]
  direction_stats$direction_match <- NULL

  direction_stats$direction_ratio <-
    direction_stats$n_direction_consistent /
    pmax(direction_stats$n_with_direction_data, 1)

  # Merge into agg_df
  agg_df <- merge(agg_df, direction_stats,
                  by = "drug", all.x = TRUE, sort = FALSE)

  n_with_data <- sum(agg_df$n_with_direction_data > 0, na.rm = TRUE)
  cat("  Drug direction data:",
      n_with_data, "drugs with mechanism info,",
      sum(agg_df$n_direction_consistent > 0, na.rm = TRUE),
      "with consistent pairs.\n")

  list(agg = agg_df, raw = raw_df)
}

# =========================================================================
# Helper: add empty direction columns
# =========================================================================

#' Add NA-filled direction columns
#'
#' Used when ChEMBL is unreachable or no direction data is available,
#' so the output CSV schema stays consistent.
#'
#' @inheritParams compute_direction_consistency
#' @return List with \code{agg} and \code{raw} elements, each with
#'   direction columns filled with \code{NA}.
#' @keywords internal
.add_empty_direction_cols <- function(agg_df, raw_df) {
  agg_df$n_with_direction_data  <- NA_integer_
  agg_df$n_direction_consistent  <- NA_integer_
  agg_df$direction_ratio        <- NA_real_

  if (!"direction_match" %in% names(raw_df)) {
    raw_df$direction_match <- NA_character_
  }

  list(agg = agg_df, raw = raw_df)
}

# =========================================================================
# Convenience: get gene direction from a gene set label
# =========================================================================

# =========================================================================
# Drug lookup -- search a drug, get targets + indications
# =========================================================================

#' Fetch drug indications from ChEMBL
#'
#' Queries \verb{GET /drug_indication.json?molecule_chembl_id=...} and returns
#' the list of indications with clinical phases.
#'
#' @param mol_id ChEMBL molecule ID (e.g. \code{"CHEMBL25"}).
#' @param max_results Maximum number of indications to return (default 50).
#' @return data.frame with columns \code{disease}, \code{mesh_heading},
#'   \code{efo_id}, \code{max_phase_for_ind}, or an empty data.frame if none.
#' @keywords internal
chembl_fetch_indications <- function(mol_id, max_results = 50) {
  if (is.null(mol_id) || is.na(mol_id) || !nzchar(mol_id)) {
    return(data.frame())
  }

  url <- paste0(CHEMBL_API_BASE, "/drug_indication.json",
                "?molecule_chembl_id=", utils::URLencode(mol_id, reserved = TRUE),
                "&limit=", max_results)

  doc <- tryCatch(
    jsonlite::fromJSON(
      httr::content(httr::GET(url, httr::timeout(20)),
                     "text", encoding = "UTF-8"),
      simplifyVector = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(doc) || is.null(doc$drug_indications) ||
      !is.data.frame(doc$drug_indications) || nrow(doc$drug_indications) == 0) {
    return(data.frame())
  }

  ind <- doc$drug_indications
  needed <- c("mesh_heading", "efo_term", "efo_id", "max_phase_for_ind")
  present <- intersect(needed, names(ind))
  if (length(present) == 0) return(data.frame())

  out <- ind[, present, drop = FALSE]

  # Normalise: prefer mesh_heading, fall back to efo_term
  if ("mesh_heading" %in% names(out)) {
    out$disease <- out$mesh_heading
  } else if ("efo_term" %in% names(out)) {
    out$disease <- out$efo_term
  } else {
    out$disease <- rep(NA_character_, nrow(out))
  }

  # Round phase for display
  if ("max_phase_for_ind" %in% names(out)) {
    phase <- as.numeric(out$max_phase_for_ind)
    out$phase_label <- ifelse(is.na(phase), "?",
                       ifelse(phase >= 4, "Approved",
                       ifelse(phase == 0, "Preclinical",
                              sprintf("Phase %g", phase))))
  }

  out[!duplicated(out$disease), , drop = FALSE]
}


#' Look up a drug: targets, mechanism, indications
#'
#' Queries ChEMBL for a drug's molecule information, mechanism of action
#' (with gene targets), and clinical indications.  Results are both printed
#' as a human-readable summary and returned invisibly as a structured list.
#'
#' @param drug_name Character; drug name (e.g. \code{"ASPIRIN"},
#'   \code{"CISPLATIN"}, \code{"TAMOXIFEN"}).
#' @param phase Filter for indications by development phase:
#'   \code{"all"} (default, show everything), \code{"approved"} (Phase 4 /
#'   approved only), or \code{"trial"} (Phase 1-3, investigational only).
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{drug}{Drug name.}
#'     \item{molecule_chembl_id}{ChEMBL molecule ID.}
#'     \item{max_phase}{Highest clinical phase reached.}
#'     \item{mechanisms}{data.frame of mechanisms (target, action_type).}
#'     \item{target_genes}{Named character vector: target_chembl_id -> gene symbol.}
#'     \item{indications}{data.frame of indications with phases.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' drug_card("ASPIRIN")
#' drug_card("CISPLATIN")
#' }
drug_card <- function(drug_name, phase = c("all", "approved", "trial")) {
  phase <- match.arg(phase)
  if (is.null(drug_name) || is.na(drug_name) || !nzchar(drug_name)) {
    stop("Please provide a drug name.")
  }

  cat("\n== Drug Card:", drug_name, "==\n\n")

  # ---- Step 1: drug name -> molecule ID ----
  cat("=> Looking up in ChEMBL...\n")
  drug_name_orig <- drug_name
  drug_name <- toupper(trimws(drug_name))
  mol_map <- chembl_map_drug_to_molecule(drug_name)
  mol_id <- mol_map[[drug_name]]

  if (is.null(mol_id)) {
    cat("  x Drug '", drug_name_orig, "' not found in ChEMBL.\n", sep = "")
    return(invisible(NULL))
  }

  # ---- Step 2: molecule info ----
  url <- paste0(CHEMBL_API_BASE, "/molecule.json?pref_name__exact=",
                utils::URLencode(drug_name, reserved = TRUE), "&limit=1")
  doc <- tryCatch(
    jsonlite::fromJSON(
      httr::content(httr::GET(url, httr::timeout(15)),
                     "text", encoding = "UTF-8"),
      simplifyVector = TRUE
    ),
    error = function(e) NULL
  )

  max_phase <- NA_integer_
  if (!is.null(doc) && !is.null(doc$molecules) && nrow(doc$molecules) > 0) {
    max_phase <- doc$molecules$max_phase[1]
  }

  cat("  + ChEMBL ID:", mol_id, "\n")
  if (!is.na(max_phase)) {
    phase_txt <- if (max_phase >= 4) "Approved"
                 else if (max_phase >= 3) sprintf("Phase %g", max_phase)
                 else if (max_phase >= 1) sprintf("Phase %g (investigational)", max_phase)
                 else "Preclinical"
    cat("  + Status:", phase_txt, "\n")
  }
  cat("\n")

  # ---- Step 3: mechanisms ----
  cat("=> Targets & Mechanism of Action:\n")
  mechs <- chembl_fetch_mechanisms(mol_id)

  if (is.data.frame(mechs) && nrow(mechs) > 0 &&
      "target_chembl_id" %in% names(mechs)) {
    # Resolve targets to gene symbols
    target_ids <- unique(mechs$target_chembl_id[
      !is.na(mechs$target_chembl_id) & nzchar(mechs$target_chembl_id)
    ])
    target_genes <- chembl_target_to_gene(target_ids)

    for (k in seq_len(nrow(mechs))) {
      tid <- mechs$target_chembl_id[k]
      gene <- if (!is.null(tid) && !is.na(tid) && tid %in% names(target_genes)) {
        target_genes[tid]
      } else NA_character_
      gene_str <- if (!is.na(gene) && nzchar(gene)) sprintf(" (%s)", gene) else ""
      moa <- mechs$mechanism_of_action[k]
      moa_str <- if (!is.na(moa) && nzchar(moa)) sprintf(" -- %s", moa) else ""
      action_str <- mechs$action_type[k]
      cat(sprintf("  * %s%s%s%s\n",
                  if (!is.na(action_str) && nzchar(action_str)) action_str else "?",
                  gene_str,
                  if (nchar(moa_str) > 0) moa_str else "",
                  if (!is.na(tid)) sprintf(" [%s]", tid) else ""))
    }
    cat("\n")
  } else {
    cat("  (No mechanism data in ChEMBL)\n\n")
    target_genes <- stats::setNames(character(0), character(0))
  }

  # ---- Step 4: indications ----
  cat("=> Indications")
  if (phase == "approved") cat(" (approved only)")
  if (phase == "trial") cat(" (investigational only)")
  cat(":\n")
  indications <- chembl_fetch_indications(mol_id)

  if (is.data.frame(indications) && nrow(indications) > 0 &&
      "disease" %in% names(indications)) {
    # Filter by phase
    if (phase != "all" && "max_phase_for_ind" %in% names(indications)) {
      phase_n <- as.numeric(indications$max_phase_for_ind)
      before <- nrow(indications)
      if (phase == "approved") {
        indications <- indications[!is.na(phase_n) & phase_n >= 4, ]
      } else { # trial
        indications <- indications[!is.na(phase_n) & phase_n >= 1 & phase_n < 4, ]
      }
      n_removed <- before - nrow(indications)
      if (n_removed > 0) cat(sprintf("  (%d indications filtered out)\n", n_removed))
    }
    # Sort by phase descending
    if ("max_phase_for_ind" %in% names(indications)) {
      phase_numeric <- as.numeric(indications$max_phase_for_ind)
      indications <- indications[order(-phase_numeric,
                                        tolower(indications$disease)), ]
    }
    for (k in seq_len(min(nrow(indications), 15))) {
      row <- indications[k, ]
      phase_txt <- if (!is.null(row$phase_label) && !is.na(row$phase_label)) {
        row$phase_label
      } else if (!is.null(row$max_phase_for_ind) && !is.na(row$max_phase_for_ind)) {
        p <- as.numeric(row$max_phase_for_ind)
        ifelse(is.na(p), "?",
               ifelse(p >= 4, "Approved",
                      sprintf("Phase %g", p)))
      } else "?"
      cat(sprintf("  * %s  (%s)\n", row$disease, phase_txt))
    }
    if (nrow(indications) > 15) {
      cat(sprintf("  ... and %d more\n", nrow(indications) - 15))
    }
    cat("\n")
  } else {
    cat("  (No indication data in ChEMBL)\n\n")
  }

  # ---- Return structured result ----
  result <- list(
    drug               = drug_name,
    molecule_chembl_id = mol_id,
    max_phase          = max_phase,
    mechanisms         = if (is.data.frame(mechs) && nrow(mechs) > 0) mechs else data.frame(),
    target_genes       = target_genes,
    indications        = if (is.data.frame(indications) && nrow(indications) > 0) indications else data.frame()
  )

  invisible(result)
}
