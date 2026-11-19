# DGIdbr 1.3.3

## Behaviour changes

* Enrichment statistics no longer turn into `NA` when several gene sets run in
  one session. Drug target counts were cached per drug name, but one name can
  map to several concept IDs (salt forms, source-specific identifiers); the
  cached path kept a single ID and every later set merged by
  `(drug, drug_concept_id)` against an incomplete table. In the bundled
  `group.csv` example this left `total_targets`, `p_value` and `fdr` as `NA` for
  86 of 2078 down-regulated hits when the up set ran first. Cached counts are
  now stored per returned concept ID.
* An unusable `direction` value (anything other than `up`/`down`, including an
  empty cell) now raises a warning that names the value and counts the unused
  rows. Before, those rows were dropped silently and the set could run with zero
  genes and no message other than `no genes, skipped.`.
* A subtype value is sanitised before it becomes an output directory (for
  example `..` or `a/b`), with a warning when the name changes; such values can
  no longer write outside `dgidb_subtype/`. Rows with an empty subtype are
  reported instead of being dropped silently.
* `approve = TRUE` errors when DGIdb returns a null approval value instead of
  silently treating the drug as unapproved.
* Drugs whose target count is missing are reported with a warning; their
  enrichment statistics stay `NA` instead of disappearing from the message log.
* `drug_card()` errors when the ChEMBL lookup itself fails (HTTP error or
  network error) instead of reporting the drug as "not found"; a failed
  indication request is reported as unavailable rather than as "no indication
  data".

## Fixes

* Duplicate `(drug, drug_concept_id)` rows returned by the target lookup are
  dropped before the merge, so a drug that DGIdb reports twice under one name no
  longer appears twice in `dgidb_hits.csv`.
* `drug_card()` accepts `CHEMBL:999` and `CHEMBL:CHEMBL999` as well as
  `CHEMBL999`.
* README and README_zh describe the real `dgidb_hits.csv` column set (9 base
  columns, 4 enrichment columns when `enrichment = TRUE`, 3 direction columns
  when scoring runs) and note that `target_genes` is `NA` for complex targets.

# DGIdbr 1.3.2

## Behaviour changes

* `DGIdbr()` and `run_gene_set()` gain a `direction` argument (default `TRUE`).
  `direction = FALSE` skips the ChEMBL mechanism-of-action scoring while keeping
  the up/down grouping, so a gene-set query no longer depends on ChEMBL being
  reachable.
* A failed primary query now raises an R error instead of writing empty tables.
  `DGIdbr()` finishes the remaining gene sets and then reports every failed set.
  Header-only `dgidb_hits.csv` / `dgidb_raw.csv` now mean a successful query that
  returned no interactions.

## Fixes

* HTTP, GraphQL, DNS and invalid-JSON failures are propagated with their status
  and message; previously a server error and a genuine zero-hit result were
  indistinguishable to the caller.
* `approve = TRUE` keeps only drugs that DGIdb marks `approved = TRUE`, and errors
  when the approval field is missing or not a boolean. Before, a failed approval
  lookup silently produced an unfiltered hits table.
* `interaction_types` (e.g. `inhibitor`) are read from the nested
  `interactionTypes` selection and kept in `dgidb_raw.csv`.
* DGIdb name lookups are paged (batches of 100 names, `pageInfo` cursors) and
  transient failures are retried up to twice. The per-request timeout defaults to
  120 s and follows `options(DGIdbr.timeout = )`; the previous fixed 20 s timeout
  with no retry aborted large gene sets.
* ChEMBL mechanism paging follows the API's `next` link instead of rebuilding
  `page=`; a batch that fails mid-way is retried on a later call and never cached
  as complete.
* `dgidb_hits.csv` and `dgidb_raw.csv` are written before ChEMBL scoring and
  rewritten afterwards, so a slow or failing ChEMBL stage cannot discard the
  DGIdb results.
* Drug hits with a missing interaction score are retained with `NA` scores
  instead of being dropped by the aggregation.
* Enrichment validates its inputs. When the background is looked up
  automatically, the sample is the number of matched DGIdb gene records rather
  than the raw input length; an explicit `background_N` uses all distinct input
  symbols. Incompatible overlap or target counts are reported and set to `NA`
  instead of yielding invalid p-values.

## Packaging and documentation

* `License: MIT + file LICENSE` with a valid `LICENSE` plus `LICENSE.md`;
  `R CMD check` reports 0 errors, 0 warnings and 0 notes.
* Added testthat coverage for API failures, approval parsing, pagination,
  aggregation and enrichment.
* README and README_zh explain how to read 1.3.2 output (error vs. zero hits,
  raw vs. filtered tables), the `direction` switch and the timeout option, and
  include a SIRT5 example (`inst/examples/sirt5.R`).
