#' Per-identifier audit of a conversion or an enrichment analysis
#'
#' Given a vector returned by \code{\link{mapGenes}()}, \code{mappingAudit()}
#' returns its decision table as a data frame: the columns \code{position},
#' \code{input_id}, \code{candidates} (a list column), \code{target},
#' \code{decision}, \code{reason} and, for a ranking, \code{score}, one row per
#' input entry in input order.
#'
#' @param x a vector returned by \code{mapGenes()}.
#' @return a data frame, or \code{NULL} with a warning when no audit is
#'   available.
#' @seealso \code{\link{mapGenes}()}
#' @export
mappingAudit <- function(x) {
  record <- if (is.atomic(x)) attr(x, "richR_mapping", exact = TRUE)
  if (!is.null(record)) return(.mapping_decision_frame(record))
  .audit_unavailable("the object carries no mapping record")
}

## The record's decision table (a plain column list) as a data frame with the
## same columns in the same order, candidates as a list column, automatic row
## names.
.mapping_decision_frame <- function(record) {
  list2DF(record$decisions)
}

.audit_unavailable <- function(why) {
  warning("mappingAudit(): mapping audit unavailable: ", why, ".", call. = FALSE)
  NULL
}
