## data-raw/kegg-refresh-helpers.R -- shared steps of data-raw/kegg.R,
## data-raw/module.R and data-raw/path.R. Sourced by those scripts; not run on
## its own.
##
## The KEGG version is derived by the installed richR's own .kegg_release()
## (the rule richR records for live KEGG labels: the release KEGG states,
## verbatim; otherwise "last update YYYY-MM-DD" from the update date KEGG
## states; otherwise NA). Install richR from this source tree before running a
## refresh so that the rule matches the code being shipped.

kegg_refresh_start <- function() {
  stopifnot(file.exists("DESCRIPTION"),
            requireNamespace("richR", quietly = TRUE),
            exists(".kegg_release", envir = asNamespace("richR"), inherits = FALSE))
  ## Bound every network call: base R connections and the httr requests
  ## KEGGREST makes.
  options(timeout = 60)
  if (requireNamespace("httr", quietly = TRUE)) {
    httr::set_config(httr::timeout(60))
  }
  kegg_release <- get(".kegg_release", envir = asNamespace("richR"))
  list(retrieved = format(Sys.Date()), release = kegg_release())
}

## The object as currently bundled, read before it is overwritten.
kegg_bundled_object <- function(file, name) {
  env <- new.env(parent = emptyenv())
  load(file, envir = env)
  env[[name]]
}

## Refuse a response that looks partial: far fewer entries than bundled.
kegg_check_not_partial <- function(new_n, old_n, what) {
  if (new_n < 0.9 * old_n) {
    stop(what, ": ", new_n, " entries retrieved against ", old_n,
         " bundled; refusing to save what looks like a partial response.", call. = FALSE)
  }
  invisible(TRUE)
}

## Print rows added / removed and labels changed between two named vectors.
kegg_print_diff <- function(old, new, what, examples = 10L) {
  added <- setdiff(names(new), names(old))
  removed <- setdiff(names(old), names(new))
  common <- intersect(names(old), names(new))
  renamed <- common[old[common] != new[common]]
  show <- function(ids) {
    if (length(ids) == 0L) return("")
    paste0(": ", paste(utils::head(ids, examples), collapse = ", "),
           if (length(ids) > examples) ", ..." else "")
  }
  cat(sprintf("== diff %s: %d -> %d rows\n", what, length(old), length(new)))
  cat(sprintf("  added %d%s\n", length(added), show(added)))
  cat(sprintf("  removed %d%s\n", length(removed), show(removed)))
  cat(sprintf("  renamed %d\n", length(renamed)))
  for (id in utils::head(renamed, examples)) {
    cat(sprintf("    %s: \"%s\" -> \"%s\"\n", id, old[[id]], new[[id]]))
  }
  invisible(list(added = added, removed = removed, renamed = renamed))
}

## Stamp the retrieval date, KEGG version and row count of `object` into
## R/data.R. The provenance line of each object has the form
##   Retrieval date (<object>): <date>; KEGG version: <version>; rebuilt by <script>.
## and the row count is the nearest preceding "@format A data.frame with N rows".
kegg_stamp_documentation <- function(object, info, n_rows, script) {
  doc_file <- file.path("R", "data.R")
  doc <- readLines(doc_file)
  line <- grep(paste0("^#'\\s+Retrieval date \\(", gsub(".", "\\.", object, fixed = TRUE),
                      "\\):"), doc)
  stopifnot(length(line) == 1L)
  doc[line] <- sprintf("#'   Retrieval date (%s): %s; KEGG version: %s; rebuilt by %s.",
                       object, info$retrieved,
                       if (is.na(info$release)) "NA" else info$release, script)
  format_lines <- grep("^#' @format A data\\.frame with [0-9,]+ rows", doc)
  format_line <- max(format_lines[format_lines < line])
  doc[format_line] <- sub("with [0-9,]+ rows", sprintf("with %d rows", n_rows),
                          doc[format_line])
  writeLines(doc, doc_file)
  invisible(TRUE)
}

kegg_source_record <- function(info, url) {
  list(retrieved = info$retrieved, release = info$release, url = url)
}
