#' Write an enrichment result to a folder of tables
#'
#' Stores a \code{richResult} or \code{GSEAResult} as a folder holding
#' \code{manifest.json} and one UTF-8 tab-separated table per data frame of the
#' object, so that \code{\link{readRich}} can rebuild it exactly.
#'
#' The folder holds \code{result.tsv} (\code{result(x)}), \code{tested.tsv}
#' (the tested table, when one is recorded; see \code{\link{testedTerms}}),
#' \code{detail.tsv} (\code{richResult}) or \code{input.tsv} (\code{GSEAResult}:
#' the ranked input as gene and score, in the order of the \code{input} slot).
#' Each table has a header row and one line per row. In character, factor and
#' list cells a backslash is written \code{\\\\}, TAB \code{\\t}, newline
#' \code{\\n} and carriage return \code{\\r}; a missing value of any class is
#' written \code{\\N}; doubles are written with 17 significant digits
#' (\code{NaN}, \code{Inf}, \code{-Inf} spelled out), which reads back
#' bit-identical; a list cell (a list of character vectors, such as a GSEA
#' \code{leadingEdge} with \code{table = FALSE}) is written as a JSON array. A
#' table with zero columns is recorded in the manifest only, without a file.
#'
#' \code{manifest.json} records the format name \code{richR-result}, the format
#' version \code{1}, the class, every other slot, the
#' \code{\link{enrichmentInfo}} metadata (or its absence), whether a tested
#' table is recorded, and for each table its file name, column names, column
#' classes, row names and row count.
#'
#' Supported column classes are character, integer, double, logical, factor and
#' lists of character vectors. Every check runs before any file is written:
#' \code{writeRich()} stops on an unsupported column, on a table that is not a
#' plain \code{data.frame} (for example the \code{grouped_df} that
#' \code{group_by()} leaves in \code{@result}; convert it with
#' \code{x@result <- as.data.frame(x@result)}), on an \code{@result} without
#' its identifier column (\code{Annot} or \code{pathway}) or \code{Term}, and
#' on an existing non-empty \code{path}. With \code{overwrite = TRUE} only a
#' folder previously written by \code{writeRich()} is replaced; its table files
#' and manifest are removed, nothing else.
#'
#' @param x a \code{richResult} or \code{GSEAResult} object
#' @param path folder to write; created when it does not exist
#' @param overwrite replace a folder previously written by \code{writeRich()}
#' @return \code{path}, invisibly
#' @seealso \code{\link{readRich}}, \code{\link{testedTerms}}
#' @importFrom jsonlite toJSON fromJSON serializeJSON unserializeJSON
#' @export
writeRich <- function(x, path, overwrite = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("writeRich(): `path` must be a single folder path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("writeRich(): `overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  ## Everything is checked and encoded in memory before the first file
  ## operation, so a refused call leaves the disk as it was.
  files <- .rich_io_encode(x)
  replace <- .rich_io_check_target(path, overwrite)
  if (replace) {
    old <- file.path(path, .rich_io_known_files)
    unlink(old[file.exists(old)])
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  ## The manifest goes last: a folder without one is never read as a result.
  for (f in names(files)) .rich_io_write_lines(files[[f]], file.path(path, f))
  invisible(path)
}

## The file names writeRich() writes; the only files overwrite = TRUE removes.
.rich_io_known_files <- c("manifest.json", "result.tsv", "tested.tsv",
                          "detail.tsv", "input.tsv")

.rich_io_format <- "richR-result"
.rich_io_format_version <- 1L

## Non-table slots of a result class (stored typed in the manifest).
.rich_io_table_slots <- c("result", "detail", "tested", "input")
.rich_io_scalar_slots <- function(cls) {
  setdiff(methods::slotNames(cls), .rich_io_table_slots)
}

.rich_io_id_column <- function(cls) {
  if (identical(cls, "GSEAResult")) "pathway" else "Annot"
}

## Validate x and return the folder content as a named list of text lines,
## manifest.json last.
.rich_io_encode <- function(x) {
  cls <- class(x)[1L]
  if (!cls %in% c("richResult", "GSEAResult")) {
    stop("writeRich() supports richResult or GSEAResult objects, not ",
         cls, ".", call. = FALSE)
  }
  tables <- list(result = x@result)
  if (.has_tested(x)) tables$tested <- x@tested
  if (identical(cls, "richResult")) {
    tables$detail <- x@detail
  } else {
    tables$input <- data.frame(
      gene = if (is.null(names(x@input))) rep(NA_character_, length(x@input))
             else names(x@input),
      score = unname(x@input), stringsAsFactors = FALSE)
  }
  for (s in names(tables)) {
    if (!identical(class(tables[[s]]), "data.frame")) {
      stop(sprintf(paste("writeRich(): @%s is a %s, not a plain data.frame;",
                         "convert it first, e.g. x@%s <- as.data.frame(x@%s)."),
                   s, paste(class(tables[[s]]), collapse = "/"), s, s),
           call. = FALSE)
    }
  }
  needed <- c(.rich_io_id_column(cls), "Term")
  lacking <- setdiff(needed, names(x@result))
  if (length(lacking)) {
    stop(sprintf("writeRich(): @result lacks the column(s) %s, which a %s needs.",
                 paste0("'", lacking, "'", collapse = ", "), cls), call. = FALSE)
  }
  ## Column classes of every table, checked before anything is encoded.
  classes <- lapply(names(tables), function(s) {
    tab <- tables[[s]]
    vapply(seq_along(tab), function(j) {
      k <- .rich_io_column_class(tab[[j]])
      if (is.na(k)) {
        stop(sprintf(paste("writeRich(): column '%s' of @%s has class %s;",
                           "supported: character, integer, double, logical,",
                           "factor and lists of character vectors."),
                     names(tab)[j], s, paste(class(tab[[j]]), collapse = "/")),
             call. = FALSE)
      }
      k
    }, character(1))
  })
  names(classes) <- names(tables)

  files <- list()
  entries <- list()
  for (s in names(tables)) {
    tab <- tables[[s]]
    file <- if (ncol(tab)) paste0(s, ".tsv") else NULL
    entry <- list(slot = s, file = file, nrow = nrow(tab), ncol = ncol(tab),
                  columns = .rich_io_column_entries(tab, classes[[s]]),
                  row_names = .rich_io_row_names(tab))
    if (identical(s, "input")) entry$names <- !is.null(names(x@input))
    entries[[length(entries) + 1L]] <- entry
    if (!is.null(file)) files[[file]] <- .rich_io_table_lines(tab, classes[[s]])
  }
  metadata <- attr(x, "richR_metadata", exact = TRUE)
  slots <- lapply(.rich_io_scalar_slots(cls), function(s) methods::slot(x, s))
  names(slots) <- .rich_io_scalar_slots(cls)
  manifest <- list(
    format = .rich_io_format,
    format_version = .rich_io_format_version,
    class = cls,
    tested_recorded = .has_tested(x),
    metadata_recorded = !is.null(metadata),
    metadata = if (is.null(metadata)) NULL else .rich_io_typed_json(metadata),
    slots = .rich_io_typed_json(slots),
    tables = entries
  )
  files[["manifest.json"]] <- as.character(jsonlite::toJSON(
    manifest, auto_unbox = TRUE, null = "null", digits = I(17), pretty = TRUE,
    json_verbatim = TRUE))
  files
}

## The class name the manifest records for a column, or NA when unsupported.
.rich_io_column_class <- function(v) {
  cls <- class(v)
  if (identical(cls, "factor")) {
    return(if (anyNA(levels(v))) NA_character_ else "factor")
  }
  if (identical(cls, "list")) {
    ok <- all(vapply(v, function(e) is.character(e) && is.null(attributes(e)),
                     logical(1)))
    return(if (ok) "list" else NA_character_)
  }
  switch(paste(cls, collapse = "/"),
         character = "character", integer = "integer", numeric = "double",
         logical = "logical", NA_character_)
}

.rich_io_column_entries <- function(tab, classes) {
  lapply(seq_along(tab), function(j) {
    entry <- list(name = names(tab)[j], class = classes[[j]])
    if (identical(classes[[j]], "factor")) entry$levels <- levels(tab[[j]])
    entry
  })
}

## The internal row.names attribute: automatic (compact), integer or character.
.rich_io_row_names <- function(tab) {
  rn <- .row_names_info(tab, type = 0L)
  if (is.character(rn)) return(list(type = "character", values = rn))
  if (length(rn) == 2L && is.na(rn[1L])) {
    return(list(type = "automatic", n = rn[2L]))
  }
  list(type = "integer", values = rn)
}

## Typed JSON (jsonlite::serializeJSON, 17 significant digits): types, NULLs,
## empty vectors and names survive exactly (plan.md B7). Embedded verbatim in
## the manifest.
.rich_io_typed_json <- function(value) {
  structure(jsonlite::serializeJSON(value, digits = I(17)), class = "json")
}

## Escape a character vector for a TSV cell: backslash, TAB, newline and CR
## become \\, \t, \n, \r; a missing value becomes the reserved token \N.
.rich_io_escape <- function(v) {
  out <- enc2utf8(as.character(v))
  hit <- which(!is.na(out) & grepl("[\\\\\t\n\r]", out, perl = TRUE))
  if (length(hit)) {
    y <- out[hit]
    y <- gsub("\\", "\\\\", y, fixed = TRUE)
    y <- gsub("\t", "\\t", y, fixed = TRUE)
    y <- gsub("\n", "\\n", y, fixed = TRUE)
    out[hit] <- gsub("\r", "\\r", y, fixed = TRUE)
  }
  out[is.na(v)] <- "\\N"
  out
}

.rich_io_cells <- function(v, class) {
  switch(class,
    character = .rich_io_escape(v),
    factor = .rich_io_escape(as.character(v)),
    integer = ifelse(is.na(v), "\\N", as.character(v)),
    logical = ifelse(is.na(v), "\\N", ifelse(v, "TRUE", "FALSE")),
    double = {
      out <- sprintf("%.17g", v)
      out[is.na(v) & !is.nan(v)] <- "\\N"
      out
    },
    list = .rich_io_escape(vapply(v, function(e) {
      as.character(jsonlite::toJSON(e, na = "null"))
    }, character(1)))
  )
}

.rich_io_table_lines <- function(tab, classes) {
  header <- paste(.rich_io_escape(names(tab)), collapse = "\t")
  if (!nrow(tab)) return(header)
  cells <- lapply(seq_along(tab), function(j) .rich_io_cells(tab[[j]], classes[[j]]))
  c(header, do.call(paste, c(unname(cells), sep = "\t")))
}

.rich_io_write_lines <- function(lines, file) {
  con <- file(file, open = "wb")
  on.exit(close(con))
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
}

## Decide what to do with `path`: FALSE = write into a new or empty folder,
## TRUE = replace a previous writeRich() folder; anything else stops.
.rich_io_check_target <- function(path, overwrite) {
  if (!file.exists(path)) return(FALSE)
  if (!dir.exists(path)) {
    stop(sprintf("writeRich(): '%s' exists and is not a folder.", path),
         call. = FALSE)
  }
  if (!length(list.files(path, all.files = TRUE, no.. = TRUE))) return(FALSE)
  if (!overwrite) {
    stop(sprintf(paste("writeRich(): folder '%s' is not empty; use",
                       "overwrite = TRUE to replace a folder written by",
                       "writeRich()."), path), call. = FALSE)
  }
  if (!.rich_io_is_folder(path)) {
    stop(sprintf(paste("writeRich(): folder '%s' was not written by writeRich()",
                       "(no manifest.json naming format '%s'); it is left",
                       "unchanged."), path, .rich_io_format), call. = FALSE)
  }
  TRUE
}

.rich_io_is_folder <- function(path) {
  f <- file.path(path, "manifest.json")
  if (!file.exists(f) || dir.exists(f)) return(FALSE)
  m <- tryCatch(.rich_io_read_json(f), error = function(e) NULL)
  is.list(m) && identical(m$format, .rich_io_format)
}

.rich_io_read_json <- function(file) {
  txt <- readLines(file, warn = FALSE, encoding = "UTF-8")
  jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = FALSE)
}
