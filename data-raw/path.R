## data-raw/path.R -- rebuild data/path.rda (object `path`).
##
## Last run for richR 0.2.0 (SPEC-RICHR-003 M3.6). Running this script replaces
## data/path.rda with a fresh KEGG snapshot, so run it only when an annotation
## refresh is intended. The result depends on live KEGG on the day it runs and
## is not byte-reproducible. How the snapshot bundled before 0.2.0 was built is
## unknown; this script reads the KEGG BRITE pathway hierarchy br08901.
##
## Run from the package root, with richR installed from this tree:
##   Rscript data-raw/path.R
## Needs network access and the jsonlite package.
##
## Target structure (what richR code reads): a data.frame with the character
## columns ko (five-digit pathway id), Level3 (pathway name), Level2
## (sub-category) and Level1 (category), with ko as the row names.
##
## Consistency (SPEC-RICHR-003 REQ-020): one row per pathway. KEGG BRITE lists
## a few pathways under several hierarchy rows; the first row in BRITE order is
## kept (the row richKEGG and richGSEA attach), so richLevel, richKEGG and
## richGSEA agree. The collapsed pathways are printed. Identifiers must be
## unique and five-digit, and Level1-Level3 non-blank, before anything is
## saved. The diff against the bundled object is printed.
##
## The object carries attr(, "richR_source") = list(retrieved, release, url);
## release is the KEGG version (see data-raw/kegg-refresh-helpers.R). The
## retrieval date, KEGG version and row count are stamped into R/data.R.

source(file.path("data-raw", "kegg-refresh-helpers.R"))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))

info <- kegg_refresh_start()
old <- kegg_bundled_object(file.path("data", "path.rda"), "path")

url <- "https://rest.kegg.jp/get/br:br08901/json"
json_file <- tempfile(fileext = ".json")
utils::download.file(url, json_file, quiet = TRUE, mode = "wb")
tree <- jsonlite::fromJSON(json_file, simplifyVector = FALSE)
stopifnot(length(tree$children) > 0L)

## Walk the Level1 -> Level2 -> Level3 hierarchy. Leaves look like
## "01100  Metabolic pathways" (id, two spaces, name).
strip_code <- function(x) sub("^[0-9]+\\s+", "", x)
rows <- list()
unparsed <- character(0)
for (l1 in tree$children) {
  for (l2 in l1$children) {
    for (l3 in l2$children) {
      if (length(l3$children) > 0L) {
        stop("unexpected BRITE depth below ", l3$name, call. = FALSE)
      }
      parts <- regmatches(l3$name, regexec("^([0-9]{5})\\s+(.*)$", l3$name))[[1L]]
      if (length(parts) != 3L) {
        unparsed <- c(unparsed, l3$name)
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        ko = parts[2L], Level3 = trimws(parts[3L]),
        Level2 = trimws(strip_code(l2$name)),
        Level1 = trimws(strip_code(l1$name)),
        stringsAsFactors = FALSE
      )
    }
  }
}
if (length(unparsed) > 0L) {
  stop("BRITE leaves without a five-digit pathway id: ",
       paste(unparsed, collapse = " | "), call. = FALSE)
}
stopifnot(length(rows) > 0L)
path <- do.call(rbind, rows)

## One row per pathway: the first in BRITE order.
repeated <- unique(path$ko[duplicated(path$ko)])
cat(sprintf("== BRITE rows collapsed to the first: %d pathway(s)\n", length(repeated)))
for (ko in repeated) {
  hits <- path[path$ko == ko, , drop = FALSE]
  cat(sprintf("  %s kept:    %s | %s | %s\n", ko, hits$Level3[1], hits$Level2[1], hits$Level1[1]))
  for (i in seq_len(nrow(hits))[-1L]) {
    cat(sprintf("  %s dropped: %s | %s | %s\n", ko, hits$Level3[i], hits$Level2[i], hits$Level1[i]))
  }
}
path <- path[!duplicated(path$ko), , drop = FALSE]
rownames(path) <- path$ko

kegg_check_not_partial(nrow(path), length(unique(old$ko)), "path")
blank <- function(x) is.na(x) | !nzchar(trimws(x))
stopifnot(identical(names(path), c("ko", "Level3", "Level2", "Level1")),
          !anyDuplicated(path$ko),
          all(grepl("^[0-9]{5}$", path$ko)),
          !any(blank(path$Level3)), !any(blank(path$Level2)), !any(blank(path$Level1)))

## Diff against the bundled object; its duplicated pathways are reported and
## compared by their first row.
old_repeated <- unique(old$ko[duplicated(old$ko)])
cat(sprintf("== bundled path: %d rows, %d pathway(s) with several rows%s\n", nrow(old),
            length(old_repeated),
            if (length(old_repeated)) paste0(": ", paste(old_repeated, collapse = ", ")) else ""))
old_first <- old[!duplicated(old$ko), , drop = FALSE]
for (level in c("Level3", "Level2", "Level1")) {
  kegg_print_diff(setNames(old_first[[level]], old_first$ko),
                  setNames(path[[level]], path$ko), paste0("path ", level))
}

attr(path, "richR_source") <- kegg_source_record(info, url)
save(path, file = file.path("data", "path.rda"), compress = "gzip")
kegg_stamp_documentation("path", info, nrow(path), "data-raw/path.R")

message("path: ", nrow(path), " pathways retrieved ", info$retrieved,
        ", KEGG version ", info$release,
        ". Re-run roxygen2::roxygenise() to refresh man/path.Rd.")
