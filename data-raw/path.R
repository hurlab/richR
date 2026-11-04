## data-raw/path.R -- rebuild data/path.rda (object `path`).
##
## NOT run as part of richR 0.2.0 M1: the bundled snapshot in data/ is left
## byte-identical. Running this script replaces data/path.rda with a fresh
## KEGG snapshot, so run it only when the annotation refresh is intended.
## The date the currently bundled snapshot was retrieved is unknown, and how it
## was built is unknown too: this script is a documented reconstruction, not
## the original. The JSON layout of the BRITE response assumed below follows the
## KEGG REST API manual (get/br:br08901/json: nested "name"/"children") and was
## NOT checked against a live response in M1 -- verify the stopifnot() guards
## pass, and compare against the bundled object, on the first real run.
##
## Run from the package root:  Rscript data-raw/path.R
## Needs network access and the jsonlite package.
##
## Target structure (what richR code reads): a data.frame with the character
## columns ko (five-digit pathway id), Level3 (pathway name), Level2
## (sub-category) and Level1 (category), with ko as the row names.

stopifnot(file.exists("DESCRIPTION"), requireNamespace("jsonlite", quietly = TRUE))

retrieved <- Sys.Date()
url <- "https://rest.kegg.jp/get/br:br08901/json"
tree <- jsonlite::fromJSON(url, simplifyVector = FALSE)

## Walk the Level1 -> Level2 -> Level3 hierarchy. Leaves look like
## "01100  Metabolic pathways" (id, two spaces, name).
rows <- list()
for (l1 in tree$children) {
  for (l2 in l1$children) {
    for (l3 in l2$children) {
      parts <- regmatches(l3$name, regexec("^([0-9]{5})\\s+(.*)$", l3$name))[[1L]]
      if (length(parts) != 3L) next
      rows[[length(rows) + 1L]] <- data.frame(
        ko = parts[2L], Level3 = parts[3L],
        Level2 = sub("^[0-9]+\\s+", "", l2$name),
        Level1 = sub("^[0-9]+\\s+", "", l1$name),
        stringsAsFactors = FALSE
      )
    }
  }
}
stopifnot(length(rows) > 0L)

path <- do.call(rbind, rows)
path <- path[!duplicated(path$ko), , drop = FALSE]
rownames(path) <- path$ko
stopifnot(all(c("ko", "Level3", "Level2", "Level1") %in% names(path)),
          all(grepl("^[0-9]{5}$", path$ko)))

save(path, file = file.path("data", "path.rda"), compress = "gzip")

## Stamp the retrieval date into the documentation of the rebuilt object.
doc <- readLines(file.path("R", "data.R"))
line <- grep("^#'\\s+Retrieval date \\(path\\):", doc)
stopifnot(length(line) == 1L)
doc[line] <- sprintf("#'   Retrieval date (path): %s (rebuilt by data-raw/path.R).", retrieved)
writeLines(doc, file.path("R", "data.R"))
message("path: ", nrow(path), " pathways retrieved ", retrieved,
        ". Re-run devtools::document() to refresh man/path.Rd.")
