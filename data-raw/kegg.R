## data-raw/kegg.R -- rebuild data/kegg.rdata (object `kegg.db`).
##
## NOT run as part of richR 0.2.0 M1: the bundled snapshot in data/ is left
## byte-identical. Running this script replaces data/kegg.rdata with a fresh
## KEGG snapshot, so run it only when the annotation refresh is intended.
## The date the currently bundled snapshot was retrieved is unknown.
##
## Run from the package root:  Rscript data-raw/kegg.R
## Needs network access and the KEGGREST package.
##
## The transformation is the one richR:::.get_kg_dat(builtin = FALSE) applies
## to the live KEGG listing: names such as "path:map01100" lose everything up
## to "map", leaving the five-digit pathway id as the row name.

stopifnot(file.exists("DESCRIPTION"), requireNamespace("KEGGREST", quietly = TRUE))

retrieved <- Sys.Date()
listing <- KEGGREST::keggList("pathway")
stopifnot(length(listing) > 0L, !is.null(names(listing)))

kegg.db <- data.frame(
  annotation = as.vector(listing),
  row.names = sub(".*map", "", names(listing)),
  stringsAsFactors = FALSE
)
stopifnot(nrow(kegg.db) > 0L, !anyDuplicated(rownames(kegg.db)))

save(kegg.db, file = file.path("data", "kegg.rdata"), compress = "xz")

## Stamp the retrieval date into the documentation of the rebuilt object.
doc <- readLines(file.path("R", "data.R"))
line <- grep("^#'\\s+Retrieval date \\(kegg\\.db\\):", doc)
stopifnot(length(line) == 1L)
doc[line] <- sprintf("#'   Retrieval date (kegg.db): %s (rebuilt by data-raw/kegg.R).", retrieved)
writeLines(doc, file.path("R", "data.R"))
message("kegg.db: ", nrow(kegg.db), " pathways retrieved ", retrieved,
        ". Re-run devtools::document() to refresh man/kegg.db.Rd.")
