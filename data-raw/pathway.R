## data-raw/pathway.R -- rebuild data/pathway.rda (object `pathway`).
##
## NOT run as part of richR 0.2.0 M1: the bundled snapshot in data/ is left
## byte-identical. Running this script replaces data/pathway.rda with a fresh
## KEGG snapshot, so run it only when the annotation refresh is intended.
## The date the currently bundled snapshot was retrieved is unknown, and so is
## its origin: it has 562 rows (kegg.db has 567) and richR code never loads it.
## This script assumes it is the same listing as kegg.db (keggList("pathway")),
## which is what its structure suggests, not something that was verified.
##
## Run from the package root:  Rscript data-raw/pathway.R
## Needs network access and the KEGGREST package.

stopifnot(file.exists("DESCRIPTION"), requireNamespace("KEGGREST", quietly = TRUE))

retrieved <- Sys.Date()
listing <- KEGGREST::keggList("pathway")
stopifnot(length(listing) > 0L, !is.null(names(listing)))

pathway <- data.frame(
  annotation = as.vector(listing),
  row.names = sub(".*map", "", names(listing)),
  stringsAsFactors = FALSE
)
stopifnot(nrow(pathway) > 0L, !anyDuplicated(rownames(pathway)))

save(pathway, file = file.path("data", "pathway.rda"), compress = "xz")

## Stamp the retrieval date into the documentation of the rebuilt object.
doc <- readLines(file.path("R", "data.R"))
line <- grep("^#'\\s+Retrieval date \\(pathway\\):", doc)
stopifnot(length(line) == 1L)
doc[line] <- sprintf("#'   Retrieval date (pathway): %s (rebuilt by data-raw/pathway.R).", retrieved)
writeLines(doc, file.path("R", "data.R"))
message("pathway: ", nrow(pathway), " pathways retrieved ", retrieved,
        ". Re-run devtools::document() to refresh man/pathway.Rd.")
