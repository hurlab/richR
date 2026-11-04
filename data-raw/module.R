## data-raw/module.R -- rebuild data/module.rda (object `module`).
##
## NOT run as part of richR 0.2.0 M1: the bundled snapshot in data/ is left
## byte-identical. Running this script replaces data/module.rda with a fresh
## KEGG snapshot, so run it only when the annotation refresh is intended.
## The date the currently bundled snapshot was retrieved is unknown.
##
## Run from the package root:  Rscript data-raw/module.R
## Needs network access and the KEGGREST package.
##
## keggList("module") returns names such as "md:M00001"; the prefix is dropped
## so the row names are the module ids ("M00001") that richR matches against.

stopifnot(file.exists("DESCRIPTION"), requireNamespace("KEGGREST", quietly = TRUE))

retrieved <- Sys.Date()
listing <- KEGGREST::keggList("module")
stopifnot(length(listing) > 0L, !is.null(names(listing)))

module <- data.frame(
  annotation = as.vector(listing),
  row.names = sub("^md:", "", names(listing)),
  stringsAsFactors = FALSE
)
stopifnot(nrow(module) > 0L, !anyDuplicated(rownames(module)),
          all(grepl("^M[0-9]{5}$", rownames(module))))

save(module, file = file.path("data", "module.rda"), compress = "xz")

## Stamp the retrieval date into the documentation of the rebuilt object.
doc <- readLines(file.path("R", "data.R"))
line <- grep("^#'\\s+Retrieval date \\(module\\):", doc)
stopifnot(length(line) == 1L)
doc[line] <- sprintf("#'   Retrieval date (module): %s (rebuilt by data-raw/module.R).", retrieved)
writeLines(doc, file.path("R", "data.R"))
message("module: ", nrow(module), " modules retrieved ", retrieved,
        ". Re-run devtools::document() to refresh man/module.Rd.")
