## data-raw/kegg.R -- rebuild data/kegg.rdata (object `kegg.db`).
##
## Last run for richR 0.2.0 (SPEC-RICHR-003 M3.6). Running this script replaces
## data/kegg.rdata with a fresh KEGG snapshot, so run it only when an
## annotation refresh is intended. The result depends on live KEGG on the day
## it runs and is not byte-reproducible.
##
## Run from the package root, with richR installed from this tree:
##   Rscript data-raw/kegg.R
## Needs network access and the KEGGREST package.
##
## The transformation is the one richR:::.get_kg_dat(builtin = FALSE) applies
## to the live KEGG listing: names such as "path:map01100" lose everything up
## to "map", leaving the five-digit pathway id as the row name. The object
## carries attr(, "richR_source") = list(retrieved, release, url); release is
## the KEGG version (see data-raw/kegg-refresh-helpers.R). The retrieval date,
## KEGG version and row count are stamped into R/data.R.

source(file.path("data-raw", "kegg-refresh-helpers.R"))
stopifnot(requireNamespace("KEGGREST", quietly = TRUE))

info <- kegg_refresh_start()
old <- kegg_bundled_object(file.path("data", "kegg.rdata"), "kegg.db")

listing <- KEGGREST::keggList("pathway")
stopifnot(length(listing) > 0L, !is.null(names(listing)))

kegg.db <- data.frame(
  annotation = as.vector(listing),
  row.names = sub(".*map", "", names(listing)),
  stringsAsFactors = FALSE
)
kegg_check_not_partial(nrow(kegg.db), nrow(old), "kegg.db")
stopifnot(!anyDuplicated(rownames(kegg.db)),
          all(grepl("^[0-9]{5}$", rownames(kegg.db))),
          !anyNA(kegg.db$annotation), all(nzchar(trimws(kegg.db$annotation))),
          "04148" %in% rownames(kegg.db))

kegg_print_diff(setNames(old$annotation, rownames(old)),
                setNames(kegg.db$annotation, rownames(kegg.db)), "kegg.db")

attr(kegg.db, "richR_source") <- kegg_source_record(info, "https://rest.kegg.jp/list/pathway")
save(kegg.db, file = file.path("data", "kegg.rdata"), compress = "xz")
kegg_stamp_documentation("kegg.db", info, nrow(kegg.db), "data-raw/kegg.R")

message("kegg.db: ", nrow(kegg.db), " pathways retrieved ", info$retrieved,
        ", KEGG version ", info$release,
        ". Re-run roxygen2::roxygenise() to refresh man/kegg.db.Rd.")
