## data-raw/module.R -- rebuild data/module.rda (object `module`).
##
## Last run for richR 0.2.0 (SPEC-RICHR-003 M3.6). Running this script replaces
## data/module.rda with a fresh KEGG snapshot, so run it only when an
## annotation refresh is intended. The result depends on live KEGG on the day
## it runs and is not byte-reproducible.
##
## Run from the package root, with richR installed from this tree:
##   Rscript data-raw/module.R
## Needs network access and the KEGGREST package.
##
## keggList("module") returns names such as "md:M00001"; the prefix is dropped
## so the row names are the module ids ("M00001") that richR matches against.
## The object carries attr(, "richR_source") = list(retrieved, release, url);
## release is the KEGG version (see data-raw/kegg-refresh-helpers.R). The
## retrieval date, KEGG version and row count are stamped into R/data.R.

source(file.path("data-raw", "kegg-refresh-helpers.R"))
stopifnot(requireNamespace("KEGGREST", quietly = TRUE))

info <- kegg_refresh_start()
old <- kegg_bundled_object(file.path("data", "module.rda"), "module")

listing <- KEGGREST::keggList("module")
stopifnot(length(listing) > 0L, !is.null(names(listing)))

module <- data.frame(
  annotation = as.vector(listing),
  row.names = sub("^md:", "", names(listing)),
  stringsAsFactors = FALSE
)
kegg_check_not_partial(nrow(module), nrow(old), "module")
stopifnot(!anyDuplicated(rownames(module)),
          all(grepl("^M[0-9]{5}$", rownames(module))),
          !anyNA(module$annotation), all(nzchar(trimws(module$annotation))))

kegg_print_diff(setNames(old$annotation, rownames(old)),
                setNames(module$annotation, rownames(module)), "module")

attr(module, "richR_source") <- kegg_source_record(info, "https://rest.kegg.jp/list/module")
save(module, file = file.path("data", "module.rda"), compress = "xz")
kegg_stamp_documentation("module", info, nrow(module), "data-raw/module.R")

message("module: ", nrow(module), " modules retrieved ", info$retrieved,
        ", KEGG version ", info$release,
        ". Re-run roxygen2::roxygenise() to refresh man/module.Rd.")
