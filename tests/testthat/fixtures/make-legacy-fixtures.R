# Regenerates the legacy fixtures used by test-tested-legacy.R and
# test-rich-io.R: result objects serialized by a richR build that predates the
# `tested` slot (richR 0.2.0 before M4, commit fca5276), plus the outputs that
# build computed from them. No network.
#
# Run from tests/testthat/fixtures/ against an installed pre-M4 build:
#   R_LIBS=<lib with the fca5276 install> Rscript --vanilla make-legacy-fixtures.R
# The committed .rds files are the source of truth for the tests; this script
# documents how they were made so a reviewer can reproduce them. It refuses to
# run against a build whose richResult class already has a `tested` slot,
# because such a build cannot produce a legacy object.

suppressPackageStartupMessages(library(richR))

if ("tested" %in% methods::slotNames("richResult") ||
    "tested" %in% methods::slotNames("GSEAResult")) {
  stop("make-legacy-fixtures.R must run against a richR build without the ",
       "`tested` slot (fca5276); the installed build has it.", call. = FALSE)
}

## ---- F-ORA: 20-gene universe, five gene sets in four roles ------------------
genes <- sprintf("G%02d", 1:20)
mk <- function(term, label, idx) {
  data.frame(GeneID = genes[idx], Term = term, Annot = label,
             stringsAsFactors = FALSE)
}
ora_annot <- rbind(
  mk("T1", "reported set", 1:5),
  mk("T2", "tested set failing the cutoff", c(6:10, 12:14)),
  mk("T3", "eligible set without overlap", 18:20),
  mk("T4", "set below minGSSize", c(1, 15)),
  mk("T5", "tested set below minSize", c(11, 15:17))
)
ora_query <- genes[c(1:7, 11)]
legacy_ora <- suppressMessages(enrich(
  ora_query, ora_annot, pvalue = 1, padj = NULL, minSize = 2, minGSSize = 3,
  keepRich = FALSE
))

## ---- F-GSEA: the characterization GSEA fixtures, fixed seed -----------------
read_fixture <- function(name, character = FALSE) {
  utils::read.csv(name, stringsAsFactors = FALSE,
                  colClasses = if (character) "character" else NA)
}
ranked <- read_fixture("gsea_ranked_a.csv")
ranked <- stats::setNames(ranked$score, ranked$GeneID)
gsea_annot <- read_fixture("gsea_annot.csv", character = TRUE)
legacy_gsea <- NULL
invisible(utils::capture.output(
  legacy_gsea <- richGSEA(ranked, gsea_annot, minSize = 3, maxSize = 20,
                          pvalue = 1, padj = NULL, seed = 42, nproc = 1)
))

stopifnot(nrow(legacy_ora@result) > 0L, nrow(legacy_gsea@result) > 0L,
          !methods::.hasSlot(legacy_ora, "tested"),
          !methods::.hasSlot(legacy_gsea, "tested"),
          identical(enrichmentInfo(legacy_ora)$schema_version, 1L),
          identical(enrichmentInfo(legacy_gsea)$schema_version, 1L))

## ---- outputs of the pre-M4 build on those objects ---------------------------
reference <- list(
  richR_version = as.character(utils::packageVersion("richR")),
  richResult = list(
    result = result(legacy_ora),
    compareResult = compareResult(list(a = legacy_ora)),
    filter = dplyr::filter(legacy_ora, Pvalue <= 0.5)@result,
    richBar = richBar(legacy_ora)$data
  ),
  GSEAResult = list(
    result = result(legacy_gsea),
    filter = dplyr::filter(legacy_gsea, pval <= 0.5)@result
  )
)

saveRDS(legacy_ora, "legacy_richResult_fca5276.rds")
saveRDS(legacy_gsea, "legacy_GSEAResult_fca5276.rds")
saveRDS(reference, "legacy_reference_fca5276.rds")
