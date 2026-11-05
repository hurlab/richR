# GSEA plot display labels with short = TRUE (SPEC-RICHR-003 REQ-011). Fork-only:
# relies on the GSEA result layout of SPEC-RICHR-002 (pathway = identifier,
# Term = label). Two pathways whose labels differ only in a whitespace run must
# stay two categories after shortening.

gsea_ws_result <- function() {
  genes <- sprintf("g%02d", 1:40)
  annot <- data.frame(
    GeneID = c(genes[1:8], genes[33:40], genes[15:22]),
    Term = rep(c("P1", "P2", "P3"), each = 8),
    Annot = rep(c("a b c d e f g h", "a  b c d e f g h", "other"), each = 8),
    stringsAsFactors = FALSE
  )
  ranks <- stats::setNames(seq(3, -3, length.out = length(genes)), genes)
  out <- NULL
  utils::capture.output(
    out <- richGSEA(ranks, annot, minSize = 1, maxSize = 50, pvalue = 1, padj = 1,
                    seed = 7, nproc = 1)
  )
  list(result = out, ranks = ranks)
}

axis_labels <- function(p) {
  pp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  out <- character()
  for (ax in c("x", "y")) {
    sc <- pp[[ax]]
    if (!is.null(sc) && is.function(sc$get_labels)) {
      out <- c(out, tryCatch(as.character(sc$get_labels()), error = function(e) character()))
    }
  }
  out[!is.na(out)]
}

test_that("GSEA plots keep whitespace-only label variants apart with short = TRUE", {
  skip_if_not_installed("fgsea")
  # Catches: shortening re-merging two GSEA pathways into one category.
  g <- gsea_ws_result()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- list(
    richLollipop = richLollipop(g$result, pvalue = 1, padj = 1, short = TRUE),
    richNES = richNES(g$result, pvalue = 1, padj = 1, short = TRUE),
    richGeneDot = richGeneDot(g$result, fc = g$ranks, pvalue = 1, padj = 1, short = TRUE),
    richECDF = richECDF(g$result, pvalue = 1, padj = 1, short = TRUE)
  )
  for (what in names(plots)) {
    cats <- axis_labels(plots[[what]])
    expect_equal(length(unique(cats[startsWith(cats, "a")])), 2, info = what)
  }
})
