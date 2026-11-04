# Characterization baseline (REQ-003): GSEA NES / ES / pval / padj / leadingEdge
# for a fixture ranked list under a fixed seed, as produced by richR 0.1.5.
# The fgsea version is part of every exact-numeric snapshot (key
# `fgsea_version`) so a diff caused by an fgsea upgrade can be told apart from
# one caused by richR. Those snapshots run only on the recorded fgsea version
# (.char_fgsea_recorded in helper-characterization.R) and are skipped with a
# message otherwise; the version-independent invariants below always run.

skip_if_not_installed("fgsea")

gsea_ranked <- function(name) {
  df <- .char_read(name)
  stats::setNames(df$score, df$GeneID)
}

gsea_annot <- function() .char_read("gsea_annot.csv", character = TRUE)

run_gsea <- function(ranked, ...) {
  # fgsea draws a text progress bar on stdout; capture it so it does not
  # interleave with the test reporter. Warnings are deliberately not muffled.
  out <- NULL
  utils::capture.output(
    out <- richGSEA(ranked, gsea_annot(), minSize = 3, maxSize = 20,
                    seed = 42, nproc = 1, ...)
  )
  out
}

pin_gsea <- function(res) {
  expect_s4_class(res, "GSEAResult")
  list(
    fgsea_version = as.character(utils::packageVersion("fgsea")),
    result = .char_pin_df(res@result),
    slots = .char_slots(res)
  )
}

# Version-independent invariants: result columns and classes, the set of
# pathways, the row count, the sign of ES and NES, and the non-fgsea slots.
# None of these depends on fgsea's permutation draws.
expect_gsea_invariants <- function(res, pathways, es_sign) {
  expect_s4_class(res, "GSEAResult")
  out <- res@result
  expect_identical(
    vapply(out, function(col) class(col)[1L], character(1)),
    c(pathway = "character", pval = "numeric", padj = "numeric",
      log2err = "numeric", ES = "numeric", NES = "numeric", size = "integer",
      leadingEdge = "character")
  )
  expect_setequal(out$pathway, pathways)
  expect_identical(nrow(out), length(pathways))
  expect_identical(sign(out$ES), unname(es_sign[out$pathway]))
  ok <- !is.na(out$NES)
  expect_identical(sign(out$NES[ok]), sign(out$ES[ok]))
  expect_identical(res@pAdjustMethod, "BH")
  expect_identical(res@genenumber, 40L)
  expect_identical(res@sep, ",")
}

test_that("richGSEA() keeps its result shape, pathway set and effect signs for ranked list A", {
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"), pvalue = 1)
  expect_gsea_invariants(
    res, c("top pathway", "down pathway", "up pathway", "mixed pathway"),
    c("top pathway" = 1, "down pathway" = -1, "up pathway" = 1, "mixed pathway" = -1)
  )
  expect_identical(res@pvalueCutoff, 1)
})

test_that("richGSEA() keeps its result shape, pathway set and effect signs for ranked list B", {
  res <- run_gsea(gsea_ranked("gsea_ranked_b.csv"), pvalue = 1)
  expect_gsea_invariants(
    res, c("top pathway", "down pathway", "up pathway", "mixed pathway"),
    c("top pathway" = 1, "down pathway" = -1, "up pathway" = 1, "mixed pathway" = 1)
  )
})

test_that("richGSEA() default pvalue filter keeps the three strong pathways on ranked list A", {
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"))
  expect_gsea_invariants(
    res, c("top pathway", "down pathway", "up pathway"),
    c("top pathway" = 1, "down pathway" = -1, "up pathway" = 1)
  )
  expect_identical(res@pvalueCutoff, 0.05)
})

test_that("richGSEA() pins NES, ES, pval, padj and leadingEdge for ranked list A", {
  .char_skip_unless_recorded_fgsea()
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"), pvalue = 1)
  expect_snapshot_value(pin_gsea(res), style = "deparse", tolerance = 1e-8)
})

test_that("richGSEA() pins the same pathways for the perturbed ranked list B", {
  .char_skip_unless_recorded_fgsea()
  res <- run_gsea(gsea_ranked("gsea_ranked_b.csv"), pvalue = 1)
  expect_snapshot_value(pin_gsea(res), style = "deparse", tolerance = 1e-8)
})

test_that("richGSEA() pins the default pvalue filter on ranked list A", {
  .char_skip_unless_recorded_fgsea()
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"))
  expect_snapshot_value(pin_gsea(res), style = "deparse", tolerance = 1e-8)
})

test_that("detail() on a fixture-built GSEAResult messages and returns the result table", {
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"), pvalue = 1)
  expect_message(returned <- detail(res), "does not have a detail slot")
  expect_identical(returned, as.data.frame(res@result))
  expect_gt(nrow(returned), 0L)
})
