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
#
# SPEC-RICHR-002 D3: `pathway` holds the annotation identifiers (T1..T5 in
# gsea_annot.csv) and the label follows in `Term`.
gsea_labels <- c(T1 = "up pathway", T2 = "down pathway", T3 = "mixed pathway",
                 T4 = "top pathway", T5 = "middle pathway")

expect_gsea_invariants <- function(res, pathways, es_sign) {
  expect_s4_class(res, "GSEAResult")
  out <- res@result
  expect_identical(
    vapply(out, function(col) class(col)[1L], character(1)),
    c(pathway = "character", pval = "numeric", padj = "numeric",
      log2err = "numeric", ES = "numeric", NES = "numeric", size = "integer",
      leadingEdge = "character", Term = "character")
  )
  expect_identical(out$Term, unname(gsea_labels[out$pathway]))
  expect_setequal(out$pathway, pathways)
  expect_identical(nrow(out), length(pathways))
  expect_identical(sign(out$ES), unname(es_sign[out$pathway]))
  ok <- !is.na(out$NES)
  expect_identical(sign(out$NES[ok]), sign(out$ES[ok]))
  expect_identical(res@pAdjustMethod, "BH")
  expect_identical(res@genenumber, 40L)
  expect_identical(res@sep, ",")
}

# SPEC-RICHR-002 D1 (REQ-015, REQ-018): pvalue = 1 keeps values equal to the
# cutoff, so T5 (pval exactly 1, ES 0) is now reported. These runs mean
# unrestricted reporting, so they pass padj = 1 with the 0.2.0 defaults.
test_that("richGSEA() keeps its result shape, pathway set and effect signs for ranked list A", {
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"), pvalue = 1, padj = 1)
  expect_gsea_invariants(
    res, c("T4", "T2", "T1", "T3", "T5"),
    c(T4 = 1, T2 = -1, T1 = 1, T3 = -1, T5 = 0)
  )
  expect_identical(res@pvalueCutoff, 1)
})

test_that("richGSEA() keeps its result shape, pathway set and effect signs for ranked list B", {
  res <- run_gsea(gsea_ranked("gsea_ranked_b.csv"), pvalue = 1, padj = 1)
  expect_gsea_invariants(
    res, c("T4", "T2", "T1", "T3", "T5"),
    c(T4 = 1, T2 = -1, T1 = 1, T3 = 1, T5 = 0)
  )
})

test_that("richGSEA() default pvalue filter keeps the three strong pathways on ranked list A", {
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"))
  expect_gsea_invariants(
    res, c("T4", "T2", "T1"),
    c(T4 = 1, T2 = -1, T1 = 1)
  )
  expect_identical(res@pvalueCutoff, 0.05)
})

test_that("richGSEA() pins NES, ES, pval, padj and leadingEdge for ranked list A", {
  .char_skip_unless_recorded_fgsea()
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"), pvalue = 1, padj = 1)
  expect_snapshot_value(pin_gsea(res), style = "deparse", tolerance = 1e-8)
})

test_that("richGSEA() pins the same pathways for the perturbed ranked list B", {
  .char_skip_unless_recorded_fgsea()
  res <- run_gsea(gsea_ranked("gsea_ranked_b.csv"), pvalue = 1, padj = 1)
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
