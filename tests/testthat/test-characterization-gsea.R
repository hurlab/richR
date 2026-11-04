# Characterization baseline (REQ-003): GSEA NES / ES / pval / padj / leadingEdge
# for a fixture ranked list under a fixed seed, as produced by richR 0.1.5.
# The fgsea version is part of every snapshot (key `fgsea_version`) so a diff
# caused by an fgsea upgrade can be told apart from one caused by richR.

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

test_that("richGSEA() pins NES, ES, pval, padj and leadingEdge for ranked list A", {
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"), pvalue = 1)
  expect_true(all(c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge") %in%
                    names(res@result)))
  expect_snapshot_value(pin_gsea(res), style = "deparse", tolerance = 1e-8)
})

test_that("richGSEA() pins the same pathways for the perturbed ranked list B", {
  res <- run_gsea(gsea_ranked("gsea_ranked_b.csv"), pvalue = 1)
  expect_snapshot_value(pin_gsea(res), style = "deparse", tolerance = 1e-8)
})

test_that("richGSEA() pins the default pvalue filter on ranked list A", {
  res <- run_gsea(gsea_ranked("gsea_ranked_a.csv"))
  expect_snapshot_value(pin_gsea(res), style = "deparse", tolerance = 1e-8)
})
