# group_by() on richResult, GSEAResult and Annot objects.
#
# The methods passed `add = add` to dplyr::group_by(), an argument dplyr
# (>= 1.0.0) made defunct, so every call failed with "The `add` argument of
# `group_by()` was deprecated in dplyr 1.0.0 and is now defunct." These tests
# catch that failure and an `add = TRUE` that replaces the existing groups
# instead of adding to them.
#
# Self-contained on purpose: inline fixtures, public API only, valid under
# testthat editions 2 and 3.

gb_plain <- function(df) {
  df <- as.data.frame(dplyr::ungroup(df), stringsAsFactors = FALSE)
  # group_by() drops row names, so compare the tables without them.
  rownames(df) <- NULL
  df
}

gb_ora <- function() {
  genes <- sprintf("G%02d", 1:20)
  annot <- data.frame(
    GeneID = genes[c(1:6, 4:11, 12:20)],
    Term = rep(c("S1", "S2", "S3"), c(6, 8, 9)),
    stringsAsFactors = FALSE
  )
  res <- suppressMessages(richR::enrich(
    genes[c(1:5, 8:10, 13)], annot, pvalue = 1, padj = NULL,
    minSize = 1, minGSSize = 1
  ))
  # A grouped view of an empty table would assert nothing.
  expect_true(nrow(res@result) >= 2)
  res
}

gb_gsea <- function() {
  genes <- sprintf("G%02d", 1:30)
  scores <- stats::setNames(seq(3, -2.8, by = -0.2), genes)
  annot <- data.frame(
    GeneID = genes[c(1:6, 12:18, 24:30)],
    Term = rep(c("P1", "P2", "P3"), c(6, 7, 7)),
    stringsAsFactors = FALSE
  )
  res <- NULL
  utils::capture.output(
    res <- richR::richGSEA(scores, annot, pvalue = 1, padj = NULL,
                           minSize = 3, maxSize = 50)
  )
  expect_true(nrow(res@result) >= 2)
  res
}

test_that("group_by() groups the result table of a richResult", {
  r <- gb_ora()
  g <- dplyr::group_by(r, Term)
  expect_true(methods::is(g, "richResult"))
  expect_true(inherits(g@result, "grouped_df"))
  expect_identical(dplyr::group_vars(g@result), "Term")
  expect_identical(gb_plain(g@result), gb_plain(r@result))

  g0 <- dplyr::group_by(r, Term, add = FALSE)
  expect_identical(dplyr::group_vars(g0@result), "Term")
  expect_identical(gb_plain(g0@result), gb_plain(r@result))
})

test_that("group_by(add = TRUE) adds to the groups of a richResult", {
  r <- gb_ora()
  g <- dplyr::group_by(dplyr::group_by(r, Term), Annot, add = TRUE)
  expect_true(methods::is(g, "richResult"))
  expect_identical(dplyr::group_vars(g@result), c("Term", "Annot"))
  expect_identical(gb_plain(g@result), gb_plain(r@result))
})

test_that("group_by() groups the result table of a GSEAResult", {
  skip_if_not_installed("fgsea")
  r <- gb_gsea()
  g <- dplyr::group_by(r, Term)
  expect_true(methods::is(g, "GSEAResult"))
  expect_true(inherits(g@result, "grouped_df"))
  expect_identical(dplyr::group_vars(g@result), "Term")
  expect_identical(gb_plain(g@result), gb_plain(r@result))
})

test_that("group_by() groups the annotation table of an Annot object", {
  a <- methods::new(
    "Annot", species = "custom", anntype = "custom", keytype = "SYMBOL",
    annot = data.frame(GeneID = c("A", "B", "C", "A"),
                       Term = c("T1", "T1", "T2", "T2"),
                       stringsAsFactors = FALSE)
  )
  g <- dplyr::group_by(a, Term)
  expect_true(methods::is(g, "Annot"))
  expect_true(inherits(g@annot, "grouped_df"))
  expect_identical(dplyr::group_vars(g@annot), "Term")
  expect_identical(gb_plain(g@annot), gb_plain(a@annot))

  g2 <- dplyr::group_by(g, GeneID, add = TRUE)
  expect_identical(dplyr::group_vars(g2@annot), c("Term", "GeneID"))
  expect_identical(gb_plain(g2@annot), gb_plain(a@annot))
})
