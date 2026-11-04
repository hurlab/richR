# Characterization baseline (REQ-004): the S4 accessor, detail(), rbind() and
# coercion surface of richResult, GSEAResult and Annot as produced by richR
# 0.1.5. These are exactly the registrations that R/misc.R and R/accessors.R
# define twice; the tests exercise whichever definition is live, so that
# removing a duplicate (M1.3) is proven equivalent or visibly different.

classes <- function() {
  list(
    richResult = .char_richResult(),
    GSEAResult = .char_gsea_result(),
    Annot = .char_annot()
  )
}

test_that("names/colnames/rownames/row.names/dim are pinned for all three classes", {
  pins <- lapply(classes(), function(x) {
    list(
      names = names(x), colnames = colnames(x), rownames = rownames(x),
      row.names = row.names(x), dim = dim(x)
    )
  })
  expect_snapshot_value(pins, style = "deparse")
})

test_that("head() and tail() print their banner and return the leading/trailing rows", {
  for (cls in names(classes())) {
    x <- classes()[[cls]]
    expect_snapshot({
      cat("class:", cls, "\n")
      print(head(x, 2))
      print(tail(x, 2))
    })
  }
})

test_that("[ and $ subset the underlying data.frame for all three classes", {
  pins <- lapply(classes(), function(x) {
    list(
      rows_cols = x[1:2, 1:2],
      first_col = x[, 1],
      single_cell = x[1, 2],
      one_row_drop_false = x[1, 1:2, drop = FALSE],
      dollar = x$Annot
    )
  })
  expect_snapshot_value(pins, style = "deparse")
  rich <- .char_richResult()
  expect_identical(rich$Term, c("apoptosis", "cell cycle", "signaling"))
  expect_identical(.char_gsea_result()$pathway, c("p1", "p2", "p3"))
  expect_identical(.char_annot()$Annot[1], "apoptosis")
})

test_that("result() and summary() are pinned for richResult and GSEAResult", {
  expect_snapshot_value(
    list(rich = result(.char_richResult()), gsea = result(.char_gsea_result())),
    style = "deparse"
  )
  expect_snapshot({
    summary(.char_richResult())
    summary(.char_gsea_result())
  })
})

test_that("detail() is pinned for richResult; for GSEAResult it messages and returns result", {
  expect_snapshot_value(detail(.char_richResult()), style = "deparse")
  # REQ-010, the one deliberate behaviour change of M1: in 0.1.5 the live
  # detail,GSEAResult method read a slot the class does not have and errored.
  # It now emits a message and returns as.data.frame(x@result).
  gsea <- .char_gsea_result()
  expect_snapshot(returned <- detail(gsea))
  expect_message(detail(gsea), "GSEAResult does not have a detail slot")
  expect_identical(returned, as.data.frame(gsea@result))
  expect_s3_class(returned, "data.frame")
})

test_that("as.data.frame() and as(<class>, 'data.frame') pin the three coercions", {
  pins <- lapply(classes(), function(x) {
    list(as.data.frame = as.data.frame(x), as = methods::as(x, "data.frame"))
  })
  expect_snapshot_value(pins, style = "deparse")
})

test_that("as(<data.frame>, 'Annot') and as(<data.frame>, 'richResult') are pinned", {
  annot_df <- data.frame(
    GeneID = c("A", "B", "C"), Term = c("T1", "T1", "T2"),
    Annot = c("first", "first", "second"), stringsAsFactors = FALSE
  )
  as_annot <- methods::as(annot_df, "Annot")
  expect_s4_class(as_annot, "Annot")

  rich_df <- .char_richResult()@result
  as_rich <- methods::as(rich_df, "richResult")
  expect_s4_class(as_rich, "richResult")
  expect_snapshot_value(list(annot = as_annot, richResult = as_rich), style = "deparse")
})

test_that("rbind() on richResult and GSEAResult objects is pinned", {
  rich_a <- .char_richResult()
  rich_b <- .char_richResult_second()
  gsea_a <- .char_gsea_result()
  gsea_b <- .char_gsea_result_second()
  expect_snapshot_value(
    list(
      rich_distinct = rbind(rich_a, rich_b),
      rich_duplicate_rownames = rbind(rich_a, rich_a),
      gsea_distinct = rbind(gsea_a, gsea_b),
      gsea_duplicate_rownames = rbind(gsea_a, gsea_a)
    ),
    style = "deparse", tolerance = 1e-8
  )
})

test_that("rbind() keeps argument names as row-name prefixes, as in 0.1.5", {
  rich_a <- .char_richResult()
  rich_b <- .char_richResult_second()
  gsea_a <- .char_gsea_result()
  gsea_b <- .char_gsea_result_second()

  rich_named <- rbind(first = rich_a, second = rich_b)
  expect_identical(
    rownames(rich_named),
    c("first.GO:0001", "first.GO:0002", "first.GO:0003",
      "second.GO:0010", "second.GO:0011")
  )
  expect_identical(
    rich_named,
    rbind(first = as.data.frame(rich_a), second = as.data.frame(rich_b))
  )

  gsea_named <- rbind(first = gsea_a, second = gsea_b)
  expect_identical(
    rownames(gsea_named),
    c("first.1", "first.2", "first.3", "second.1", "second.2")
  )
  expect_identical(
    gsea_named,
    rbind(first = as.data.frame(gsea_a), second = as.data.frame(gsea_b))
  )
})

test_that("rbind methods called with no arguments return NULL, as in 0.1.5", {
  # base::rbind() never dispatches with zero arguments, so the methods are
  # called directly. In 0.1.5 do.call(rbind, list()) returned NULL; a
  # bindROWS(objects[[1L]], ...) body errors with "subscript out of bounds".
  expect_null(richR:::rbind.richResult())
  expect_null(richR:::rbind.GSEAResult())
})
