# Characterization baseline (REQ-001, REQ-002): the full ORA result surface of
# enrich(), richGO() and richKEGG() on small committed fixtures, as produced by
# richR 0.1.5. Every test asserts a non-empty result before snapshotting so an
# empty sweep can never pass as a pin.

ora_query <- function() .char_read("ora_query.csv")$GeneID
ora_annot <- function() .char_read("ora_annot.csv", character = TRUE)

pin_ora <- function(res) {
  expect_s4_class(res, "richResult")
  .char_nonempty(res@result)
  .char_nonempty(res@detail)
  list(
    result = .char_pin_df(res@result),
    detail = .char_pin_df(res@detail),
    slots = .char_slots(res)
  )
}

test_that("enrich() pins the full ORA result surface on the fixture annotation", {
  res <- suppressMessages(enrich(
    ora_query(), ora_annot(), pvalue = 0.5, minSize = 2, minGSSize = 3,
    keepRich = TRUE
  ))
  expect_snapshot_value(pin_ora(res), style = "deparse", tolerance = 1e-8)
})

test_that("enrich() pins the detail table when x is a DEG data.frame", {
  degs <- .char_read("ora_degs.csv")
  x <- data.frame(lfc = degs$lfc, row.names = degs$GeneID)
  res <- suppressMessages(enrich(
    x, ora_annot(), pvalue = 0.5, minSize = 2, minGSSize = 3, keepRich = TRUE
  ))
  expect_s4_class(res, "richResult")
  .char_nonempty(res@detail)
  expect_true("lfc" %in% names(res@detail))
  expect_snapshot_value(
    list(result = .char_pin_df(res@result), detail = .char_pin_df(res@detail)),
    style = "deparse", tolerance = 1e-8
  )
})

test_that("richGO() pins the ORA surface on offline GO-style annotation", {
  annot <- ora_annot()
  go <- data.frame(
    GeneID = annot$GeneID, GOALL = annot$Term, Annot = annot$Annot,
    ONTOLOGYALL = "BP", stringsAsFactors = FALSE
  )
  res <- suppressMessages(richGO(
    ora_query(), go, ontology = "BP", pvalue = 0.5, minSize = 2, minGSSize = 3
  ))
  expect_snapshot_value(pin_ora(res), style = "deparse", tolerance = 1e-8)
})

test_that("richKEGG() pins the ORA surface with bundled KEGG names and levels", {
  res <- suppressMessages(richKEGG(
    ora_query(), .char_read("kegg_annot.csv", character = TRUE),
    pvalue = 0.5, minSize = 2, minGSSize = 3
  ))
  expect_true(all(c("ko", "Level3", "Level2", "Level1") %in% names(res@result)))
  expect_snapshot_value(pin_ora(res), style = "deparse", tolerance = 1e-8)
})
