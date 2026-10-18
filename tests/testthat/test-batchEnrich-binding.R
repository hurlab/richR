batch_binding_annotation <- function() {
  data.frame(
    GeneID = LETTERS[1:6],
    Annot = rep(c("T1", "T2", "T3"), each = 2),
    Term = rep(c("Term one", "Term two", "Term three"), each = 2)
  )
}

test_that("batchEnrich binds differing columns from real native ORA results", {
  annot <- batch_binding_annotation()
  first <- suppressMessages(richR::enrich(
    c("A", "B"), annot, pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  ))
  second <- suppressMessages(richR::enrich(
    c("C", "D"), annot, pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  ))
  first@result$extra_first <- "first"
  second@result$extra_second <- 7L

  batched <- NULL
  testthat::with_mocked_bindings(
    enrich = function(x, object, ...) {
      if (identical(x, c("A", "B"))) first else second
    },
    {
      batched <- suppressMessages(richR::batchEnrich(
        list(first = c("A", "B"), second = c("C", "D")), annot,
        pvalue = 0.1, minSize = 1, minGSSize = 1,
        keepRich = FALSE, verbose = FALSE
      ))
    },
    .package = "richR"
  )

  expect_equal(batched$combined$Pvalue, c(first@result$Pvalue, second@result$Pvalue))
  expect_identical(batched$combined$extra_first, c("first", NA_character_))
  expect_identical(batched$combined$extra_second, c(NA_integer_, 7L))
})
