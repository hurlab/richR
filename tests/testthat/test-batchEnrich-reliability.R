batch_test_annotation <- function() {
  data.frame(
    GeneID = LETTERS[1:6],
    Annot = rep(c("T1", "T2", "T3"), each = 2),
    Term = rep(c("Term one", "Term two", "Term three"), each = 2)
  )
}

test_that("batchEnrich combines exact returned rows and reports every group outcome", {
  annot <- batch_test_annotation()
  groups <- list(success = c("A", "B"), empty = c("A", "C"), failed = character())
  result <- expect_warning(suppressMessages(richR::batchEnrich(
    groups, annot, pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE, verbose = FALSE
  )), "failed for 'failed'")

  expect_named(result, c("results", "combined", "status"))
  expect_named(result$results, c("success", "empty"))
  expect_equal(nrow(result$results$success@result), 1L)
  expect_equal(nrow(result$results$empty@result), 0L)
  expect_equal(nrow(result$combined), sum(vapply(result$results, function(x) nrow(x@result), integer(1))))
  expect_identical(levels(result$combined$group), c("success", "empty"))
  expect_equal(result$combined$Pvalue, result$results$success@result$Pvalue)
  expect_true(all(result$combined$Pvalue > 0.05))

  status <- result$status
  expect_identical(names(status), c("group", "status", "n_terms", "error"))
  expect_identical(status$group, c("success", "empty", "failed"))
  expect_identical(status$status, c("success", "empty", "failed"))
  expect_type(status$group, "character")
  expect_type(status$status, "character")
  expect_type(status$n_terms, "integer")
  expect_type(status$error, "character")
  expect_identical(status$n_terms, c(1L, 0L, NA_integer_))
  expect_identical(status$error[1:2], c(NA_character_, NA_character_))
  expect_match(status$error[[3]], "empty")
})

test_that("batchEnrich returns valid all-empty reports and errors when all fail", {
  annot <- batch_test_annotation()
  empty <- suppressMessages(richR::batchEnrich(
    list(first = "A", second = "C"), annot, pvalue = 0.05,
    minSize = 1, minGSSize = 1, keepRich = FALSE, verbose = FALSE
  ))

  expect_named(empty$results, c("first", "second"))
  expect_true(all(vapply(empty$results, function(x) nrow(x@result) == 0L, logical(1))))
  expect_equal(nrow(empty$combined), 0L)
  expect_identical(levels(empty$combined$group), c("first", "second"))
  expect_identical(empty$status$group, c("first", "second"))
  expect_identical(empty$status$status, c("empty", "empty"))
  expect_identical(empty$status$n_terms, c(0L, 0L))
  expect_identical(empty$status$error, c(NA_character_, NA_character_))

  expect_warning(expect_error(richR::batchEnrich(
    list(first = character(), second = "unmapped"), annot, verbose = FALSE
  ), "all enrichment analyses failed"), "failed for 'first'")
})

test_that("batchEnrich reports a group-specific invalid-input failure", {
  result <- expect_warning(suppressMessages(richR::batchEnrich(
    list(valid = c("A", "B"), invalid = NA_character_),
    batch_test_annotation(), pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE, verbose = FALSE
  )), "failed for 'invalid'")
  expect_identical(names(result$results), "valid")
  expect_identical(result$status$status, c("success", "failed"))
  expect_match(result$status$error[[2]], "NA|missing|only")
})

test_that("batchEnrich combines exact p-value and adjusted-p results without a second cutoff", {
  annot <- batch_test_annotation()
  groups <- list(first = c("A", "B"), second = c("C", "D"))
  pvalue_batch <- suppressMessages(richR::batchEnrich(
    groups, annot, pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE, verbose = FALSE
  ))
  expect_equal(vapply(pvalue_batch$results, function(x) nrow(x@result), integer(1)),
               c(first = 1L, second = 1L))
  expect_equal(vapply(pvalue_batch$results, function(x) x@result$Pvalue[[1]], numeric(1)),
               c(first = 1 / 15, second = 1 / 15), tolerance = 1e-12)
  expect_equal(nrow(pvalue_batch$combined), 2L)
  expect_identical(as.character(pvalue_batch$combined$group), c("first", "second"))
  expect_equal(pvalue_batch$combined$Pvalue, c(1 / 15, 1 / 15), tolerance = 1e-12)
  standalone <- richR::compareResult(pvalue_batch$results)
  expect_equal(nrow(standalone), 0L)

  padj_batch <- suppressMessages(richR::batchEnrich(
    groups, annot, padj = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE, verbose = FALSE
  ))
  expect_equal(nrow(padj_batch$combined), 2L)
  expect_equal(padj_batch$combined$Padj, c(1 / 15, 1 / 15), tolerance = 1e-12)
})
