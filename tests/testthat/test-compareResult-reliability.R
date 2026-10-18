batch_test_annotation <- function() {
  data.frame(
    GeneID = LETTERS[1:6],
    Annot = rep(c("T1", "T2", "T3"), each = 2),
    Term = rep(c("Term one", "Term two", "Term three"), each = 2)
  )
}

batch_test_result <- function(genes, cutoff = 0.1) {
  suppressMessages(richR::enrich(
    genes, batch_test_annotation(), pvalue = cutoff, minSize = 1,
    minGSSize = 1, keepRich = FALSE
  ))
}

test_that("compareResult preserves unnamed and distinct explicit group identities", {
  first <- batch_test_result(c("A", "B"))
  second <- batch_test_result(c("C", "D"))

  unnamed <- richR::compareResult(list(first, second), pvalue = 0.1)
  expect_identical(unique(unnamed$group), c("Group 1", "Group 2"))
  expect_identical(unnamed$Annot, c("T1", "T2"))

  expect_error(
    richR::compareResult(setNames(list(first, second), c("same", "same"))),
    "same"
  )
  expect_error(richR::compareResult(setNames(list(first, second), c("ok", NA))),
               "missing|position")
  expect_error(richR::compareResult(setNames(list(first, second), c("ok", "  "))),
               "blank|position")
  expect_error(richR::compareResult(list()), "non-empty")
  expect_error(richR::compareResult("not a list"), "list")
  expect_error(richR::compareResult(list(bad = 42)), "bad.*richResult.*data.frame.*NULL")
})

test_that("compareResult keeps its strict cutoff policy and real extra columns", {
  first <- batch_test_result(c("A", "B"))
  second <- batch_test_result(c("C", "D"))
  first_df <- first@result
  first_df$extra_first <- "first"
  second_df <- second@result
  second_df$extra_second <- 7L

  strict <- richR::compareResult(list(first = first, second = second))
  expect_equal(nrow(strict), 0L)
  retained <- richR::compareResult(list(first = first, second = second), pvalue = 0.1)
  expect_equal(retained$Pvalue, c(1 / 15, 1 / 15), tolerance = 1e-12)
  expect_identical(retained$group, c("first", "second"))

  with_extra <- richR::compareResult(
    list(first = first_df, second = second_df), pvalue = 0.1
  )
  expect_identical(with_extra$extra_first, c("first", NA_character_))
  expect_identical(with_extra$extra_second, c(NA_integer_, 7L))
})

test_that("include-all comparison preserves group axes without synthetic terms", {
  first <- batch_test_result(c("A", "B"))
  second <- batch_test_result(c("C", "D"))
  empty <- batch_test_result(c("A", "B"), cutoff = 0.05)

  all_nonempty <- richR::compareResult(
    list(first = first, second = second), pvalue = 0.1, include.all = TRUE
  )
  expect_identical(levels(all_nonempty$group), c("first", "second"))
  expect_equal(nrow(all_nonempty), nrow(first@result) + nrow(second@result))
  expect_false(any(all_nonempty$Pvalue == 1))

  mixed <- richR::compareResult(
    list(first = first, empty = empty, missing = NULL), pvalue = 0.1,
    include.all = TRUE
  )
  expect_identical(levels(mixed$group), c("first", "empty", "missing"))
  expect_identical(unique(as.character(mixed$group)), "first")
  expect_equal(nrow(mixed), nrow(first@result))
  expect_type(mixed$group, "integer")

  all_empty <- richR::compareResult(
    list(empty = empty, missing = NULL), include.all = TRUE
  )
  expect_equal(nrow(all_empty), 0L)
  expect_identical(levels(all_empty$group), c("empty", "missing"))
  expect_type(all_empty$Annot, "character")
  expect_identical(typeof(all_empty$Annotated), typeof(empty@result$Annotated))
})

test_that("include-all comparison preserves standalone significance cutoffs", {
  first <- batch_test_result(c("A", "B"))
  result <- richR::compareResult(list(first = first), include.all = TRUE)
  expect_equal(nrow(result), 0L)
  expect_identical(levels(result$group), "first")
})

test_that("include-all plots retain empty axes and draw only real results", {
  first <- batch_test_result(c("A", "B"))
  empty <- batch_test_result(c("A", "B"), cutoff = 0.05)
  compared <- richR::compareResult(
    list(first = first, empty = empty, missing = NULL), pvalue = 0.1,
    include.all = TRUE
  )

  plot <- richR::richCompareDot(compared, include.all = TRUE)
  built <- ggplot2::ggplot_build(plot)
  labels <- built$layout$panel_params[[1]]$x$get_labels()
  expect_identical(labels, c("first", "empty", "missing"))
  expect_equal(nrow(built$data[[1]]), nrow(first@result))

  filtered_plot <- richR::richCompareDot(
    compared, include.all = TRUE, pathway = "T1", level = "Term"
  )
  filtered_labels <- ggplot2::ggplot_build(filtered_plot)$layout$panel_params[[1]]$x$get_labels()
  expect_identical(filtered_labels, c("first", "empty", "missing"))

  dropped_plot <- richR::richCompareDot(
    compared, pvalue = 0.1, pathway = "T1", level = "Term"
  )
  dropped_build <- ggplot2::ggplot_build(dropped_plot)
  expect_identical(dropped_build$layout$panel_params[[1]]$x$get_labels(), "first")
  expect_equal(nrow(dropped_build$data[[1]]), 1L)

  all_empty <- richR::compareResult(list(empty = empty, missing = NULL), include.all = TRUE)
  empty_plot <- richR::richCompareDot(all_empty, include.all = TRUE)
  empty_build <- ggplot2::ggplot_build(empty_plot)
  expect_identical(empty_build$layout$panel_params[[1]]$x$get_labels(),
                   c("empty", "missing"))
  expect_equal(nrow(empty_build$data[[1]]), 0L)

  output <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(output, plot = empty_plot, width = 5, height = 4)
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
})

test_that("compareResult uses strict pvalue and adjusted-pvalue cutoffs", {
  boundary <- data.frame(
    Annot = c("p-boundary", "padj-boundary", "inside"),
    Term = c("P", "Q", "R"),
    Annotated = c(2, 2, 2), Significant = c(1, 1, 1),
    Pvalue = c(0.05, 0.01, 0.049), Padj = c(0.01, 1, 0.99),
    RichFactor = c(0.5, 0.5, 0.5), GeneID = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  pfiltered <- richR::compareResult(list(g = boundary), pvalue = 0.05)
  expect_identical(pfiltered$Annot, "inside")
  padj_filtered <- richR::compareResult(list(g = boundary), padj = 1)
  expect_identical(padj_filtered$Annot, c("p-boundary", "inside"))
})

test_that("all NULL comparisons return a typed empty frame and retain plot labels", {
  all_null <- richR::compareResult(list(first = NULL, second = NULL), include.all = TRUE)
  expect_equal(nrow(all_null), 0L)
  expect_identical(levels(all_null$group), c("first", "second"))
  expect_type(all_null$Pvalue, "double")
  expect_type(all_null$Annot, "character")
  plot <- richR::richCompareDot(all_null, include.all = TRUE)
  expect_identical(ggplot2::ggplot_build(plot)$layout$panel_params[[1]]$x$get_labels(),
                   c("first", "second"))
})
