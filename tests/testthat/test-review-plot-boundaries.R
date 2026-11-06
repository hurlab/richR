test_that("heatmap raw data and rendering tolerate undefined clustering", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  x <- .m61_fp3()
  one <- x; one@result <- one@result[1, , drop = FALSE]
  flat <- x; flat@result$Pvalue <- rep(0.02, 3)
  cases <- list(list(one, one), list(x), list(flat, flat), list(x, x))
  for (inputs in cases) {
    tab <- richHeatmap(inputs, returnData = TRUE)
    expect_equal(tab$Group_1, inputs[[1]]@result$Pvalue)
    p <- richHeatmap(inputs)
    expect_s3_class(p, "ggplot")
    expect_equal(sort(p$data$Pvalue), sort(unlist(tab[-1], use.names = FALSE)))
    expect_silent(ggplot2::ggplotGrob(p))
    if (length(inputs) == 1L || identical(inputs[[1]], flat) || nrow(tab) == 1L)
      expect_identical(levels(p$data$Term), tab$Term)
  }
})

test_that("genuine empty dot results signal no result and return invisible NULL", {
  x <- .m61_fp3(); x@result <- x@result[FALSE, , drop = FALSE]
  out <- NULL
  expect_output(out <- withVisible(richDot(x)), "No Pathway enrichment results were found!")
  expect_null(out$value)
  expect_false(out$visible)
})
