# AC-001 .. AC-004 (SPEC-RICHR-006 REQ-001, IC-1): a table plot draws exactly
# the rows its selection keeps and, when that drops rows, signals exactly one
# `richR_plot_selection` message per call naming the counts and the rule.

# The richR_plot_selection conditions signalled while evaluating `expr`, and
# the value; every message is muffled.
.m61_sel <- function(expr) {
  conds <- list()
  value <- withCallingHandlers(expr, message = function(m) {
    if (inherits(m, "richR_plot_selection")) conds[[length(conds) + 1L]] <<- m
    invokeRestart("muffleMessage")
  })
  list(value = value, conds = conds)
}

# Drawn identifiers, read from plot data and built layer data ("term A1" ...).
.m61_terms <- function(p) {
  b <- ggplot2::ggplot_build(p)$data
  tm <- c(if (is.data.frame(p$data) && "Term" %in% names(p$data)) as.character(p$data$Term),
          unlist(lapply(b, function(l) if ("label" %in% names(l)) as.character(l$label))))
  sort(unique(grep("^term ", tm, value = TRUE)))
}

.m61_expect_one <- function(s, fun, has, info) {
  expect_length(s$conds, 1L)
  if (length(s$conds) != 1L) return(invisible())
  cond <- s$conds[[1]]
  expect_identical(class(cond), c("richR_plot_selection", "message", "condition"), info = info)
  msg <- conditionMessage(cond)
  expect_true(startsWith(msg, paste0(fun, ": drawing ")), info = paste(info, msg))
  for (h in has) expect_true(grepl(h, msg, fixed = TRUE), info = paste(info, "missing", h, "in", msg))
}

test_that("single-result plots draw the rows passing the defaults and say how many were dropped (AC-001)", {
  # Catches: a plot drawing fewer terms than its input holds without telling the
  # caller (R-1), a count or rule that misstates the selection, a message when
  # nothing was dropped, and any change to which rows the defaults select.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- .m61_fp6(); d <- x@result
  expect_identical(d$Annot[d$Pvalue <= 0.05 & d$Padj <= 0.05], c("A1", "A2", "A3"))
  three <- paste("term", c("A1", "A2", "A3")); six <- paste("term", paste0("A", 1:6))
  funs <- c("richBar", "richDot", "richNetplot", "richNetwork", "richLollipop", "richCircle",
            "richScatter", "richTermSim", "richGeneDot", "richGeneHeat", "richGeneBar")
  for (fun in funs) {
    f <- get(fun, envir = asNamespace("richR"))
    for (input in c("richResult", "data.frame")) {
      obj <- if (input == "richResult") x else d
      extra <- if (fun == "richNetwork" && input == "data.frame") list(gene = x@gene) else list()
      info <- paste(fun, input)
      set.seed(1)
      s <- .m61_sel(do.call(f, c(list(obj), extra)))
      expect_identical(.m61_terms(s$value), three, info = info)
      .m61_expect_one(s, fun, c("3 of 6 rows", "Pvalue <= 0.05", "Padj <= 0.05", "input order"), info)
      set.seed(1)
      s <- .m61_sel(do.call(f, c(list(obj, pvalue = 1, padj = NULL), extra)))
      expect_identical(.m61_terms(s$value), six, info = paste(info, "pvalue = 1"))
      expect_length(s$conds, 0L)
    }
  }
})

test_that("top truncation is reported with its order rule, and top = Inf draws every row (AC-002)", {
  # Catches: positional and directional truncation drawn silently (R-7), a
  # heat-map count in rows instead of identifiers, and an order rule other
  # than the documented one.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- .m61_fp6(); r2 <- .m61_r2()
  s <- .m61_sel(richBar(x, pvalue = 1, padj = NULL, top = 2))
  expect_identical(as.character(s$value$data$Annot), c("A1", "A2"))
  .m61_expect_one(s, "richBar", c("2 of 6 rows", "input order"), "richBar top = 2")
  s <- .m61_sel(richBar(x, pvalue = 1, padj = NULL, top = Inf))
  expect_identical(nrow(s$value$data), 6L)
  expect_length(s$conds, 0L)

  s <- .m61_sel(richVolcano(x, top = 2))
  expect_setequal(as.character(s$value$data$Annot), c("A1", "A2"))
  .m61_expect_one(s, "richVolcano", c("2 of 6 rows", "input order"), "richVolcano top = 2")
  s <- .m61_sel(richVolcano(x, top = Inf))
  expect_identical(nrow(s$value$data), 6L)
  expect_length(s$conds, 0L)

  s <- .m61_sel(richHeatmap(list(x, r2), top = 2))
  expect_identical(sort(unique(as.character(s$value$data$Term))), paste("term", c("A1", "A2", "B9")))
  .m61_expect_one(s, "richHeatmap", c("3 of 7 identifiers"), "richHeatmap top = 2")
  s <- .m61_sel(richHeatmap(list(x, r2), top = Inf))
  expect_length(unique(as.character(s$value$data$Term)), 7L)
  expect_length(s$conds, 0L)
})

test_that("richNES and richECDF report their directional top selection (AC-002)", {
  # Catches: only the most positive (richNES) or most negative (richECDF) NES
  # drawn without notice (R-7).
  .char_skip_unless_recorded_fgsea()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  g <- .m61_gsea(testthat::test_path("fixtures"))
  expect_identical(nrow(g@result), 5L)
  s <- .m61_sel(richNES(g, pvalue = 1, padj = 1, top = 2))
  expect_equal(round(s$value$data$NES, 2), c(2.38, 2.17))
  .m61_expect_one(s, "richNES", c("2 of 5 rows", "NES, descending"), "richNES top = 2")
  s <- .m61_sel(richNES(g, pvalue = 1, padj = 1, top = Inf))
  expect_identical(nrow(s$value$data), 5L)
  expect_length(s$conds, 0L)

  s <- .m61_sel(richECDF(g, pvalue = 1, padj = 1, top = 2))
  expect_equal(round(s$value$data$NES, 2), c(-2.4, -0.54))
  .m61_expect_one(s, "richECDF", c("2 of 5 rows", "NES, ascending"), "richECDF top = 2")
  s <- .m61_sel(richECDF(g, pvalue = 1, padj = 1, top = Inf))
  expect_identical(nrow(s$value$data), 5L)
  expect_length(s$conds, 0L)
})

test_that("richCompareDot reports re-thresholding a comparison table (AC-003)", {
  # Catches: a comparison table built with explicit cutoffs re-thresholded by
  # the plot defaults without notice (R-10).
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- .m61_fp6()
  cmp <- compareResult(list(a = x, b = x), pvalue = 1, padj = NULL)
  expect_identical(nrow(cmp), 12L)
  s <- .m61_sel(richCompareDot(cmp))
  expect_identical(nrow(s$value$data), 6L)
  .m61_expect_one(s, "richCompareDot", c("6 of 12 rows", "Pvalue <= 0.05", "Padj <= 0.05"),
                  "richCompareDot defaults")
  for (args in list(list(pvalue = 1, padj = NULL), list(include.all = TRUE))) {
    s <- .m61_sel(do.call(richCompareDot, c(list(cmp), args)))
    expect_identical(nrow(s$value$data), 12L)
    expect_length(s$conds, 0L)
  }
})

test_that("richNetmap counts reported rows, not padding, and signals once (AC-004)", {
  # Catches: padding rows (R-8, R-14: 100 rows received, 92 all-NA) or the
  # per-group top counted as input, and a second message from the nested
  # richNetwork call.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- .m61_fp6(); r2b <- .m61_r2b()
  set.seed(1)
  s <- .m61_sel(richNetmap(list(a = x, b = r2b)))
  expect_identical(.m61_terms(s$value), paste("term", c("A1", "A2", "A3", "B1", "B2")))
  .m61_expect_one(s, "richNetmap", c("5 of 8 rows"), "richNetmap defaults")
  if (length(s$conds)) {
    msg <- conditionMessage(s$conds[[1]])
    expect_false(grepl("of 100", msg, fixed = TRUE))
    expect_false(grepl("of 50 ", msg, fixed = TRUE))
  }
  set.seed(1)
  s <- .m61_sel(richNetmap(list(a = x, b = r2b), pvalue = 1, padj = NULL))
  expect_length(.m61_terms(s$value), 8L)
  expect_length(s$conds, 0L)
})
