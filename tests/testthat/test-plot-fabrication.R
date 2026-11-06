# AC-005 .. AC-008 (SPEC-RICHR-006 REQ-002 .. REQ-005): plots draw only values
# the result carries, titled by what they are.

# Discrete-axis label of each row of a layer's data (the scale's limits, in
# level order, indexed by the mapped position).
.m61_axis_values <- function(p, layer, aes = "x") {
  b <- ggplot2::ggplot_build(p)
  sc <- if (aes == "x") b$layout$panel_scales_x[[1]] else b$layout$panel_scales_y[[1]]
  sc$get_limits()[as.integer(b$data[[layer]][[aes]])]
}

# Index of the layers drawn with a given geom.
.m61_layers_of <- function(p, geom) {
  which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
}

# Warnings signalled while rendering a plot.
.m61_render_warnings <- function(p) {
  w <- character()
  withCallingHandlers(ggplot2::ggplotGrob(p), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  w
}

.m61_rgb <- function(col) grDevices::col2rgb(col)[, 1]

test_that("GSEA lollipop, gene-bar and scatter plots draw NES, not a constant (AC-005)", {
  # Catches: GSEA pathways drawn as a constant RichFactor = 1 or log2 FE = 0
  # under an ORA title (R-3, R-16).
  skip_if_not_installed("fgsea")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  g <- .m61_gsea(testthat::test_path("fixtures"))
  res <- g@result
  nes <- stats::setNames(res$NES, res$pathway)
  size <- stats::setNames(res$size, res$pathway)
  expect_gt(length(unique(round(nes, 6))), 1L)

  for (fun in c("richLollipop", "richGeneBar")) {
    p <- suppressMessages(get(fun)(g, pvalue = 1, padj = 1))
    layer <- .m61_layers_of(p, if (fun == "richLollipop") "GeomPoint" else "GeomBar")[1]
    term <- .m61_axis_values(p, layer, "x")
    pathway <- p$data$Annot[match(term, p$data$Term)]
    expect_setequal(pathway, names(nes))
    expect_equal(ggplot2::ggplot_build(p)$data[[layer]]$y, unname(nes[pathway]),
                 tolerance = 1e-12, info = fun)
    expect_identical(ggplot2::get_labs(p)$y, "NES", info = fun)
    expect_false(any(c("RichFactor", "Annotated") %in% names(p$data)), info = fun)
  }

  p <- suppressMessages(richScatter(g, pvalue = 1, padj = 1))
  ld <- ggplot2::layer_data(p, .m61_layers_of(p, "GeomPoint")[1])
  pathway <- p$data$Annot
  expect_setequal(pathway, names(nes))
  expect_equal(ld$y, unname(nes[pathway]), tolerance = 1e-12)
  expect_equal(ld$x, unname(log10(size[pathway])), tolerance = 1e-12)
  expect_identical(ggplot2::get_labs(p)$y, "NES")
  expect_identical(ggplot2::get_labs(p)$x, "log10(Gene-set size)")
  expect_false(any(c("RichFactor", "Annotated", "FoldEnrichment", "log2_fe") %in% names(p$data)))
})

test_that("volcano and scatter titles name the effect they draw when FoldEnrichment is absent (AC-006)", {
  # Catches: a rich factor or a constant presented as fold enrichment (R-4,
  # R-15), and results saved before richR 0.1.2 (no FoldEnrichment) broken by
  # an error.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  with_fe <- .m61_fp6()@result
  no_fe <- with_fe; no_fe$FoldEnrichment <- NULL
  none <- no_fe; none$Significant <- NULL; none$Annotated <- NULL
  rich <- stats::setNames(log2(with_fe$Significant / with_fe$Annotated), with_fe$Annot)
  expect_equal(round(unname(rich), 3), c(-1.737, -1.515, -1.322, -1.152, -1, -0.862))

  v <- richVolcano(with_fe, pvalue = 1, padj = NULL)
  expect_identical(ggplot2::get_labs(v)$x, "log2(Fold Enrichment)")
  expect_identical(ggplot2::get_labs(v)$colour, "log2(FE)")
  s <- richScatter(with_fe, pvalue = 1, padj = NULL)
  expect_identical(ggplot2::get_labs(s)$y, "log2(Fold Enrichment)")

  v <- richVolcano(no_fe, pvalue = 1, padj = NULL)
  expect_equal(stats::setNames(v$data$effect, v$data$Annot)[names(rich)], rich, tolerance = 1e-12)
  expect_identical(ggplot2::get_labs(v)$x, "log2(Rich Factor = Significant/Annotated)")
  expect_identical(ggplot2::get_labs(v)$colour, "log2(Rich Factor)")
  s <- richScatter(no_fe, pvalue = 1, padj = NULL)
  ld <- ggplot2::layer_data(s, .m61_layers_of(s, "GeomPoint")[1])
  expect_equal(stats::setNames(ld$y, s$data$Annot)[names(rich)], rich, tolerance = 1e-12)
  expect_identical(ggplot2::get_labs(s)$y, "log2(Rich Factor = Significant/Annotated)")

  v <- richVolcano(none, pvalue = 1, padj = NULL)
  expect_true(all(v$data$effect == 0))
  expect_identical(nrow(v$data), 6L)
  expect_identical(ggplot2::get_labs(v)$x, "No effect column: all points at 0")
  expect_identical(ggplot2::get_labs(v)$colour, "No effect column")
})
