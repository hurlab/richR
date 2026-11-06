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
