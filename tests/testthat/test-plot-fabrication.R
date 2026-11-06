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

test_that("a member gene without a fold change is drawn grey55, never as fold change 0 (AC-007)", {
  # Catches: a gene without a fold change drawn as fold change 0 - white, the
  # same as "no change" (R-5, R-17, R-23, R-24); such a gene dropped at render
  # time (an NA size, R-23); and "member without a fold change" drawn the same
  # as "not a member" (audit N3).
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- .m61_fp6(); fc5 <- .m61_fc5()
  grey55 <- .m61_rgb("grey55")
  expect_identical(grey55, .m61_rgb("#8C8C8C"))

  # dot leg: A1 = g01-g07, fc holds g01-g05
  pd <- richGeneDot(x, fc = fc5, pvalue = 1, padj = NULL, top = 1)
  d <- pd$data
  expect_setequal(d$Gene, sprintf("g%02d", 1:7))
  expect_identical(nrow(d), 7L)
  missing <- d$Gene %in% c("g06", "g07")
  expect_true(all(is.na(d$fc_value[missing])))
  expect_equal(unname(d$fc_value[!missing]), unname(fc5[d$Gene[!missing]]))
  pts <- do.call(rbind, lapply(.m61_layers_of(pd, "GeomPoint"), function(i) {
    ld <- ggplot2::layer_data(pd, i)
    data.frame(gene = .m61_axis_values(pd, i, "x"), size = ld$size, colour = ld$colour,
               stringsAsFactors = FALSE)
  }))
  expect_identical(nrow(pts), 7L)
  expect_true(all(is.finite(pts$size)))
  expect_setequal(pts$gene, sprintf("g%02d", 1:7))
  for (g in c("g06", "g07")) {
    expect_identical(.m61_rgb(pts$colour[pts$gene == g]), grey55, info = g)
  }
  expect_false(any(grepl("Removed [0-9]+ rows", .m61_render_warnings(pd))))
  cap <- ggplot2::get_labs(pd)$caption
  expect_true(is.character(cap) && grepl("grey55", cap, fixed = TRUE) &&
                grepl("no fold change supplied", cap, fixed = TRUE), info = cap)

  # heat leg: A1 = g01-g07, A2 = g02-g09
  ph <- richGeneHeat(x, fc = fc5, pvalue = 1, padj = NULL, top = 2)
  m <- ph$data
  expect_identical(nrow(m), 18L)
  key <- paste(m$Gene, sub("^term ", "", m$Term))
  nonmember <- is.na(m$present)
  member_missing <- !nonmember & !(m$Gene %in% names(fc5))
  with_fc <- !nonmember & !member_missing
  expect_setequal(key[nonmember], c("g08 A1", "g09 A1", "g01 A2"))
  expect_setequal(key[member_missing], c("g06 A1", "g07 A1", "g06 A2", "g07 A2", "g08 A2", "g09 A2"))
  expect_identical(sum(with_fc), 9L)
  expect_true(all(is.na(m$fc_value[member_missing])))
  # drawn fill of a cell: the fill of the last tile layer with a tile there
  drawn <- stats::setNames(rep(NA_character_, nrow(m)), key)
  for (i in .m61_layers_of(ph, "GeomTile")) {
    ld <- ggplot2::layer_data(ph, i)
    at <- paste(.m61_axis_values(ph, i, "x"), sub("^term ", "", .m61_axis_values(ph, i, "y")))
    drawn[at] <- ld$fill
  }
  expect_false(anyNA(drawn))
  rgb_of <- function(cells) lapply(drawn[key[cells]], .m61_rgb)
  for (v in rgb_of(nonmember)) expect_identical(v, .m61_rgb("grey90"))
  for (v in rgb_of(member_missing)) expect_identical(v, grey55)
  expect_false(identical(.m61_rgb("grey90"), grey55))
  expect_false(any(vapply(rgb_of(with_fc), identical, logical(1), grey55)))
  expect_false(any(grepl("Removed [0-9]+ rows", .m61_render_warnings(ph))))
  cap <- ggplot2::get_labs(ph)$caption
  expect_true(is.character(cap) && grepl("grey55", cap, fixed = TRUE) &&
                grepl("no fold change supplied", cap, fixed = TRUE) &&
                grepl("not in the term", cap, fixed = TRUE), info = cap)
})

test_that("richHeatmap draws an absent term as missing and names the column it draws (AC-008)", {
  # Catches: an absent term shown as p = 1 (R-6), a legend naming the other
  # p-value column (R-13), and a row-order change.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- .m61_fp6(); r2 <- .m61_r2()
  tab <- richHeatmap(list(x, r2), returnData = TRUE)
  cell <- function(term, group) tab[tab$Term == paste("term", term), group]
  expect_true(is.na(cell("A2", "Group_2")))
  expect_true(is.na(cell("B9", "Group_1")))
  expect_identical(cell("A1", "Group_1"), 1e-4)
  expect_identical(cell("A1", "Group_2"), 0.02)
  expect_identical(cell("A2", "Group_1"), 1e-3)
  expect_identical(cell("B9", "Group_2"), 0.03)

  ref <- readRDS(testthat::test_path("fixtures", "plot_reference_m6base.rds"))$plots
  for (use in c(FALSE, TRUE)) {
    col <- if (use) "Padj" else "Pvalue"
    p <- richHeatmap(list(x, r2), usePadj = use)
    info <- paste("usePadj =", use)
    expect_identical(ggplot2::get_labs(p)$fill, paste0("-log10(", col, ")"), info = info)
    d <- p$data
    expect_true(col %in% names(d), info = info)
    drawn <- stats::setNames(d[[col]], paste(sub("^term ", "", d$Term), d$variable))
    expect_true(is.na(drawn[["A2 Group_2"]]), info = info)
    expect_true(is.na(drawn[["B9 Group_1"]]), info = info)
    expect_identical(drawn[["A1 Group_1"]], x@result[["A1", col]], info = info)
    expect_identical(drawn[["B9 Group_2"]], r2@result[["B9", col]], info = info)
    # absent cells are drawn in the scale's missing-value colour
    b <- ggplot2::ggplot_build(p)
    na_col <- b$plot$scales$get_scales("fill")$na.value
    ld <- b$data[[1]]
    at <- paste(sub("^term ", "", .m61_axis_values(p, 1, "y")), .m61_axis_values(p, 1, "x"))
    # A2-A6 are absent from r2, B9 from F-P6
    absent <- c(paste0("A", 2:6, " Group_2"), "B9 Group_1")
    expect_true(all(absent %in% at), info = info)
    for (k in absent) {
      expect_true(is.na(drawn[[k]]), info = paste(info, k))
      expect_identical(.m61_rgb(ld$fill[at == k]), .m61_rgb(na_col), info = paste(info, k))
    }
    for (k in setdiff(at, absent)) {
      expect_false(identical(.m61_rgb(ld$fill[at == k]), .m61_rgb(na_col)), info = paste(info, k))
    }
    # row order as on the M6.1 base
    r <- ref[[if (use) "richHeatmap_ac008_q" else "richHeatmap_ac008_p"]]
    expect_identical(levels(d$Term), levels(r$data$Term), info = info)
  }
})
