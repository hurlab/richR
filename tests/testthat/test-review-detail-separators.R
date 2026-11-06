# SPEC-RICHR-010: separator metacharacters are literal gene delimiters.
.review010_annot <- function() {
  data.frame(GeneID = paste0("G", 1:60),
             Term = rep(paste0("T", 1:6), each = 10))
}
.review010_enrich <- function(x, sep = ",", annot = .review010_annot(), ...) {
  suppressMessages(enrich(x, annot, pvalue = 1, padj = 1,
                         minSize = 1, minGSSize = 1, sep = sep, ...))
}

test_that("literal separators preserve ORA detail identities and score joins", {
  genes <- paste0("G", c(1, 2, 11, 12))
  scored <- data.frame(score = 1:4, row.names = genes)
  control <- .review010_enrich(scored)
  for (sep in c(",", "|", ".")) {
    vector <- .review010_enrich(genes, sep)
    result <- .review010_enrich(scored, sep)
    expect_setequal(vector@detail$GeneID, genes)
    expect_setequal(result@detail$GeneID, genes)
    expect_equal(result@detail$score, scored[result@detail$GeneID, "score"], info = sep)
    # Only the serialized separator and its recorded provenance differ.
    expect_identical(result@result[setdiff(names(result@result), "GeneID")],
                     control@result[setdiff(names(control@result), "GeneID")])
    expect_identical(result@tested[setdiff(names(result@tested), "GeneID")],
                     control@tested[setdiff(names(control@tested), "GeneID")])
    expect_identical(strsplit(result@tested$GeneID, sep, fixed = TRUE),
                     strsplit(control@tested$GeneID, ",", fixed = TRUE))
    expect_identical(result@gene, control@gene)
    expect_identical(result@sep, sep)
  }
})

test_that("analysis row names control detail scores without mutating inputs or standalone joins", {
  genes <- paste0("G", c(1, 2, 11, 12))
  scored <- data.frame(gene = rev(genes), score = 1:4, row.names = genes)
  original <- scored
  control_input <- scored["score"]
  control <- .review010_enrich(control_input)
  for (sep in c(",", "|", ".")) {
    result <- .review010_enrich(scored, sep)
    expect_identical(scored, original)
    expect_equal(result@detail$score, scored[result@detail$GeneID, "score"])
    expect_equal(result@detail, control@detail)
    standalone <- getdetail(result@result, scored, sep = sep)
    expect_equal(standalone$score, scored$score[match(standalone$GeneID, scored$gene)])
  }
  # A supplied universe remains the analysis population.
  universe <- paste0("G", 1:30)
  result <- .review010_enrich(scored, universe = universe)
  expected <- .review010_enrich(control_input, universe = universe)
  expect_identical(result@result, expected@result)
  expect_identical(result@tested, expected@tested)
  expect_identical(enrichmentInfo(result), enrichmentInfo(expected))
})

test_that("gene plots keep complete IDs for literal separators", {
  genes <- paste0("G", c(1, 2, 11, 12))
  for (sep in c(",", "|", ".")) {
    result <- .review010_enrich(genes, sep)
    bar <- suppressMessages(richGeneBar(result, pvalue = 1, padj = NULL))
    dot <- suppressMessages(richGeneDot(result, pvalue = 1, padj = NULL))
    expect_setequal(unlist(strsplit(bar$data$gene_label, ", ", fixed = TRUE)), genes)
    expect_setequal(dot$data$Gene, genes)
  }
})

.review010_overlap <- function(sep) {
  annot <- rbind(
    data.frame(GeneID = paste0("G", c(1, 2, 3, 2, 3, 4, 6, 3, 5)),
               Term = rep(c("A", "B", "C"), c(3, 4, 2))),
    data.frame(GeneID = paste0("G", 6:20), Term = "BG"))
  .review010_enrich(paste0("G", 1:5), sep, annot)
}

for (.review010_fn in c("richNetwork", "richNetplot")) {
test_that(paste(.review010_fn, "preserves memberships for literal separators"), {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  control <- .review010_overlap(",")
  fn <- .review010_fn
  {
    draw <- function(x) {
      set.seed(10)
      suppressWarnings(suppressMessages(get(fn)(x, pvalue = 1, padj = NULL)))
    }
    expected <- draw(control)
    for (sep in c("|", ".")) {
      actual <- draw(.review010_overlap(sep))
      expect_equal(actual$data, expected$data, info = paste(fn, sep))
    }
    if (fn == "richNetplot") {
      expect_setequal(expected$data$label, c(paste0("G", 1:5), "A", "B", "C"))
    }
  }
})
}

test_that("TermSim similarities describe the original gene sets", {
  weights <- function(sep) {
    p <- suppressMessages(richTermSim(.review010_overlap(sep), pvalue = 1, padj = NULL))
    unlist(lapply(p$layers, function(layer) {
      if ("similarity" %in% names(layer$data)) layer$data$similarity else numeric()
    }))
  }
  expected <- weights(",")
  expect_equal(sort(unname(expected)), c(0.25, 0.25, 0.5))
  for (sep in c("|", ".")) expect_equal(weights(sep), expected)
})

test_that("dataframe UpSet keeps literal gene memberships and intersections", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # Capture the public plot's component data instead of rendering a raster.
  draw <- function(sep, column) {
    frames <- lapply(c(paste("G1", "G2", sep = sep), paste("G2", "G3", sep = sep)),
                     function(ids) stats::setNames(data.frame(ids), column))
    names(frames) <- c("a", "b")
    parts <- testthat::with_mocked_bindings(
      richUpset(frames, sep = sep, set.size.show = FALSE),
      plot_grid = function(...) list(...), .package = "richR")
    lapply(parts, function(part) if (inherits(part, "ggplot")) part$data else part)
  }
  for (column in c("GeneID", "leadingEdge")) {
    expected <- draw(",", column)
    for (sep in c("|", ".")) expect_equal(draw(sep, column), expected)
    extracted <- richR:::.extract_gene_lists(
      list(a = data.frame(GeneID = "G1|G2"), b = data.frame(GeneID = "G2|G3")), "|")
    expect_identical(extracted, list(a = c("G1", "G2"), b = c("G2", "G3")))
  }
})
