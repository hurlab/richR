# Fixtures for the plot tests of SPEC-RICHR-006 (M6.1, pure plots):
# test-plot-selection.R, test-plot-fabrication.R and test-plot-reference.R.
# fixtures/make-plot-reference.R sources this file too, so the AC-009 reference
# and the test that compares against it run the same calls on the same inputs.
# Inputs are built with methods::new() and the public API only.

# F-P6 builder: an ORA result with one row per identifier. Term `id` has the
# genes g(i)..g(2i+5): A1 = g01-g07, A2 = g02-g09, A3 = g03-g11.
.m61_mk <- function(ids, p, padj) {
  n <- length(ids)
  methods::new("richResult",
    result = data.frame(Annot = ids, Term = paste("term", ids), Annotated = rep(20, n),
      Significant = 5 + seq_len(n), RichFactor = (5 + seq_len(n)) / 20,
      FoldEnrichment = 2 + seq_len(n) / 10, zscore = 3, Pvalue = p, Padj = padj,
      GeneID = vapply(seq_len(n), function(i) paste(sprintf("g%02d", seq(i, i + 5 + i)), collapse = ","), ""),
      Level2 = "L2", stringsAsFactors = FALSE, row.names = ids),
    detail = data.frame(), pvalueCutoff = 1, pAdjustMethod = "BH", padjCutoff = 1,
    genenumber = 30, organism = "custom", ontology = "", gene = sprintf("g%02d", 1:30),
    keytype = "", sep = ",")
}

# F-P6: six rows; A1-A3 pass the default cutoffs (Pvalue <= 0.05, Padj <= 0.05).
.m61_fp6 <- function() {
  .m61_mk(paste0("A", 1:6), p = c(1e-4, 1e-3, 0.01, 0.04, 0.2, 0.4),
          padj = c(1e-3, 0.01, 0.03, 0.08, 0.4, 0.6))
}

# F-P6 restricted to the rows that pass the defaults (A1-A3).
.m61_fp3 <- function() {
  x <- .m61_fp6()
  x@result <- x@result[1:3, ]
  x
}

# r2: one identifier shared with F-P6 (A1), one not (B9).
.m61_r2 <- function() .m61_mk(c("A1", "B9"), p = c(0.02, 0.03), padj = c(0.04, 0.05))

# r2b: identifiers disjoint from F-P6, genes overlapping it; both rows pass.
.m61_r2b <- function() .m61_mk(c("B1", "B2"), p = c(0.002, 0.003), padj = c(0.004, 0.005))

# A second result over A1-A3 with different p-values (a two-result comparison
# whose rows all pass the defaults and whose heat-map rows vary).
.m61_fp3b <- function() {
  .m61_mk(paste0("A", 1:3), p = c(0.03, 2e-3, 1e-5), padj = c(0.04, 0.01, 1e-4))
}

# Fold changes for every gene g01-g30, and fc5: for g01-g05 only.
.m61_fc30 <- function() stats::setNames(seq(2, -2, length.out = 30), sprintf("g%02d", 1:30))
.m61_fc5 <- function() .m61_fc30()[1:5]

# The ranked list of F-GSEA, from the committed characterization fixture.
.m61_ranks <- function(dir) {
  rk <- utils::read.csv(file.path(dir, "gsea_ranked_a.csv"), stringsAsFactors = FALSE)
  stats::setNames(rk$score, rk$GeneID)
}

# F-GSEA: 5 pathways (needs fgsea; callers guard it).
.m61_gsea <- function(dir) {
  ann <- utils::read.csv(file.path(dir, "gsea_annot.csv"), stringsAsFactors = FALSE)
  out <- NULL
  utils::capture.output(out <- suppressWarnings(suppressMessages(
    richGSEA(.m61_ranks(dir), ann, pvalue = 1, padj = 1, minSize = 3, maxSize = 50, seed = 1))))
  out
}

# F-GSEA restricted to the pathways that pass the default cutoffs.
.m61_gsea_sig <- function(dir) {
  g <- .m61_gsea(dir)
  g@result <- g@result[!is.na(g@result$pval) & g@result$pval <= 0.05 &
                       !is.na(g@result$padj) & g@result$padj <= 0.05, ]
  g
}

# Plot data, layer data and titles of a ggplot; a data.frame (returnData = TRUE)
# is kept as it is; a list of ggplots (the parts of a combined plot) is recorded
# part by part.
.m61_snapshot <- function(p) {
  if (is.data.frame(p)) return(list(table = p))
  if (!inherits(p, "ggplot") && is.list(p)) return(lapply(p, .m61_snapshot))
  list(data = p$data,
       layers = lapply(seq_along(p$layers), function(i) ggplot2::layer_data(p, i)),
       labs = ggplot2::get_labs(p))
}

# The parts handed to cowplot::plot_grid() by a combined plot.
.m61_parts <- function(draw, package) {
  grab <- function(..., plotlist = NULL) {
    parts <- c(list(...), plotlist)
    parts[vapply(parts, function(x) inherits(x, "ggplot"), logical(1))]
  }
  testthat::with_mocked_bindings(draw(), plot_grid = grab, .package = package)
}

# The AC-009 reference calls: every plot of set P (each method) plus
# richClusterDot and richUpset, on inputs whose rows all pass the default
# selection and carry every drawn column. Each entry is a function that draws.
# `gsea = FALSE` leaves out the calls that need fgsea.
.m61_reference_calls <- function(dir, gsea = TRUE) {
  x3 <- .m61_fp3(); d3 <- x3@result; fc30 <- .m61_fc30()
  cmp <- compareResult(list(a = x3, b = .m61_fp3b()))
  seeded <- function(f) function() { set.seed(1); f() }
  calls <- list(
    richBar_rr = function() richBar(x3),
    richBar_df = function() richBar(d3),
    richDot_rr = function() richDot(x3),
    richDot_df = function() richDot(d3),
    richNetplot_rr = seeded(function() richNetplot(x3)),
    richNetplot_df = seeded(function() richNetplot(d3)),
    richNetwork_rr = seeded(function() richNetwork(x3)),
    richNetwork_df = seeded(function() richNetwork(d3, gene = x3@gene)),
    richNetmap = seeded(function() richNetmap(list(a = x3, b = .m61_r2b()))),
    richLollipop_rr = function() richLollipop(x3),
    richLollipop_df = function() richLollipop(d3),
    richCircle_rr = function() richCircle(x3),
    richCircle_df = function() richCircle(d3),
    richScatter_rr = function() richScatter(x3),
    richScatter_df = function() richScatter(d3),
    richTermSim_rr = function() richTermSim(x3),
    richTermSim_df = function() richTermSim(d3),
    richGeneDot_rr_fc = function() richGeneDot(x3, fc = fc30),
    richGeneDot_df_fc = function() richGeneDot(d3, fc = fc30),
    richGeneDot_rr_nofc = function() richGeneDot(x3),
    richGeneHeat_rr_fc = function() richGeneHeat(x3, fc = fc30),
    richGeneHeat_df_fc = function() richGeneHeat(d3, fc = fc30),
    richGeneHeat_rr_nofc = function() richGeneHeat(x3),
    richGeneBar_rr = function() richGeneBar(x3),
    richGeneBar_df = function() richGeneBar(d3),
    richVolcano_rr = function() richVolcano(x3),
    richVolcano_df = function() richVolcano(d3),
    richCompareDot = function() richCompareDot(cmp),
    richHeatmap_two = function() richHeatmap(list(x3, .m61_fp3b())),
    richHeatmap_ac008_data = function() richHeatmap(list(.m61_fp6(), .m61_r2()), returnData = TRUE),
    richHeatmap_ac008_p = function() richHeatmap(list(.m61_fp6(), .m61_r2()), usePadj = FALSE),
    richHeatmap_ac008_q = function() richHeatmap(list(.m61_fp6(), .m61_r2()), usePadj = TRUE),
    richClusterDot_ora = function() .m61_parts(function() richClusterDot(cmp), "cowplot"),
    richUpset = function() .m61_parts(function() richUpset(list(a = x3, b = .m61_r2b()),
                                                          set.size.show = FALSE), "richR")
  )
  if (isTRUE(gsea)) {
    g <- .m61_gsea_sig(dir); gd <- g@result; rk <- .m61_ranks(dir)
    calls <- c(calls, list(
      gsea_richNetplot = seeded(function() richNetplot(g)),
      gsea_richNetwork = seeded(function() richNetwork(g)),
      gsea_richLollipop = function() richLollipop(g),
      gsea_richCircle = function() richCircle(g),
      gsea_richNES = function() richNES(g),
      gsea_richNES_df = function() richNES(gd),
      gsea_richScatter = function() richScatter(g),
      gsea_richECDF = function() richECDF(g),
      gsea_richECDF_df = function() richECDF(gd),
      gsea_richTermSim = function() richTermSim(g),
      gsea_richGeneDot_fc = function() richGeneDot(g, fc = rk),
      gsea_richGeneDot_nofc = function() richGeneDot(g),
      gsea_richGeneHeat_fc = function() richGeneHeat(g, fc = rk),
      gsea_richGeneHeat_nofc = function() richGeneHeat(g),
      gsea_richGeneBar = function() richGeneBar(g),
      gsea_richVolcano = function() richVolcano(g),
      gsea_richClusterDot = function() .m61_parts(function() richClusterDot(g), "cowplot")
    ))
  }
  calls
}

# Run each reference call and record what it draws.
.m61_record <- function(calls) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  lapply(calls, function(f) .m61_snapshot(suppressMessages(f())))
}
