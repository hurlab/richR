# GSEA gene-set identity (SPEC-RICHR-002 D3; REQ-003 .. REQ-006).
#
# Gene sets are keyed by identifier and labelled afterwards: `@result$pathway`
# holds the identifier and `@result$Term` the label. Distinct identifiers that
# share a label stay distinct gene sets.

identity_scores <- function() {
  genes <- sprintf("g%02d", 1:40)
  stats::setNames(seq(3, -3, length.out = length(genes)), genes)
}

identity_gsea <- function(annot, ...) {
  out <- NULL
  utils::capture.output(
    out <- richGSEA(identity_scores(), annot, minSize = 1, maxSize = 50,
                    pvalue = 1, padj = NULL, seed = 7, nproc = 1, ...)
  )
  out
}

genes_of <- function(idx) sprintf("g%02d", idx)

# ---- REQ-003: memberships and labels -----------------------------------------

test_that("two identifiers sharing a label stay two gene sets with their own sizes", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(
    GeneID = genes_of(c(1:5, 33:40)),
    Term = rep(c("S1", "S2"), c(5, 8)),
    Annot = "shared label",
    stringsAsFactors = FALSE
  )
  res <- identity_gsea(annot)@result
  expect_setequal(res$pathway, c("S1", "S2"))
  expect_identical(res$size[match(c("S1", "S2"), res$pathway)], c(5L, 8L))
  expect_identical(res$Term, rep("shared label", 2))
  expect_identical(colnames(res)[1:9],
                   c("pathway", "pval", "padj", "log2err", "ES", "NES", "size",
                     "leadingEdge", "Term"))
})

test_that("an identifier without a label is its own Term", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(
    GeneID = genes_of(c(1:5, 33:40)),
    Term = rep(c("L1", "N1"), c(5, 8)),
    Annot = rep(c("labelled", NA), c(5, 8)),
    stringsAsFactors = FALSE
  )
  res <- identity_gsea(annot)@result
  expect_identical(res$Term[res$pathway == "N1"], "N1")
  expect_identical(res$Term[res$pathway == "L1"], "labelled")
})

test_that("conflicting labels keep the first label in row order and warn once", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(
    GeneID = genes_of(c(1:5, 33:40)),
    Term = rep(c("C1", "OK"), c(5, 8)),
    Annot = c("first label", "first label", "second label", "second label",
              "first label", rep("fine", 8)),
    stringsAsFactors = FALSE
  )
  expect_warning(res <- identity_gsea(annot)@result, "C1")
  expect_identical(res$Term[res$pathway == "C1"], "first label")
  expect_identical(res$Term[res$pathway == "OK"], "fine")
})

test_that("the identifier column is found by name when the label column comes first", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(
    GeneID = genes_of(c(1:5, 33:40)),
    Annot = rep(c("alpha", "beta"), c(5, 8)),
    Term = rep(c("ID_A", "ID_B"), c(5, 8)),
    stringsAsFactors = FALSE
  )
  res <- identity_gsea(annot)@result
  expect_setequal(res$pathway, c("ID_A", "ID_B"))
  expect_identical(res$Term[match(c("ID_A", "ID_B"), res$pathway)], c("alpha", "beta"))
})

test_that("an unnamed custom frame uses column 2 as the identifier", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(
    gene = genes_of(c(1:5, 33:40)),
    set = rep(c("X1", "X2"), c(5, 8)),
    stringsAsFactors = FALSE
  )
  res <- identity_gsea(annot)@result
  expect_setequal(res$pathway, c("X1", "X2"))
  expect_identical(res$Term, res$pathway)
})

test_that("the bundled GSEA fixture yields T-identifiers with names in Term", {
  skip_if_not_installed("fgsea")
  annot <- utils::read.csv(test_path("fixtures", "gsea_annot.csv"),
                           colClasses = "character")
  ranked <- utils::read.csv(test_path("fixtures", "gsea_ranked_a.csv"))
  out <- NULL
  utils::capture.output(
    out <- richGSEA(stats::setNames(ranked$score, ranked$GeneID), annot,
                    minSize = 3, maxSize = 20, pvalue = 1, seed = 42, nproc = 1)
  )
  res <- out@result
  expect_true(all(res$pathway %in% paste0("T", 1:5)))
  labels <- unique(annot[, c("Term", "Annot")])
  expect_identical(res$Term, labels$Annot[match(res$pathway, labels$Term)])
})

# ---- REQ-004: KEGG join ------------------------------------------------------

kegg_path <- function() {
  env <- new.env(parent = emptyenv())
  utils::data("path", package = "richR", envir = env)
  env$path
}

expect_kegg_levels <- function(res, keys) {
  path <- kegg_path()
  ref <- path[match(keys, path$ko), , drop = FALSE]
  expect_identical(as.character(res$Level3), as.character(ref$Level3))
  expect_identical(as.character(res$Level2), as.character(ref$Level2))
  expect_identical(as.character(res$Level1), as.character(ref$Level1))
}

test_that("KEGG levels join 5-digit identifiers on ko", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(GeneID = genes_of(c(1:5, 33:40)),
                      PATH = rep(c("00010", "00020"), c(5, 8)),
                      stringsAsFactors = FALSE)
  res <- identity_gsea(annot, KEGG = TRUE)@result
  expect_kegg_levels(res, res$pathway)
})

test_that("KEGG levels join organism-prefixed identifiers on their 5 digits", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(GeneID = genes_of(c(1:5, 33:40)),
                      PATH = rep(c("hsa00010", "hsa00020"), c(5, 8)),
                      stringsAsFactors = FALSE)
  res <- identity_gsea(annot, KEGG = TRUE)@result
  expect_setequal(res$pathway, c("hsa00010", "hsa00020"))
  expect_kegg_levels(res, sub("^hsa", "", res$pathway))
})

test_that("KEGG levels join a pathway name on Level3", {
  skip_if_not_installed("fgsea")
  path <- kegg_path()
  annot <- data.frame(GeneID = genes_of(c(1:5, 33:40)),
                      Term = rep(path$Level3[1:2], c(5, 8)),
                      stringsAsFactors = FALSE)
  res <- identity_gsea(annot, KEGG = TRUE)@result
  expect_kegg_levels(res, path$ko[match(res$pathway, path$Level3)])
})

# ---- hand-built D3 result used by the consumer and selector tests -------------
# Identifiers 04110 (a KEGG-style id), A1 and A2 (sharing the label "shared")
# and U1 (a unique label).

selector_annot <- function() {
  data.frame(
    GeneID = genes_of(c(1:5, 6:10, 11:15, 33:40)),
    Term = rep(c("04110", "A1", "A2", "U1"), c(5, 5, 5, 8)),
    Annot = rep(c("Cell cycle", "shared", "shared", "unique label"), c(5, 5, 5, 8)),
    stringsAsFactors = FALSE
  )
}

selector_result <- function(pval = c(0.001, 0.002, 0.003, 0.004)) {
  fc <- identity_scores()
  res <- data.frame(
    pathway = c("04110", "A1", "A2", "U1"), pval = pval, padj = pval * 2,
    log2err = rep(0.5, 4), ES = c(0.9, 0.8, 0.7, -0.9),
    NES = c(1.9, 1.8, 1.7, -1.9), size = c(5L, 5L, 5L, 8L),
    leadingEdge = c("g01", "g06", "g11", "g40"),
    Term = c("Cell cycle", "shared", "shared", "unique label"),
    stringsAsFactors = FALSE
  )
  methods::new(
    "GSEAResult", result = res, pvalueCutoff = 1, pAdjustMethod = "BH",
    padjCutoff = numeric(), genenumber = length(fc), organism = character(),
    gene = names(fc), input = fc, ontology = character(), keytype = "",
    sep = ","
  )
}

# ---- REQ-005: consumers ------------------------------------------------------

consumer_annot <- function() {
  data.frame(
    GeneID = genes_of(c(1:5, 6:10, 33:40)),
    Term = rep(c("04110", "S2", "S3"), c(5, 5, 8)),
    Annot = rep(c("Cell cycle", "Second set", "Third set"), c(5, 5, 8)),
    stringsAsFactors = FALSE
  )
}

consumer_result <- function(table = TRUE) {
  identity_gsea(consumer_annot(), table = table)
}

legacy_result <- function() {
  res <- consumer_result()
  res@result$Term <- NULL
  res
}

test_that("getGenes() returns the same genes for a pathway identifier and its label", {
  skip_if_not_installed("fgsea")
  for (table in c(TRUE, FALSE)) {
    res <- consumer_result(table)
    by_id <- getGenes(res, term = "S3")
    by_label <- getGenes(res, term = "Third set")
    expect_true(length(by_id) > 0)
    expect_identical(by_label, by_id)
  }
})

test_that("show() prints the Term label, or pathway for a legacy result", {
  skip_if_not_installed("fgsea")
  expect_output(show(consumer_result()), "Third set")
  out <- utils::capture.output(show(legacy_result()))
  expect_true(any(grepl("S3", out)))
  expect_false(any(grepl("Third set", out)))
})

test_that("ORA-shaped frames keep the identifier as Annot and the label as Term", {
  skip_if_not_installed("fgsea")
  res <- consumer_result()
  frame <- richR:::.gsea_to_ora_frame(res)
  expect_identical(frame$Annot, res@result$pathway)
  expect_identical(frame$Term, res@result$Term)
  legacy <- richR:::.gsea_to_ora_frame(legacy_result())
  expect_identical(legacy$Annot, legacy_result()@result$pathway)
  expect_identical(legacy$Term, legacy_result()@result$pathway)
})

test_that(".standardise_gsea_cols() keeps pathway as identifier and Term as label", {
  skip_if_not_installed("fgsea")
  res <- consumer_result()
  std <- richR:::.standardise_gsea_cols(res@result)
  expect_identical(std$pathway, res@result$pathway)
  expect_identical(std$Term, res@result$Term)
})

test_that("richClusterDot(method = 'auto') takes the GSEA branch for a D3 result", {
  skip_if_not_installed("fgsea")
  route <- function(x) {
    testthat::with_mocked_bindings(
      richClusterDot(x, method = "auto"),
      .ggcluster_gsea = function(data, ...) "gsea",
      .ggcluster_enrich = function(data, ...) "enrich",
      .package = "richR"
    )
  }
  expect_identical(route(consumer_result()), "gsea")
  # a hand-built GSEAResult carrying both pathway and Term
  expect_identical(route(selector_result()), "gsea")
  expect_identical(route(legacy_result()), "gsea")
})

test_that("the GSEA cluster plot labels pathways by Term", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("cowplot")
  # capture the main panel handed to cowplot::plot_grid()
  pathway_labels <- function(x) {
    parts <- testthat::with_mocked_bindings(
      richClusterDot(x, method = "gsea"),
      plot_grid = function(...) list(...),
      .package = "cowplot"
    )
    p <- parts[[1]]
    expect_s3_class(p, "ggplot")
    unlist(lapply(p$layers, function(l) {
      if (is.data.frame(l$data) && "label" %in% names(l$data)) as.character(l$data$label)
    }))
  }
  labels <- pathway_labels(consumer_result())
  expect_true(all(c("Cell cycle", "Second set", "Third set") %in% labels))
  expect_false(any(c("04110", "S2", "S3") %in% labels))
  # a hand-built result whose identifiers differ from its labels
  labels <- pathway_labels(selector_result())
  expect_true(all(c("Cell cycle", "shared", "unique label") %in% labels))
  expect_false(any(c("04110", "A1", "A2", "U1") %in% labels))
  # a legacy result without Term shows the pathway
  expect_true("S3" %in% pathway_labels(legacy_result()))
})

test_that("network adapters build their frames by name", {
  skip_if_not_installed("fgsea")
  for (table in c(TRUE, FALSE)) {
    res <- consumer_result(table)
    netplot_frame <- testthat::with_mocked_bindings(
      richNetplot(res),
      ggrich_internal = function(object, ...) object,
      .package = "richR"
    )
    network_frame <- testthat::with_mocked_bindings(
      richNetwork(res),
      ggnetwork_internal = function(object, ...) object,
      .package = "richR"
    )
    for (frame in list(netplot_frame, network_frame)) {
      expect_identical(frame$Annot, res@result$pathway)
      expect_identical(frame$Term, res@result$Term)
      expect_identical(frame$Pvalue, res@result$pval)
      expect_identical(frame$Padj, res@result$padj)
      expect_true(all(nzchar(frame$GeneID)))
      expect_true(is.character(frame$GeneID))
    }
    expect_identical(network_frame$Significant, res@result$size)
  }
})

test_that("getPathwayStats() and searchPathways() expose and match the Term label", {
  skip_if_not_installed("fgsea")
  res <- consumer_result()
  stats_df <- getPathwayStats(res)
  expect_identical(stats_df$Term, res@result$Term)
  expect_identical(stats_df$pathway, res@result$pathway)
  hits <- searchPathways(res, "third")
  expect_identical(unname(hits), "S3")
  expect_identical(names(hits), "Third set")
  expect_identical(unname(searchPathways(res, "^S2$")), "S2")
})

# ---- REQ-006: selectors and curve facets --------------------------------------

curve_pathways <- function(selector) {
  out <- richGSEAcurve(selector_annot(), selector_result(), pathways = selector,
                       return_data = TRUE)
  sort(unique(as.character(out$data$pathway)))
}

test_that("richGSEAcurve() resolves identifiers, prefixed KEGG identifiers and unique labels", {
  skip_if_not_installed("fgsea")
  expect_identical(curve_pathways("04110"), "04110")
  expect_identical(curve_pathways("hsa04110"), "04110")
  expect_identical(curve_pathways("unique label"), "U1")
  expect_identical(curve_pathways("Cell cycle"), "04110")
})

test_that("richGSEAcurve() errors on a shared label and on an unknown selector", {
  skip_if_not_installed("fgsea")
  expect_error(curve_pathways("shared"), "A1, A2")
  expect_error(curve_pathways("no such pathway"), "no such pathway")
})

test_that("richGSEAcurve() facets show Term while the data keep pathway", {
  skip_if_not_installed("fgsea")
  out <- richGSEAcurve(selector_annot(), selector_result(),
                       pathways = c("04110", "U1"), return_data = TRUE)
  key <- unique(out$data[, c("Group", "pathway")])
  expect_identical(key$Group[match(c("04110", "U1"), key$pathway)],
                   c("Cell cycle", "unique label"))
})

test_that("richGSEAcurve() pattern matching looks at pathway and Term", {
  skip_if_not_installed("fgsea")
  out <- richGSEAcurve(selector_annot(), selector_result(),
                       pathway_pattern = "unique|04110", return_data = TRUE)
  expect_setequal(unique(as.character(out$data$pathway)), c("04110", "U1"))
})

test_that("getGenes() resolves selectors like the curve functions", {
  skip_if_not_installed("fgsea")
  res <- selector_result()
  expect_identical(getGenes(res, term = "hsa04110"), getGenes(res, term = "04110"))
  expect_identical(getGenes(res, term = "Cell cycle"), "g01")
  expect_error(getGenes(res, term = "shared"), "A1, A2")
  expect_error(getGenes(res, term = "no such pathway"), "no such pathway")
})

test_that("compareGSEA() resolves pathway selectors over all groups and warns on filtered ones", {
  skip_if_not_installed("fgsea")
  groups <- list(
    G1 = selector_result(),
    G2 = selector_result(pval = c(0.001, 0.002, 0.003, 0.9))
  )
  p <- compareGSEA(groups, selector_annot(), pathway = c("hsa04110", "unique label"),
                   pvalue = 0.05)
  key <- unique(p$data[, c("group", "Group", "pathway")])
  expect_setequal(paste(key$group, key$pathway), c("G1 04110", "G2 04110", "G1 U1"))
  expect_identical(unique(key$Group[key$pathway == "04110"]), "Cell cycle")
  expect_error(compareGSEA(groups, selector_annot(), pathway = "shared"), "A1, A2")
  expect_error(compareGSEA(groups, selector_annot(), pathway = "nothing"), "nothing")
  filtered <- list(G1 = selector_result(pval = c(0.001, 0.002, 0.003, 0.9)),
                   G2 = selector_result(pval = c(0.001, 0.002, 0.003, 0.9)))
  expect_warning(
    compareGSEA(filtered, selector_annot(), pathway = c("04110", "U1"), pvalue = 0.05),
    "U1"
  )
})

test_that("compareGSEA() facets show Term while the data keep pathway", {
  skip_if_not_installed("fgsea")
  p <- compareGSEA(list(G1 = selector_result()), selector_annot(), pvalue = 0.05)
  key <- unique(p$data[, c("Group", "pathway")])
  expect_identical(key$Group[match(c("04110", "U1"), key$pathway)],
                   c("Cell cycle", "unique label"))
  # the two identifiers sharing the label "shared" stay two facets
  expect_identical(length(unique(key$Group[key$pathway %in% c("A1", "A2")])), 2L)
})

# ---- Post-review: distinct pathways sharing a label stay distinct in plots ----
# selector_result(): A1 and A2 are two pathways labelled "shared" with different
# p-values; every plot must draw them as two categories with their own values.

collision_labels <- c("Cell cycle", "shared [A1]", "shared [A2]", "unique label")

test_that("GSEA plot frames give pathways sharing a label distinct display labels", {
  res <- selector_result()
  frame <- richR:::.gsea_to_ora_frame(res)
  expect_identical(frame$Annot, res@result$pathway)
  expect_identical(frame$Term, collision_labels)
  std <- richR:::.standardise_gsea_cols(res@result)
  expect_identical(std$pathway, res@result$pathway)
  expect_identical(std$Term, collision_labels)
})

test_that("GSEA plot methods draw one category per pathway when labels collide", {
  res <- selector_result()
  plot_terms <- function(p) sort(unique(as.character(p$data$Term)))
  expect_identical(plot_terms(richLollipop(res)), sort(collision_labels))
  expect_identical(plot_terms(richScatter(res)), sort(collision_labels))
  expect_identical(plot_terms(richGeneBar(res)), sort(collision_labels))
  expect_identical(plot_terms(richCircle(res)), sort(collision_labels))
  expect_identical(plot_terms(richNES(res)), sort(collision_labels))
  expect_identical(plot_terms(richECDF(res)), sort(collision_labels))
  expect_identical(plot_terms(richVolcano(res)), sort(collision_labels))
})

test_that("gene plots key significance by pathway, not by the shared label", {
  res <- selector_result()
  expected <- stats::setNames(-log10(res@result$padj), collision_labels)
  per_term <- function(d) {
    d <- d[!is.na(d$neg_log_p), ]
    vapply(names(expected), function(t) {
      v <- unique(d$neg_log_p[as.character(d$Term) == t])
      if (length(v) == 1L) v else NA_real_
    }, numeric(1))
  }
  expect_equal(per_term(richGeneDot(res)$data), expected)
  expect_equal(per_term(richGeneHeat(res)$data), expected)
})

test_that("gene plots map each gene row to its own term row's p-value", {
  # two ORA rows sharing a Term label; each gene belongs to one row only
  df <- data.frame(Annot = c("X", "Y"), Term = c("dup", "dup"),
                   Annotated = c(5, 5), Significant = c(1, 1),
                   Pvalue = c(0.001, 0.002), Padj = c(0.01, 0.02),
                   GeneID = c("ga", "gb"), stringsAsFactors = FALSE)
  dot <- richGeneDot(df, pvalue = 1, padj = 1)$data
  expect_equal(dot$neg_log_p[match(c("ga", "gb"), dot$Gene)], -log10(c(0.01, 0.02)))
})

# ---- Post-review: richGSEAplot() resolves term like the other selectors -------

# A GSEAResult computed before gene sets were keyed by identifier: `pathway`
# holds the label and there is no Term column.
label_keyed_result <- function() {
  res <- consumer_result()
  annot <- consumer_annot()
  res@result$pathway <- annot$Annot[match(res@result$pathway, annot$Term)]
  res@result$Term <- NULL
  res
}

test_that("richGSEAplot() draws the same gene set for an identifier and its label", {
  skip_if_not_installed("fgsea")
  res <- consumer_result()
  by_id <- richGSEAplot("04110", consumer_annot(), res)
  expect_no_warning(by_label <- richGSEAplot("Cell cycle", consumer_annot(), res))
  expect_identical(by_label$data, by_id$data)
  by_kegg <- richGSEAplot("hsa04110", consumer_annot(), res)
  expect_identical(by_kegg$data, by_id$data)
})

test_that("richGSEAplot() supports a label-keyed result from before 0.2.0", {
  skip_if_not_installed("fgsea")
  by_id <- richGSEAplot("04110", consumer_annot(), consumer_result())
  expect_no_warning(legacy <- richGSEAplot("Cell cycle", consumer_annot(),
                                           label_keyed_result()))
  expect_identical(legacy$data, by_id$data)
})

test_that("richGSEAplot() rejects an ambiguous or unknown term", {
  skip_if_not_installed("fgsea")
  expect_error(richGSEAplot("shared", consumer_annot(), selector_result()),
               "matches several pathways: A1, A2")
  expect_error(richGSEAplot("nope", consumer_annot(), consumer_result()),
               "'nope' was not found")
})
