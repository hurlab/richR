# Characterization baseline (REQ-004): richLevel(), compareResult(),
# compareGSEA() (padj = NULL, and non-NULL padj after the SPEC-RICHR-002 B1 fix)
# and richCluster(), as produced by richR 0.1.5. Known defects are pinned as
# they stand (REQ-005): a later fix must show up as a snapshot change.

compare_query_a <- function() .char_read("ora_query.csv")$GeneID
compare_query_b <- function() sprintf("G%02d", c(3, 4, 5, 9, 10, 11, 12, 16, 17, 18, 19, 20))
compare_annot <- function() .char_read("ora_annot.csv", character = TRUE)

compare_ora <- function(genes) {
  suppressMessages(enrich(
    genes, compare_annot(), pvalue = 0.5, padj = NULL, minSize = 2, minGSSize = 3,
    keepRich = TRUE
  ))
}

test_that("richLevel() pins Level2 and Level3 enrichment on bundled KEGG levels", {
  kegg <- .char_read("kegg_annot.csv", character = TRUE)
  ko <- data.frame(GeneID = kegg$GeneID, PATH = kegg$Term, stringsAsFactors = FALSE)
  level2 <- suppressMessages(richLevel(
    compare_query_a(), ko, level = "Level2", pvalue = 0.5, padj = NULL, minSize = 2, minGSSize = 3
  ))
  level3 <- suppressMessages(richLevel(
    compare_query_a(), ko, level = "Level3", pvalue = 0.5, padj = NULL, minSize = 2, minGSSize = 3
  ))
  expect_s4_class(level2, "richResult")
  expect_s4_class(level3, "richResult")
  expect_snapshot_value(
    list(
      level2 = list(result = .char_pin_df(level2@result), detail = .char_pin_df(level2@detail)),
      level3 = list(result = .char_pin_df(level3@result), detail = .char_pin_df(level3@detail))
    ),
    style = "deparse", tolerance = 1e-8
  )
})

test_that("compareResult() pins group stacking, cutoffs, include.all and extra input types", {
  first <- compare_ora(compare_query_a())
  second <- compare_ora(compare_query_b())
  .char_nonempty(first@result)
  .char_nonempty(second@result)

  by_pvalue <- compareResult(list(first = first, second = second), pvalue = 0.5, padj = NULL)
  by_padj <- compareResult(list(first = first, second = second), pvalue = 1, padj = 0.2)
  include_all <- compareResult(
    list(first = first, second = second, none = NULL), pvalue = 0.5, padj = NULL, include.all = TRUE
  )
  from_frames <- compareResult(
    list(first = first@result, second = second@result), pvalue = 0.5, padj = NULL
  )
  expect_identical(levels(include_all$group), c("first", "second", "none"))
  expect_snapshot_value(
    list(
      by_pvalue = .char_pin_df(by_pvalue),
      by_padj = .char_pin_df(by_padj),
      include_all = .char_pin_df(include_all),
      from_frames = .char_pin_df(from_frames)
    ),
    style = "deparse", tolerance = 1e-8
  )
})

gsea_pair <- function() {
  skip_if_not_installed("fgsea")
  annot <- .char_read("gsea_annot.csv", character = TRUE)
  run <- function(name) {
    df <- .char_read(name)
    out <- NULL
    utils::capture.output(
      out <- richGSEA(stats::setNames(df$score, df$GeneID), annot, minSize = 3,
                      maxSize = 20, pvalue = 1, padj = 1, seed = 42, nproc = 1)
    )
    out
  }
  list(annot = annot, results = list(A = run("gsea_ranked_a.csv"), B = run("gsea_ranked_b.csv")))
}

test_that("compareGSEA() with padj = NULL keeps the shape of the enrichment-curve data", {
  fx <- gsea_pair()
  plot_obj <- compareGSEA(fx$results, fx$annot, padj = NULL)
  expect_s3_class(plot_obj, "ggplot")
  curve <- .char_nonempty(plot_obj$data)
  # SPEC-RICHR-002 D3: facets (Group) show the Term label; the identifier is
  # kept in a `pathway` column.
  expect_identical(
    vapply(curve, function(col) class(col)[1L], character(1)),
    c(x = "numeric", y = "numeric", Group = "character", group = "character",
      pathway = "character")
  )
  expect_setequal(curve$Group, c("top pathway", "down pathway", "up pathway"))
  expect_setequal(curve$pathway, c("T4", "T2", "T1"))
  expect_setequal(curve$group, c("A", "B"))
  expect_identical(length(plot_obj$layers), 6L)
})

test_that("compareGSEA() with padj = NULL pins the enrichment-curve data behind the plot", {
  .char_skip_unless_recorded_fgsea()
  fx <- gsea_pair()
  plot_obj <- compareGSEA(fx$results, fx$annot, padj = NULL)
  expect_s3_class(plot_obj, "ggplot")
  .char_nonempty(plot_obj$data)
  expect_snapshot_value(
    list(
      fgsea_version = as.character(utils::packageVersion("fgsea")),
      curve = .char_pin_df(plot_obj$data),
      layers = length(plot_obj$layers)
    ),
    style = "deparse", tolerance = 1e-8
  )
})

test_that("compareGSEA() with a non-NULL padj applies both cutoffs to the fixture pair", {
  # SPEC-RICHR-002 B1: a non-NULL `padj` used to fail with "object 'Padj' not
  # found"; it now applies the raw and the adjusted cutoff together. On this
  # fixture every pathway with pval < 0.05 also has padj < 0.05, so the curve
  # data equal the padj = NULL curve.
  fx <- gsea_pair()
  with_padj <- compareGSEA(fx$results, fx$annot, padj = 0.05)
  without_padj <- compareGSEA(fx$results, fx$annot, padj = NULL)
  expect_s3_class(with_padj, "ggplot")
  curve <- .char_nonempty(with_padj$data)
  expect_identical(curve, without_padj$data)
  expect_snapshot_value(
    list(drawn = sort(unique(paste(curve$group, curve$Group, sep = "@")))),
    style = "deparse"
  )
})

test_that("richCluster() pins kappa clusters for data.frame and richResult input", {
  terms <- .char_read("cluster_terms.csv")
  as_rich <- function(df) {
    genes <- unique(unlist(strsplit(df$GeneID, ",", fixed = TRUE)))
    methods::new(
      "richResult", result = df, detail = data.frame(), pvalueCutoff = 0.1,
      pAdjustMethod = "BH", padjCutoff = 0.1, genenumber = length(genes),
      organism = "", ontology = "", gene = genes, keytype = "", sep = ","
    )
  }
  by_annot <- richCluster(terms, gene = NULL, useTerm = FALSE, cutoff = 0.1,
                          overlap = 0, minSize = 2)
  by_term <- richCluster(terms, gene = NULL, useTerm = TRUE, cutoff = 0.1,
                         overlap = 0, minSize = 2)
  from_rich <- richCluster(as_rich(terms), gene = NULL, useTerm = FALSE, cutoff = 0.1,
                           overlap = 0, minSize = 2)
  expect_snapshot_value(
    list(
      by_annot = .char_pin_df(by_annot),
      by_term = .char_pin_df(by_term),
      from_rich = .char_pin_df(from_rich)
    ),
    style = "deparse", tolerance = 1e-8
  )
})
