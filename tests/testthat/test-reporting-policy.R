# Reporting-cutoff policy (SPEC-RICHR-002 D1, REQ-015 ... REQ-022).
#
# One significance predicate decides every reported row: a row is kept when
# its raw p-value is non-NA and <= pvalue and, when padj is non-NULL, its
# adjusted p-value is non-NA and <= padj. Expected selections are written out
# by hand here, never copied from richR output.

# Synthetic rows: equality at 0.05, conflicting raw/adjusted values, NA raw p,
# NA adjusted p and p = 1.
policy_rows <- function() {
  data.frame(
    Annot = c("EQ", "RAW", "ADJ", "NARAW", "NAADJ", "ONE"),
    Term = c("equal", "raw only", "adjusted only", "NA raw", "NA adjusted", "one"),
    Annotated = rep(10, 6),
    Significant = rep(3, 6),
    RichFactor = rep(0.3, 6),
    Pvalue = c(0.05, 0.04, 0.06, NA, 0.01, 1),
    Padj = c(0.05, 0.06, 0.04, 0.01, NA, 1),
    GeneID = rep("g1,g2,g3", 6),
    stringsAsFactors = FALSE
  )
}

policy_keep <- function(pvalue, padj) {
  richR:::.sig_keep(policy_rows(), pvalue = pvalue, padj = padj)
}

test_that("the predicate keeps equality at both cutoffs and drops NA rows", {
  keep <- policy_keep(0.05, 0.05)
  expect_identical(keep, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_false(anyNA(keep))
})

test_that("the predicate with padj = NULL filters on the raw p-value only", {
  keep <- policy_keep(0.05, NULL)
  # NA adjusted p is dropped only while the adjusted cutoff is active.
  expect_identical(keep, c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_false(anyNA(keep))
})

test_that("the predicate applies both cutoffs when they conflict", {
  expect_identical(policy_keep(0.04, 0.06), c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(policy_keep(0.06, 0.04), c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
})

test_that("pvalue = 1 and padj = 1 keep every row with non-NA p-values, including p = 1", {
  expect_identical(policy_keep(1, 1), c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))
})

test_that("the predicate takes explicit column names and handles zero rows", {
  gsea_rows <- data.frame(pathway = c("a", "b"), pval = c(0.05, NA),
                          padj = c(0.05, 0.01))
  expect_identical(richR:::.sig_keep(gsea_rows, 0.05, 0.05, p_col = "pval",
                                     padj_col = "padj"), c(TRUE, FALSE))
  empty <- policy_rows()[0, ]
  expect_identical(richR:::.sig_keep(empty, 0.05, 0.05), logical(0))
  expect_identical(nrow(richR:::.filter_ora_result(
    empty, 0.05, 0.05, 1, 500, 1, 500, TRUE)), 0L)
})

test_that("the ORA filter uses the predicate and never creates NA rows", {
  rows <- policy_rows()
  out <- richR:::.filter_ora_result(rows, pvalue = 0.05, padj = 0.05, minSize = 1,
                                    maxSize = 500, minGSSize = 1, maxGSSize = 500,
                                    keepRich = TRUE)
  expect_identical(out$Annot, "EQ")
  unrestricted <- richR:::.filter_ora_result(rows, pvalue = 1, padj = 1, minSize = 1,
                                             maxSize = 500, minGSSize = 1,
                                             maxGSSize = 500, keepRich = TRUE)
  expect_false(anyNA(unrestricted$Annot))
  expect_setequal(unrestricted$Annot, c("EQ", "RAW", "ADJ", "ONE"))
  raw_only <- richR:::.filter_ora_result(rows, pvalue = 0.05, padj = NULL, minSize = 1,
                                         maxSize = 500, minGSSize = 1,
                                         maxGSSize = 500, keepRich = TRUE)
  expect_setequal(raw_only$Annot, c("EQ", "RAW", "NAADJ"))
})

test_that("the plot filter uses the predicate and never creates NA rows", {
  rows <- policy_rows()
  out <- richR:::.filter_and_top(rows, top = 10, pvalue = 0.05, padj = 0.05)
  expect_identical(out$Annot, "EQ")
  p <- richR::richDot(rows, pvalue = 1, padj = 1)
  expect_s3_class(p, "ggplot")
  expect_false(anyNA(p$data$Annot))
  expect_setequal(p$data$Annot, c("EQ", "RAW", "ADJ", "ONE"))
})

# GSEA result built directly: rows at the cutoffs and NA rows.
# getPathways() returns identifiers; filterPathways() returns the result rows.
policy_gsea <- function() {
  genes <- sprintf("g%02d", 1:12)
  fc <- stats::setNames(seq(2, -2, length.out = length(genes)), genes)
  res <- data.frame(
    pathway = c("EQ", "RAW", "ADJ", "NARAW", "NAADJ", "ONE"),
    pval = c(0.05, 0.04, 0.06, NA, 0.01, 1),
    padj = c(0.05, 0.06, 0.04, 0.01, NA, 1),
    log2err = rep(0.5, 6),
    ES = c(0.9, 0.8, -0.7, 0.6, 0.3, -0.5),
    NES = c(1.9, 1.7, -1.5, 1.3, 0.6, -1.1),
    size = rep(4L, 6),
    leadingEdge = rep("g01", 6),
    Term = c("equal", "raw only", "adjusted only", "NA raw", "NA adjusted", "one"),
    stringsAsFactors = FALSE
  )
  methods::new(
    "GSEAResult", result = res, pvalueCutoff = 1, pAdjustMethod = "BH",
    padjCutoff = 1, genenumber = length(fc), organism = character(),
    gene = names(fc), input = fc, ontology = character(), keytype = "",
    sep = ","
  )
}

test_that("getPathways() and filterPathways() keep exactly the predicate's rows", {
  x <- policy_gsea()
  # Defaults: padj_cutoff = 0.05, pval_cutoff = 0.05 -> both apply, equality kept.
  expect_identical(richR::getPathways(x), "EQ")
  expect_identical(richR::filterPathways(x)$pathway, "EQ")
  # padj_cutoff = NULL -> raw p-value only.
  expect_setequal(richR::getPathways(x, padj_cutoff = NULL),
                  c("EQ", "RAW", "NAADJ"))
  expect_setequal(richR::filterPathways(x, padj_cutoff = NULL)$pathway,
                  c("EQ", "RAW", "NAADJ"))
  # pval_cutoff < padj_cutoff -> both apply.
  expect_identical(richR::getPathways(x, padj_cutoff = 0.06, pval_cutoff = 0.04),
                   "RAW")
  expect_identical(richR::filterPathways(x, padj_cutoff = 0.06, pval_cutoff = 0.04)$pathway,
                   "RAW")
  expect_false(anyNA(richR::getPathways(x, padj_cutoff = 1, pval_cutoff = 1)))
  expect_setequal(richR::getPathways(x, padj_cutoff = 1, pval_cutoff = 1),
                  c("EQ", "RAW", "ADJ", "ONE"))
})

test_that("richGSEAcurve() selects pathways with the predicate", {
  skip_if_not_installed("fgsea")
  x <- policy_gsea()
  annot <- data.frame(GeneID = rep(sprintf("g%02d", 1:12), 1),
                      Term = rep(c("EQ", "RAW", "ADJ"), each = 4),
                      stringsAsFactors = FALSE)
  curve <- richR::richGSEAcurve(annot, x, pvalue = 0.05, padj = 0.05, return_data = TRUE)
  expect_identical(unname(curve$pathways), "EQ")
})

policy_fgsea <- function(...) {
  data.frame(
    pathway = c("EQ", "HIGH", "NA", "LOW"),
    pval = c(0.05, 0.06, NA, 0.01),
    NES = c(1, 1, 0, -1),
    size = c(2L, 2L, 2L, 2L),
    leadingEdge = I(list("A", "B", "C", "D"))
  )
}

policy_gsea_run <- function(...) {
  stats <- c(A = 2, B = 1, C = -1, D = -2)
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = c("EQ", "HIGH", "NA", "LOW"))
  testthat::with_mocked_bindings(
    fgseaMultilevel = policy_fgsea,
    richR::richGSEA(stats, annot, minSize = 1, table = FALSE, padj.method = "none", ...),
    .package = "richR"
  )
}

test_that("richGSEA() adjusts first and keeps equality at the cutoffs", {
  res <- policy_gsea_run(pvalue = 0.05, padj = 0.05)
  # padj.method = "none" leaves padj equal to pval, so the boundary row stays.
  expect_setequal(res@result$pathway, c("EQ", "LOW"))
  expect_equal(res@result$padj, res@result$pval)
})

test_that("pvalue = 1 and padj = 1 report every row with non-NA p-values", {
  gsea <- policy_gsea_run(pvalue = 1, padj = 1)
  expect_setequal(gsea@result$pathway, c("EQ", "HIGH", "LOW"))
  expect_false(anyNA(gsea@result$pathway))
  frame <- policy_rows()
  ora <- richR:::.filter_ora_result(frame, pvalue = 1, padj = 1, minSize = 1,
                                    maxSize = 500, minGSSize = 1, maxGSSize = 500,
                                    keepRich = TRUE)
  expected <- frame$Annot[!is.na(frame$Pvalue) & !is.na(frame$Padj)]
  expect_setequal(ora$Annot, expected)
  cmp <- richR::compareResult(list(A = frame), pvalue = 1, padj = 1)
  expect_setequal(cmp$Annot, expected)
})

test_that("stored cutoffs follow the supplied values", {
  annot <- data.frame(GeneID = paste0("G", 1:12),
                      Term = rep(c("T1", "T2", "T3"), each = 4))
  ora_null <- suppressMessages(richR::enrich(paste0("G", 1:4), annot, pvalue = 1,
                                             padj = NULL, minGSSize = 1))
  expect_identical(ora_null@padjCutoff, numeric(0))
  expect_null(richR::enrichmentInfo(ora_null)$options$padj)
  ora_set <- suppressMessages(richR::enrich(paste0("G", 1:4), annot, pvalue = 1,
                                            padj = 0.1, minGSSize = 1))
  expect_identical(ora_set@padjCutoff, 0.1)
  expect_identical(richR::enrichmentInfo(ora_set)$options$padj, 0.1)
  gsea_null <- policy_gsea_run(pvalue = 1, padj = NULL)
  expect_identical(gsea_null@padjCutoff, numeric(0))
  expect_null(richR::enrichmentInfo(gsea_null)$options$padj)
  gsea_set <- policy_gsea_run(pvalue = 1, padj = 0.1)
  expect_identical(gsea_set@padjCutoff, 0.1)
  expect_identical(richR::enrichmentInfo(gsea_set)$options$padj, 0.1)
})

# Comparison rows: a small raw p-value with an adjusted p-value of 1, and the
# reverse.
policy_compare_rows <- function() {
  data.frame(
    Annot = c("RAWSMALL", "ADJSMALL"), Term = c("raw small", "adjusted small"),
    Annotated = c(10, 10), Significant = c(3, 3), Pvalue = c(0.01, 1),
    Padj = c(1, 0.01), GeneID = c("g1", "g2"), stringsAsFactors = FALSE
  )
}

test_that("compareResult() has no Padj < 1 or Pvalue < 1 guard", {
  x <- list(A = policy_compare_rows())
  raw_only <- richR::compareResult(x, padj = NULL)
  expect_identical(raw_only$Annot, "RAWSMALL")
  adjusted <- richR::compareResult(x, pvalue = 1, padj = 0.05)
  expect_identical(adjusted$Annot, "ADJSMALL")
  boundary <- richR::compareResult(list(A = policy_rows()), pvalue = 0.05, padj = 0.05)
  expect_identical(boundary$Annot, "EQ")
})

policy_compare_gsea <- function() {
  genes <- sprintf("g%02d", 1:8)
  fc <- stats::setNames(seq(2, -2, length.out = length(genes)), genes)
  res <- data.frame(
    pathway = c("RAWSMALL", "ADJSMALL"), pval = c(0.01, 1), padj = c(1, 0.01),
    log2err = c(0.5, 0.5), ES = c(0.8, -0.6), NES = c(1.7, -1.2), size = c(4L, 4L),
    leadingEdge = c("g01", "g08"), stringsAsFactors = FALSE
  )
  methods::new(
    "GSEAResult", result = res, pvalueCutoff = 1, pAdjustMethod = "BH",
    padjCutoff = 1, genenumber = length(fc), organism = character(),
    gene = names(fc), input = fc, ontology = character(), keytype = "",
    sep = ","
  )
}

test_that("compareGSEA() has no Padj < 1 or Pvalue < 1 guard", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(GeneID = sprintf("g%02d", 1:8),
                      Term = rep(c("RAWSMALL", "ADJSMALL"), each = 4),
                      stringsAsFactors = FALSE)
  x <- list(A = policy_compare_gsea())
  raw_only <- richR::compareGSEA(x, annot, padj = NULL)
  expect_identical(unique(as.character(raw_only$data$Group)), "RAWSMALL")
  adjusted <- richR::compareGSEA(x, annot, pvalue = 1, padj = 0.05)
  expect_identical(unique(as.character(adjusted$data$Group)), "ADJSMALL")
})

test_that("richCompareDot(include.all = TRUE) skips the predicate", {
  x <- data.frame(
    Annot = c("a", "b", "c", "d"), Term = c("ta", "tb", "tc", "td"),
    Annotated = c(10, 10, 10, 10), Significant = c(2, 3, 4, 5),
    Pvalue = c(0.01, 1, NA, 0.2), Padj = c(0.02, 1, 0.5, NA),
    GeneID = c("g1", "g2", "g3", "g4"),
    group = factor(c("A", "A", "B", "B"), levels = c("A", "B", "C")),
    stringsAsFactors = FALSE
  )
  p <- richR::richCompareDot(x, include.all = TRUE)
  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$data), nrow(x))
  expect_setequal(as.character(p$data$Term), x$Term)
  # Without include.all the predicate applies: only row a passes 0.05.
  q <- richR::richCompareDot(x, pvalue = 0.05, padj = 0.05)
  expect_identical(as.character(q$data$Term), "ta")
})

# Shared validation rule: finite numeric scalar in (0, 1].
policy_bad_cutoffs <- list(0, -0.1, 1.1, NA, c(0.05, 0.1), "0.05")
policy_cutoff_error <- "must be between 0 and 1: one finite numeric scalar in \\(0, 1\\]"

test_that("ORA, GSEA, compare and plot entry points share one cutoff rule", {
  annot <- data.frame(GeneID = paste0("G", 1:12),
                      Term = rep(c("T1", "T2", "T3"), each = 4))
  genes <- paste0("G", 1:4)
  frame <- policy_rows()
  gsea <- policy_gsea()
  for (bad in policy_bad_cutoffs) {
    expect_error(richR::enrich(genes, annot, pvalue = bad, minGSSize = 1),
                 policy_cutoff_error)
    expect_error(richR::enrich(genes, annot, padj = bad, minGSSize = 1),
                 policy_cutoff_error)
    expect_error(policy_gsea_run(pvalue = bad), policy_cutoff_error)
    expect_error(policy_gsea_run(padj = bad), policy_cutoff_error)
    expect_error(richR::compareResult(list(A = frame), pvalue = bad), policy_cutoff_error)
    expect_error(richR::compareResult(list(A = frame), padj = bad), policy_cutoff_error)
    expect_error(richR::richDot(frame, pvalue = bad), policy_cutoff_error)
    expect_error(richR::richDot(frame, padj = bad), policy_cutoff_error)
    expect_error(richR::getPathways(gsea, pval_cutoff = bad), policy_cutoff_error)
    expect_error(richR::getPathways(gsea, padj_cutoff = bad), policy_cutoff_error)
  }
  expect_s4_class(suppressMessages(richR::enrich(genes, annot, pvalue = 1, padj = 1,
                                                 minGSSize = 1)), "richResult")
  expect_s4_class(policy_gsea_run(pvalue = 1, padj = 1), "GSEAResult")
  expect_s3_class(richR::compareResult(list(A = frame), pvalue = 1, padj = 1), "data.frame")
  expect_s3_class(richR::richDot(frame, pvalue = 1, padj = 1), "ggplot")
  expect_type(richR::getPathways(gsea, pval_cutoff = 1, padj_cutoff = 1), "character")
})

# Defaults (REQ-014): every exported function, S4 generic and S4 method that
# takes pvalue and padj defaults to pvalue = 0.05 and padj = 0.05.
policy_cutoff_defaults <- function(fun) {
  fm <- formals(fun)
  if (!"padj" %in% names(fm)) return(NULL)
  list(pvalue = fm$pvalue, padj = fm$padj)
}

test_that("exported functions, generics and methods default to pvalue = padj = 0.05", {
  ns <- asNamespace("richR")
  checked <- 0L
  for (nm in getNamespaceExports("richR")) {
    fun <- get(nm, envir = ns)
    if (!is.function(fun)) next
    funs <- list(fun)
    if (methods::isGeneric(nm, where = ns)) {
      funs <- c(funs, lapply(methods::findMethods(nm, where = ns),
                             methods::unRematchDefinition))
    }
    for (f in funs) {
      d <- policy_cutoff_defaults(f)
      if (is.null(d)) next
      checked <- checked + 1L
      expect_identical(d$padj, 0.05, info = nm)
      expect_identical(d$pvalue, 0.05, info = nm)
    }
  }
  expect_gt(checked, 50L)
})

test_that("ORA under the default cutoffs applies the adjusted cutoff; padj = NULL is raw-only", {
  query <- .char_read("ora_query.csv")$GeneID
  annot <- .char_read("ora_annot.csv", character = TRUE)
  # The fixture's smallest Padj is 0.0690, so padj = 0.05 reports nothing.
  default <- suppressMessages(richR::enrich(query, annot, pvalue = 0.5, minSize = 2,
                                            minGSSize = 3, keepRich = TRUE))
  expect_identical(nrow(default@result), 0L)
  expect_identical(default@padjCutoff, 0.05)
  raw_only <- suppressMessages(richR::enrich(query, annot, pvalue = 0.5, padj = NULL,
                                             minSize = 2, minGSSize = 3, keepRich = TRUE))
  expect_gt(nrow(raw_only@result), 0L)
  expect_identical(raw_only@padjCutoff, numeric(0))
})

# ---- Post-review: empty selections and cutoff validation ---------------------

policy_curve_annot <- function() {
  data.frame(GeneID = sprintf("g%02d", 1:12),
             Term = rep(c("EQ", "RAW", "ADJ"), each = 4),
             stringsAsFactors = FALSE)
}

test_that("richGSEAcurve() stops clearly when no pathway passes the cutoffs", {
  skip_if_not_installed("fgsea")
  expect_error(richR::richGSEAcurve(policy_curve_annot(), policy_gsea(),
                                    pvalue = 0.01, padj = 0.01),
               "no pathway passes the cutoffs")
})

test_that("compareGSEA() stops clearly when no pathway passes the cutoffs", {
  skip_if_not_installed("fgsea")
  annot <- data.frame(GeneID = sprintf("g%02d", 1:8),
                      Term = rep(c("RAWSMALL", "ADJSMALL"), each = 4),
                      stringsAsFactors = FALSE)
  x <- list(A = policy_compare_gsea())
  expect_error(richR::compareGSEA(x, annot, pvalue = 1e-12, padj = 1e-12),
               "no pathway passes the cutoffs")
  # the selected pathway is removed by the cutoffs in every group
  expect_error(suppressWarnings(richR::compareGSEA(x, annot, pathway = "ADJSMALL",
                                                   padj = NULL)),
               "no pathway passes the cutoffs")
})

test_that("richGSEAcurve() validates cutoffs for explicit selections too", {
  x <- policy_gsea()
  for (bad in policy_bad_cutoffs) {
    expect_error(richR::richGSEAcurve(policy_curve_annot(), x, pathways = "EQ",
                                      pvalue = bad), policy_cutoff_error)
    expect_error(richR::richGSEAcurve(policy_curve_annot(), x, pathway_pattern = "EQ",
                                      padj = bad), policy_cutoff_error)
  }
})

test_that("richVolcano() and its aliases validate cutoffs", {
  frame <- policy_rows()
  gsea <- policy_gsea()
  for (bad in policy_bad_cutoffs) {
    expect_error(richR::richVolcano(frame, pvalue = bad), policy_cutoff_error)
    expect_error(richR::richVolcano(frame, padj = bad), policy_cutoff_error)
    expect_error(richR::ggvolcano(gsea, pvalue = bad), policy_cutoff_error)
    expect_error(richR::ggVolcano(gsea, padj = bad), policy_cutoff_error)
  }
  expect_s3_class(richR::richVolcano(frame, pvalue = 1, padj = 1), "ggplot")
  expect_s3_class(richR::richVolcano(gsea, padj = NULL), "ggplot")
})

# ---- Post-review: show() and summary() report the stored cutoffs -------------

test_that("show() states the inclusive rule with the stored cutoffs", {
  ora <- .char_richResult()
  expect_output(show(ora), "Significant:  3 terms (pvalue <= 0.05)", fixed = TRUE)
  ora@padjCutoff <- 0.1
  expect_output(show(ora), "(pvalue <= 0.05, padj <= 0.1)", fixed = TRUE)
  gsea <- policy_gsea()
  expect_output(show(gsea), "Significant:  6 pathways (pvalue <= 1, padj <= 1)",
                fixed = TRUE)
})

test_that("summary(<GSEAResult>) counts the rows passing the stored cutoffs", {
  # stored pvalue = padj = 1: EQ, RAW, ADJ and ONE pass; the NA rows do not
  expect_output(summary(policy_gsea()), "Total significant biological term is: 4",
                fixed = TRUE)
  # stored pvalue = 0.05 and no adjusted cutoff: p1 and p2 pass
  expect_output(summary(.char_gsea_result()), "Total significant biological term is: 2",
                fixed = TRUE)
})
