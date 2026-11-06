# SPEC-RICHR-012: exact leading-edge identifiers and comparison group identity.
review_identity_groups <- function() {
  fc <- stats::setNames(seq(2.3, -2.3, length.out = 24), sprintf("g%02d", 1:24))
  result <- data.frame(pathway = c("P1", "P2"), pval = c(.01, .03),
                       padj = c(.02, .09), NES = c(1.9, 1.7),
                       leadingEdge = c("g01", "g05"))
  a <- methods::new("GSEAResult", result = result, gene = names(fc), input = fc, sep = ",")
  b <- a
  b@result$NES <- c(-2.4, -3.1)
  list(A = a, B = b)
}

review_identity_annot <- function() {
  data.frame(GeneID = sprintf("g%02d", 1:24), Term = rep(paste0("P", 1:6), each = 4))
}

review_identity_labels <- function(p) {
  layers <- Filter(function(layer) "annotateText" %in% names(layer$data), p$layers)
  as.data.frame(layers[[1]]$data)
}

test_that("joined leading edges preserve exact spaced identifiers and literal separators", {
  annot <- read.csv(test_path("fixtures", "gsea_annot.csv"), colClasses = "character")
  ranked <- read.csv(test_path("fixtures", "gsea_ranked_a.csv"))
  annot$GeneID <- sub("G", "G ", annot$GeneID, fixed = TRUE)
  scores <- stats::setNames(ranked$score, sub("G", "G ", ranked$GeneID, fixed = TRUE))
  run <- function(table, sep) suppressMessages(richGSEA(
    scores, annot, seed = 42, nproc = 1, minSize = 3, maxSize = 20,
    pvalue = 1, padj = NULL, table = table, sep = sep))
  listed <- run(FALSE, "|")
  for (sep in c(",", "|", "::")) {
    joined <- run(TRUE, sep)
    expect_identical(joined@sep, sep)
    expect_identical(strsplit(joined@result$leadingEdge, sep, fixed = TRUE),
                     unname(listed@result$leadingEdge))
    expect_true(all(unlist(strsplit(joined@result$leadingEdge, sep, fixed = TRUE)) %in% names(scores)))
    expect_identical(joined@result[setdiff(names(joined@result), "leadingEdge")],
                     listed@result[setdiff(names(listed@result), "leadingEdge")])
    expect_identical(joined@gene, listed@gene)
    expect_identical(joined@input, listed@input)
  }
})

test_that("comparison rejects ambiguous explicit group names before reading results", {
  for (bad in list(c("A", "A"), c("A", ""), c("A", " "), c("A", NA_character_))) {
    expect_error(compareGSEA(stats::setNames(list(NULL, NULL), bad), NULL),
                 "x has (duplicate|missing or blank) name")
  }
})

test_that("comparison keeps original group statistics and positional optional genes", {
  groups <- review_identity_groups()
  original <- groups
  for (unnamed in c(FALSE, TRUE)) {
    x <- if (unnamed) unname(groups) else groups[c("B", "A")]
    expected_names <- if (unnamed) c("Group 1", "Group 2") else c("B", "A")
    for (genes in list(NULL, lapply(x, function(g) stats::setNames(g@input, g@gene)),
                       unname(lapply(x, function(g) stats::setNames(g@input, g@gene))),
                       stats::setNames(lapply(x, function(g) stats::setNames(g@input, g@gene)),
                                       c("same", "same")))) {
      p <- compareGSEA(x, review_identity_annot(), gene = genes,
                       pvalue = 1, padj = NULL, label = TRUE)
      labels <- review_identity_labels(p)
      expect_setequal(labels$group, expected_names)
      for (i in seq_along(x)) {
        rows <- labels[labels$group == expected_names[i], ]
        expect_equal(rows$NES, x[[i]]@result$NES[match(rows$pathway, x[[i]]@result$pathway)])
      }
    }
  }
  expect_identical(groups, original)
})
