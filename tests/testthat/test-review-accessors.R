# Regression controls for SPEC-RICHR-013; no shared fixture state.
.review_access_result <- function() {
  methods::new("richResult",
    result = data.frame(
      Annot = c("A", "B", "C", "D"),
      Term = c("B", "unique label", "shared label", "shared label"),
      Pvalue = c(0.001, 0.02, 0.03, 0.04),
      Padj = c(0.004, 0.04, 0.04, 0.04),
      Significant = c(3L, 2L, 1L, 1L), Annotated = c(4L, 4L, 5L, 7L),
      GeneID = c("G1,G2,G3", "G5,G6", "G7", "G8")),
    gene = paste0("G", c(1:3, 5:8)), genenumber = 7,
    pvalueCutoff = 0.05, padjCutoff = 0.1, pAdjustMethod = "BH",
    organism = "human", ontology = "custom", keytype = "SYMBOL", sep = ",",
    detail = data.frame(background = 20L), tested = data.frame(background = 20L))
}

test_that("ORA identifiers precede labels and label-only matches are unique", {
  x <- .review_access_result()
  expect_identical(getGenes(x, "B"), c("G5", "G6"))
  expect_identical(getGenes(x, "unique label"), c("G5", "G6"))
  expect_error(getGenes(x, "shared label"), "[Aa]mbiguous|several")
  expect_error(getGenes(x, "absent"), "not found")
  expect_identical(getGenes(x), x@gene)
})

for (delimiter in c("|", ".")) {
  test_that(paste("ORA and joined GSEA split literal delimiter", delimiter), {
    genes <- c("gene one", "gene+two", "gene_three")
    x <- .review_access_result()
    x@result$GeneID[2] <- paste(genes, collapse = delimiter)
    x@sep <- delimiter
    expect_identical(getGenes(x, "B", sep = delimiter), genes)
    # The public default stays comma, even when stored metadata says otherwise.
    expect_identical(getGenes(x, "B"), paste(genes, collapse = delimiter))
    g <- methods::new("GSEAResult", result = data.frame(
      pathway = "hsa00010", Term = "glycolysis",
      leadingEdge = paste(genes, collapse = delimiter)), gene = genes, sep = delimiter)
    expect_identical(getGenes(g, "00010", sep = delimiter), genes)
    expect_identical(getGenes(g, "glycolysis", sep = delimiter), genes)
    g@result$leadingEdge <- I(list(genes))
    expect_identical(getGenes(g, "hsa00010", sep = delimiter), genes)
  })
}

test_that("comma splitting and explicit overrides retain their existing contract", {
  x <- .review_access_result()
  expect_identical(getGenes(x, "A"), c("G1", "G2", "G3"))
  x@sep <- "|"
  expect_identical(getGenes(x, "A", sep = ","), c("G1", "G2", "G3"))
  g <- methods::new("GSEAResult", result = data.frame(
    pathway = "p1", Term = "one", leadingEdge = "G1,G2"))
  expect_identical(getGenes(g, "p1"), c("G1", "G2"))
})

for (columns in list(c("Annot", "Term"), "Term", "Pvalue", "Padj", character())) {
  test_that(paste("show handles valid richResult columns", paste(columns, collapse = ",")), {
    original <- .review_access_result()
    x <- dplyr::select(original, dplyr::all_of(columns))
    expect_true(methods::validObject(x))
    before <- x
    output <- capture.output(returned <- methods::show(x))
    expect_match(paste(output, collapse = "\n"), "richResult object")
    if ("Term" %in% columns) expect_match(paste(output, collapse = "\n"), "unique label")
    expect_identical(returned, before)
    expect_identical(x, before)
    for (slot in setdiff(methods::slotNames(original), "result")) {
      expect_identical(methods::slot(x, slot), methods::slot(original, slot))
    }
    expect_identical(original, .review_access_result())
  })
}
