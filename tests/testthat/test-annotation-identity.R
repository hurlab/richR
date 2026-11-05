# Annotation identity (SPEC-RICHR-003 REQ-001 .. REQ-005).
#
# One identifier rule for every builder, ORA and GSEA: the gene is column 1;
# the identifier is the first of GOALL, PATH, Module, Term present (else column
# 2); the label is the column named Annot. ORA results report the identifier in
# `Annot` and the label in `Term`.

id_genes <- sprintf("g%02d", 1:40)

# Three gene sets S1..S3 with labels distinct from their identifiers.
id_sets <- function() {
  data.frame(
    GeneID = c(id_genes[1:6], id_genes[7:14], id_genes[30:37]),
    id = rep(c("S1", "S2", "S3"), c(6, 8, 8)),
    label = rep(c("set one", "set two", "set three"), c(6, 8, 8)),
    stringsAsFactors = FALSE
  )
}

id_query <- function() id_genes[c(1:4, 7:10, 30:33)]

id_ranks <- function() stats::setNames(seq(3, -3, length.out = length(id_genes)), id_genes)

id_ora <- function(annot, ...) {
  suppressMessages(enrich(id_query(), annot, pvalue = 1, padj = 1, minSize = 1,
                          minGSSize = 1, ...))
}

id_gsea <- function(annot) {
  out <- NULL
  utils::capture.output(
    out <- richGSEA(id_ranks(), annot, minSize = 1, maxSize = 50, pvalue = 1,
                    padj = 1, seed = 7, nproc = 1)
  )
  out
}

# The four layouts of AC-001, all carrying the same memberships.
id_layouts <- function() {
  s <- id_sets()
  list(
    label_first = data.frame(GeneID = s$GeneID, Annot = s$label, Term = s$id,
                             stringsAsFactors = FALSE),
    id_first = data.frame(GeneID = s$GeneID, Term = s$id, Annot = s$label,
                          stringsAsFactors = FALSE),
    go = data.frame(GeneID = s$GeneID, GOALL = s$id, ONTOLOGYALL = "BP",
                    Annot = s$label, stringsAsFactors = FALSE),
    unnamed = data.frame(V1 = s$GeneID, V2 = s$id, stringsAsFactors = FALSE)
  )
}

test_that("ORA keys every annotation layout by its identifier, never by its label", {
  # Catches: enrich() keying a (GeneID, Annot, Term) frame by the labels in column 2.
  for (layout in names(id_layouts())) {
    res <- id_ora(id_layouts()[[layout]])@result
    expect_setequal(res$Annot, c("S1", "S2", "S3"))
    expect_identical(res$Annotated[match(c("S1", "S2", "S3"), res$Annot)], c(6L, 8L, 8L))
    expected_terms <- if (layout == "unnamed") res$Annot else
      c(S1 = "set one", S2 = "set two", S3 = "set three")[res$Annot]
    expect_identical(res$Term, unname(expected_terms), info = layout)
  }
})

test_that("ORA and GSEA report the same identifier set for every layout", {
  skip_if_not_installed("fgsea")
  # Catches: the same frame keyed differently by ORA and GSEA.
  for (layout in names(id_layouts())) {
    annot <- id_layouts()[[layout]]
    ora_ids <- id_ora(annot)@result$Annot
    gsea_ids <- id_gsea(annot)@result$pathway
    expect_setequal(ora_ids, c("S1", "S2", "S3"))
    expect_setequal(gsea_ids, ora_ids)
  }
})

conflict_annot <- function() {
  data.frame(
    GeneID = c(id_genes[1:6], id_genes[7:14]),
    Term = rep(c("T1", "T2"), c(6, 8)),
    Annot = c("a", "a", "b", "b", "a", "b", rep("fine", 8)),
    stringsAsFactors = FALSE
  )
}

conflict_warnings <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, conflicts = grep("conflicting labels", msgs, value = TRUE))
}

test_that("ORA warns once on conflicting labels and keeps the first in row order", {
  # Catches: a silent pick among conflicting labels.
  out <- conflict_warnings(id_ora(conflict_annot()))
  expect_length(out$conflicts, 1L)
  expect_match(out$conflicts, "T1")
  expect_false(grepl("T2", out$conflicts))
  # The message names the column convention, so a frame whose identifier and
  # label columns are swapped can be recognised from it.
  expect_match(out$conflicts, "Annot", fixed = TRUE)
  expect_match(out$conflicts, "Term (or GOALL, PATH, Module", fixed = TRUE)
  res <- out$value@result
  expect_identical(res$Term[res$Annot == "T1"], "a")
  expect_identical(res$Term[res$Annot == "T2"], "fine")
})

test_that("GSEA warns once on the same conflicting labels and keeps the same label", {
  skip_if_not_installed("fgsea")
  out <- conflict_warnings(id_gsea(conflict_annot()))
  expect_length(out$conflicts, 1L)
  expect_match(out$conflicts, "T1")
  res <- out$value@result
  expect_identical(res$Term[res$pathway == "T1"], "a")
})

test_that("buildAnnot counts and enrich keys GO sets by GOALL, not by column position", {
  # Catches: one ontology code reported (and tested) as a single gene set.
  go <- data.frame(
    GeneID = id_genes[1:9], ONTOLOGYALL = "BP",
    GOALL = rep(c("GO:0000001", "GO:0000002", "GO:0000003"), each = 3),
    Annot = rep(c("alpha", "beta", "gamma"), each = 3),
    stringsAsFactors = FALSE
  )
  msgs <- character()
  annot <- withCallingHandlers(
    testthat::with_mocked_bindings(
      buildAnnot(species = "human", keytype = "SYMBOL", anntype = "GO"),
      .makeGOdata = function(...) go, .package = "richR"
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("9 genes, 3 terms", msgs)))
  expect_identical(colnames(annot@annot), colnames(go))
  res <- suppressMessages(enrich(id_genes[c(1, 4, 7)], annot, pvalue = 1, padj = 1,
                                 minSize = 1, minGSSize = 1))@result
  expect_setequal(res$Annot, c("GO:0000001", "GO:0000002", "GO:0000003"))
  expect_identical(res$Term[match("GO:0000002", res$Annot)], "beta")
})

test_that("richGO keys by GOALL when column 2 is not the identifier", {
  # Catches: richGO testing the ontology code as one gene set.
  ids <- c("GO:0006915", "GO:0007049")
  go <- data.frame(GeneID = id_genes[1:10], ONTOLOGYALL = "BP",
                   GOALL = rep(ids, each = 5), stringsAsFactors = FALSE)
  res <- suppressMessages(richGO(id_genes[c(1:3, 6:8)], go, ontology = "BP",
                                 pvalue = 1, padj = 1, minSize = 1, minGSSize = 1))@result
  expect_setequal(res$Annot, ids)
  expect_identical(res$Term[match(ids, res$Annot)],
                   unname(AnnotationDbi::Term(GO.db::GOTERM[ids])))
})

test_that("richKEGG takes labels from the KEGG term list, then the frame Annot, then the ID", {
  # Catches: an identifier shown where a label exists.
  kegg <- data.frame(
    GeneID = id_genes[1:9],
    PATH = rep(c("00010", "99999", "99998"), each = 3),
    Annot = rep(c("frame label for 00010", "Custom pathway", ""), each = 3),
    stringsAsFactors = FALSE
  )
  res <- suppressMessages(richKEGG(id_genes[c(1, 4, 7)], kegg, pvalue = 1, padj = 1,
                                   minSize = 1, minGSSize = 1))@result
  expect_identical(res$Term[res$Annot == "00010"], "Glycolysis / Gluconeogenesis")
  expect_identical(res$Term[res$Annot == "99999"], "Custom pathway")
  expect_identical(res$Term[res$Annot == "99998"], "99998")
})
