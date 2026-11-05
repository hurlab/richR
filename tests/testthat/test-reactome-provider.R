# Reactome identity against the real reactome.db provider (Suggests; local, no
# network). Guards the in-memory maps of test-reactome-identity.R against
# drifting from the shape of the real maps. reactome.db is loaded, never attached.
# Fork-only: also holds the result-label checks (ORA Term, GSEA Term) that the
# upstream-portable test-reactome-identity.R deliberately leaves out.

reactome_provider_annot <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- suppressMessages(buildAnnot(species = "human", keytype = "ENTREZID",
                                            anntype = "Reactome"))
    }
    cache
  }
})

test_that("buildAnnot(Reactome) works from the real provider without attaching it", {
  skip_if_not_installed("reactome.db")
  expect_false("package:reactome.db" %in% search())
  annot <- reactome_provider_annot()@annot
  expect_identical(colnames(annot), c("GeneID", "Term", "Annot"))
  expect_true(all(grepl("^R-HSA-[0-9]+$", annot$Term)))
  expect_false(any(startsWith(annot$Annot, "Homo sapiens: ")))
  expect_false("package:reactome.db" %in% search())
})

test_that("every human pathway with a gene is its own identifier", {
  skip_if_not_installed("reactome.db")
  # Catches: real Reactome identifiers lost to a shared pathway name.
  # Independent oracle: human pathways by their stable-identifier prefix
  # (R-HSA-), not by the species name the implementation filters on.
  id2gene <- AnnotationDbi::as.list(reactome.db::reactomePATHID2EXTID)
  human <- grep("^R-HSA-[0-9]+$", names(id2gene), value = TRUE)
  with_genes <- human[vapply(id2gene[human], function(g) any(!is.na(g)), logical(1))]
  terms <- unique(reactome_provider_annot()@annot$Term)
  expect_identical(sort(terms), sort(unique(with_genes)))
  if (identical(as.character(utils::packageVersion("reactome.db")), "1.96.0")) {
    # R-1 (SPEC-RICHR-003): 2845 human pathways with genes, 18 of them sharing 9 names.
    expect_identical(length(terms), 2845L)
  }
})

test_that("real same-name Reactome pathways stay distinct in ORA and GSEA", {
  skip_if_not_installed("reactome.db")
  # Catches: real same-name Reactome pathways merged into one gene set.
  annot <- reactome_provider_annot()
  frame <- annot@annot
  labels <- unique(frame[, c("Term", "Annot")])
  shared_name <- labels$Annot[duplicated(labels$Annot)][1]
  ids <- sort(labels$Term[labels$Annot == shared_name])
  expect_true(length(ids) >= 2L)
  sub <- frame[frame$Term %in% ids, ]
  ora <- suppressWarnings(suppressMessages(enrich(
    unique(sub$GeneID), sub, pvalue = 1, padj = 1, minSize = 1, minGSSize = 1
  )))@result
  expect_identical(sort(ora$Annot), ids)
  expect_true(all(ora$Term == shared_name))
  skip_if_not_installed("fgsea")
  genes <- unique(c(sub$GeneID, sprintf("bg%03d", 1:40)))
  ranks <- stats::setNames(seq(3, -3, length.out = length(genes)), genes)
  gsea <- NULL
  utils::capture.output(
    gsea <- richGSEA(ranks, sub, minSize = 1, maxSize = 1000, pvalue = 1, padj = 1,
                     seed = 7, nproc = 1)@result
  )
  expect_identical(sort(gsea$pathway), ids)
  expect_true(all(gsea$Term == shared_name))
})

test_that("same-name Reactome pathways keep their identifier and name in ORA and GSEA results", {
  # Fork-only labelling (ORA labels from SPEC-RICHR-003 REQ-002, GSEA identity
  # from SPEC-RICHR-002); offline, with in-memory maps, so it runs without
  # reactome.db. Catches: same-name pathways merged, or results showing the
  # identifier where the pathway name exists.
  maps <- list(
    reactomePATHID2EXTID = list("R-HSA-1" = as.character(101:105),
                                "R-HSA-2" = as.character(106:112)),
    reactomePATHID2NAME = list("R-HSA-1" = "Homo sapiens: Apoptosis",
                               "R-HSA-2" = "Homo sapiens: Apoptosis")
  )
  annot <- testthat::with_mocked_bindings(
    suppressMessages(buildAnnot(species = "human", keytype = "ENTREZID",
                                anntype = "Reactome")),
    .reactome_map = function(name) maps[[name]], .package = "richR"
  )
  ora <- suppressMessages(enrich(c("101", "102", "106", "107"), annot, pvalue = 1,
                                 padj = 1, minSize = 1, minGSSize = 1))@result
  expect_identical(sort(ora$Annot), c("R-HSA-1", "R-HSA-2"))
  expect_identical(ora$Term, c("Apoptosis", "Apoptosis"))
  skip_if_not_installed("fgsea")
  genes <- c(as.character(101:112), sprintf("bg%02d", 1:30))
  ranks <- stats::setNames(seq(3, -3, length.out = length(genes)), genes)
  gsea <- NULL
  utils::capture.output(
    gsea <- richGSEA(ranks, annot, minSize = 1, maxSize = 50, pvalue = 1, padj = 1,
                     seed = 7, nproc = 1)@result
  )
  expect_identical(sort(gsea$pathway), c("R-HSA-1", "R-HSA-2"))
  expect_identical(gsea$Term, c("Apoptosis", "Apoptosis"))
})
