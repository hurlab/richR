# Reactome identity against the real reactome.db provider (Suggests; local, no
# network). Guards the in-memory maps of test-reactome-identity.R against
# drifting from the shape of the real maps. reactome.db is loaded, never attached.

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
  id2gene <- AnnotationDbi::as.list(reactome.db::reactomePATHID2EXTID)
  id2name <- AnnotationDbi::as.list(reactome.db::reactomePATHID2NAME)
  first_name <- vapply(id2name, function(x) as.character(x)[1], character(1))
  human <- names(first_name)[startsWith(first_name, "Homo sapiens: ")]
  with_genes <- human[human %in% names(id2gene)]
  with_genes <- with_genes[vapply(id2gene[with_genes],
                                  function(g) any(!is.na(g)), logical(1))]
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
})
