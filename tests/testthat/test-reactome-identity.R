# Reactome annotation identity (REQ-006 .. REQ-008), offline: the reactome.db
# maps are replaced by small in-memory maps through the .reactome_map accessor,
# so these tests need neither reactome.db nor the network. Upstream-portable:
# they assert only the annotation frame and ORA gene-set keys, never result
# labels (the fork-only label checks live in test-reactome-provider.R).

reactome_maps <- function(name) {
  maps <- list(
    reactomePATHID2EXTID = list(
      "R-HSA-1" = c("101", "102", "103", "104", "105"),
      "R-HSA-2" = c("106", "107", "108", "109", "110", "111", "112"),
      "R-HSA-5" = c("120", "121", "122"),
      "R-MMU-4" = c("201", "202", "203")
    ),
    reactomePATHID2NAME = list(
      "R-HSA-1" = "Homo sapiens: Apoptosis",
      "R-HSA-2" = "Homo sapiens: Apoptosis",
      "R-HSA-3" = "Homo sapiens: Pathway without genes",
      "R-HSA-5" = "Homo sapiens: Signaling by WNT",
      "R-MMU-4" = "Mus musculus: Apoptosis"
    )
  )
  maps[[name]]
}

build_reactome <- function() {
  testthat::with_mocked_bindings(
    suppressMessages(buildAnnot(species = "human", keytype = "ENTREZID",
                                anntype = "Reactome")),
    .reactome_map = reactome_maps, .package = "richR"
  )
}

test_that("Reactome annotation keeps the pathway identifier and the plain name", {
  # Catches: the identifier dropped from the Reactome annotation.
  annot <- build_reactome()@annot
  expect_identical(colnames(annot), c("GeneID", "Term", "Annot"))
  expect_true(all(grepl("^R-HSA-[0-9]+$", annot$Term)))
  expect_false(any(grepl("Homo sapiens: ", annot$Annot, fixed = TRUE)))
  expect_identical(annot$Annot[annot$Term == "R-HSA-5"][1], "Signaling by WNT")
})

test_that("Reactome pathways that share a name stay separate identifiers", {
  # Catches: Reactome identifiers lost to a shared name, a gene-less or
  # other-species pathway leaking into the annotation.
  annot <- build_reactome()@annot
  expect_identical(sort(unique(annot$Term)), c("R-HSA-1", "R-HSA-2", "R-HSA-5"))
  expect_identical(sort(annot$GeneID[annot$Term == "R-HSA-1"]),
                   c("101", "102", "103", "104", "105"))
  expect_identical(length(annot$GeneID[annot$Term == "R-HSA-2"]), 7L)
})

test_that("enrich() tests same-name Reactome pathways as two gene sets", {
  # Catches: same-name Reactome pathways merged into one gene set in ORA.
  annot <- build_reactome()
  res <- suppressMessages(enrich(c("101", "102", "106", "107", "120"), annot,
                                 pvalue = 1, padj = 1, minSize = 1, minGSSize = 1))@result
  apoptosis <- res[res$Annot %in% c("R-HSA-1", "R-HSA-2"), ]
  expect_identical(sort(apoptosis$Annot), c("R-HSA-1", "R-HSA-2"))
  expect_equal(apoptosis$Annotated[match(c("R-HSA-1", "R-HSA-2"), apoptosis$Annot)],
               c(5, 7))
})

test_that("a Reactome identifier with several names takes the first and warns once", {
  # Catches: a silent pick among an identifier's names (REQ-003 rule).
  several_names <- function(name) {
    map <- reactome_maps(name)
    if (name == "reactomePATHID2EXTID") map[["R-HSA-7"]] <- c("130", "131")
    if (name == "reactomePATHID2NAME") {
      map[["R-HSA-7"]] <- c("Homo sapiens: First name", "Homo sapiens: Second name")
    }
    map
  }
  msgs <- character()
  annot <- withCallingHandlers(
    testthat::with_mocked_bindings(
      suppressMessages(buildAnnot(species = "human", keytype = "ENTREZID",
                                  anntype = "Reactome")),
      .reactome_map = several_names, .package = "richR"
    ),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )@annot
  expect_identical(unique(annot$Annot[annot$Term == "R-HSA-7"]), "First name")
  named <- grep("R-HSA-7", msgs, value = TRUE)
  expect_identical(length(named), 1L)
  expect_false(any(grepl("R-HSA-1|R-HSA-2|R-HSA-5", named)))
})
