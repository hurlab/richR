# Live providers (SPEC-RICHR-003 REQ-022): the default suite mocks KEGG REST and
# msigdbr; these tests run the real services, so a change in what they return
# (and so a drift of the mocks) is caught. They reach the network and run only
# with RICHR_LIVE_TESTS=true.

skip_unless_live <- function() {
  skip_if_not(identical(Sys.getenv("RICHR_LIVE_TESTS"), "true"),
              "live provider tests run only with RICHR_LIVE_TESTS=true")
}

test_that("live KEGG labels and version have the shapes richR parses", {
  skip_unless_live()
  labels <- richR:::.kegg_labels(builtin = FALSE, query_release = TRUE)
  expect_identical(labels$source, "KEGGREST")
  expect_null(labels$reason)
  expect_true(all(grepl("^[0-9]{5}$", rownames(labels$data))))
  expect_false(is.na(labels$data["04148", "annotation"]))
  # The KEGG version is the stated release or the stated update date (REQ-014).
  expect_match(labels$release,
               "^(Release .+|last update [0-9]{4}-[0-9]{2}-[0-9]{2})$")
})

test_that("a live KEGG annotation records KEGGREST and the KEGG version", {
  skip_unless_live()
  skip_if_not_installed("org.Hs.eg.db")
  kegg <- suppressMessages(richR::buildAnnot(species = "human", keytype = "ENTREZID",
                                             anntype = "KEGG", builtin = FALSE))
  expect_true(nrow(kegg@annot) > 1000L)
  expect_true(all(grepl("^[0-9]{5}$", kegg@annot$PATH)))
  info <- richR::enrichmentInfo(kegg)$provenance$annotation
  expect_identical(info$label_source$source, "KEGGREST")
  expect_identical(info$source_version, info$label_source$release)
  expect_match(info$source_version,
               "^(Release .+|last update [0-9]{4}-[0-9]{2}-[0-9]{2})$")
})

test_that("a live MSigDB annotation records the msigdbr database version", {
  skip_unless_live()
  hallmark <- suppressMessages(utils::capture.output(
    annot <- richR::buildMSIGDB(species = "human", keytype = "SYMBOL", anntype = "HALLMARK")))
  expect_true(nrow(annot@annot) > 1000L)
  expect_identical(as.character(annot@annot$Term), as.character(annot@annot$Annot))
  expect_true(length(unique(annot@annot$Term)) >= 40L)
  info <- richR::enrichmentInfo(annot)$provenance$annotation
  expect_identical(info$provider_package, "msigdbr")
  expect_false(is.na(info$source_version))
})
