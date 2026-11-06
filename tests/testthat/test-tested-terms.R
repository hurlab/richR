# The tested table (SPEC-RICHR-004 REQ-001 .. REQ-004): what an analysis
# tested, not only what it reported. Fixtures F-ORA and F-GSEA are built in
# helper-tested.R with the public API.

# ---- ORA (AC-001, AC-002) ------------------------------------------------------

test_that("ORA records every tested term with overlap > 0 and flags the reported ones", {
  # Catches tested-but-not-significant terms discarded, size- or minSize-
  # filtered terms reported as untested, and a reported flag that disagrees
  # with @result.
  r <- .fora_reference()
  expect_identical(r@result$Annot, "T1")

  tt <- testedTerms(r)
  .fora_roles(tt)
  expect_type(tt$reported, "logical")
  expect_identical(tt$Annot[tt$reported], "T1")
  expect_false(any(c("T3", "T4") %in% tt$Annot))

  shared <- names(r@result)
  expect_identical(names(tt), c(shared, "reported"))
  expect_identical(tt[tt$reported, shared, drop = FALSE], r@result)

  info <- enrichmentInfo(r)
  expect_identical(info$schema_version, 2L)
  expect_identical(info$testing$zero_overlap_count, 0L)
  expect_identical(info$testing$tested_count, 3L)
  expect_identical(nrow(tt) + info$testing$zero_overlap_count,
                   info$testing$tested_count)
})

test_that("ORA under padj.family = 'all' counts zero-overlap terms and keeps the family's Padj", {
  # Catches a tested table built after the filters, a zero-overlap count that
  # does not reconcile with the adjusted family, and a tested Padj recomputed
  # over the stored rows only.
  r <- .fora_reference(padj.family = "all")
  tt <- testedTerms(r)
  .fora_roles(tt)
  expect_identical(tt$Annot[tt$reported], "T1")

  info <- enrichmentInfo(r)
  expect_identical(info$testing$zero_overlap_count, 1L)
  expect_identical(info$testing$tested_count, 4L)
  expect_identical(nrow(tt) + info$testing$zero_overlap_count,
                   info$testing$tested_count)
  # T3 (overlap 0) enters the adjusted family with p = 1.
  expect_identical(tt$Padj,
                   stats::p.adjust(c(tt$Pvalue, 1), "BH")[seq_len(nrow(tt))])
})

test_that("richKEGG records its tested table with the KEGG level columns", {
  # Catches an ORA entry point that bypasses the tested record or drops the
  # extra columns its @result carries.
  k <- suppressMessages(richKEGG(
    .char_read("ora_query.csv")$GeneID,
    .char_read("kegg_annot.csv", character = TRUE),
    pvalue = 0.05, padj = NULL, minSize = 2, minGSSize = 3
  ))
  tt <- testedTerms(k)
  expect_gt(nrow(k@result), 0)
  expect_true(any(!tt$reported))
  expect_setequal(tt$Annot[tt$reported], k@result$Annot)
  expect_identical(names(tt), c(names(k@result), "reported"))
  expect_true(all(c("ko", "Level3", "Level2", "Level1") %in% names(tt)))
  expect_identical(nrow(tt), enrichmentInfo(k)$testing$tested_count)
})

test_that("an ORA that tested nothing records an empty table, not an unavailable one", {
  # Catches an empty analysis reported as tested-unavailable (NULL).
  r <- .fora_reference(minGSSize = 50)
  expect_identical(enrichmentInfo(r)$testing$tested_count, 0L)
  expect_silent(tt <- testedTerms(r))
  expect_true(is.data.frame(tt))
  expect_identical(nrow(tt), 0L)
  expect_true(all(c("Annot", "Term", "Pvalue", "Padj", "reported") %in% names(tt)))
})

test_that("testedTerms() rejects objects that are not enrichment results", {
  annot <- new("Annot", species = "x", anntype = "custom", keytype = "SYMBOL",
               annot = data.frame(GeneID = "A", Term = "T"))
  expect_error(testedTerms(annot), "richResult or GSEAResult")
  expect_error(testedTerms(data.frame(Annot = "T1")), "richResult or GSEAResult")
})
