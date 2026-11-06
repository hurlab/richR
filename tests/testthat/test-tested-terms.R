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

# ---- GSEA (AC-003) -------------------------------------------------------------

# Pathway sizes after intersecting with the ranked genes, computed without
# fgsea: the family fgsea tests under [minSize, maxSize].
gsea_sizes <- function(annot, ranked) {
  inside <- annot[annot$GeneID %in% names(ranked), ]
  sizes <- table(inside$Term)
  stats::setNames(as.integer(sizes), names(sizes))
}

check_gsea_tested <- function(g, sizes, maxSize) {
  tt <- testedTerms(g)
  family <- names(sizes)[sizes >= 3 & sizes <= maxSize]
  expect_setequal(tt$pathway, family)
  expect_identical(tt$size, unname(sizes[tt$pathway]))
  expect_identical(tt$padj, stats::p.adjust(tt$pval, "BH"))
  expect_identical(tt$reported, tt$pathway %in% g@result$pathway)
  expect_gt(nrow(g@result), 0)
  expect_true(any(!tt$reported))
  expect_identical(names(tt), c(names(g@result), "reported"))
  expect_identical(tt[tt$reported, names(g@result), drop = FALSE], g@result)
  info <- enrichmentInfo(g)
  expect_identical(info$schema_version, 2L)
  expect_identical(info$testing$zero_overlap_count, 0L)
  expect_identical(nrow(tt), info$testing$tested_count)
  tt
}

test_that("GSEA records every pathway fgsea returned and flags the reported ones", {
  # Catches GSEA keeping only significant pathways, a tested table built from
  # every annotation pathway instead of what fgsea tested, a tested padj that
  # differs from the one the reported rows used, and tested rows missing the
  # columns @result carries.
  skip_if_not_installed("fgsea")
  annot <- .char_read("gsea_annot.csv", character = TRUE)
  sizes <- gsea_sizes(annot, .fgsea_ranked())
  expect_identical(sizes[c("T1", "T2", "T3", "T4", "T5")],
                   c(T1 = 8L, T2 = 8L, T3 = 7L, T4 = 5L, T5 = 10L))

  g <- .fgsea_run(minSize = 3, maxSize = 9, pvalue = 0.05, table = TRUE)
  tt <- check_gsea_tested(g, sizes, maxSize = 9)
  expect_false("T5" %in% tt$pathway)
  expect_type(tt$leadingEdge, "character")

  gl <- .fgsea_run(minSize = 3, maxSize = 9, pvalue = 0.05, table = FALSE)
  ttl <- check_gsea_tested(gl, sizes, maxSize = 9)
  expect_type(ttl$leadingEdge, "list")
})

test_that("GSEA with KEGG = TRUE records the KEGG level columns in the tested table", {
  # Catches tested rows missing the KEGG level columns richGSEA adds to
  # @result after filtering.
  skip_if_not_installed("fgsea")
  annot <- .char_read("gsea_annot.csv", character = TRUE)
  ids <- c(T1 = "00010", T2 = "00061", T3 = "00190", T4 = "00230", T5 = "00970")
  annot$Term <- unname(ids[annot$Term])
  sizes <- gsea_sizes(annot, .fgsea_ranked())

  g <- .fgsea_run(annot = annot, minSize = 3, maxSize = 9, pvalue = 0.05,
                  KEGG = TRUE)
  tt <- check_gsea_tested(g, sizes, maxSize = 9)
  expect_true(all(c("Level3", "Level2", "Level1") %in% names(tt)))
  expect_false(anyNA(tt$Level1))
})

# ---- derived views (AC-005 verb leg) -------------------------------------------

test_that("dplyr verbs on a richResult change @result only and keep the analysis record", {
  # Catches a verb that also prunes the tested table (losing the "tested, not
  # significant" record) or rewrites the analysis options in the metadata.
  r <- .fora_reference(pvalue = 1)
  expect_gt(nrow(r@result), 1)
  tested <- testedTerms(r)
  info <- enrichmentInfo(r)
  plain <- r@result
  views <- list(
    filter = list(dplyr::filter(r, Annot != "T1"),
                  dplyr::filter(plain, Annot != "T1")),
    arrange = list(dplyr::arrange(r, dplyr::desc(Pvalue)),
                   dplyr::arrange(plain, dplyr::desc(Pvalue))),
    slice = list(dplyr::slice(r, 0), dplyr::slice(plain, 0)),
    mutate = list(dplyr::mutate(r, Pvalue = 1), dplyr::mutate(plain, Pvalue = 1)),
    select = list(dplyr::select(r, Annot, Term, Pvalue),
                  dplyr::select(plain, Annot, Term, Pvalue)),
    rename = list(dplyr::rename(r, P = Pvalue), dplyr::rename(plain, P = Pvalue)),
    # The method's .drop default is FALSE.
    group_by = list(dplyr::group_by(r, Term),
                    dplyr::group_by(plain, Term, .drop = FALSE))
  )
  for (verb in names(views)) {
    v <- views[[verb]][[1]]
    expect_true(methods::is(v, "richResult"), info = verb)
    expect_identical(v@result, views[[verb]][[2]], info = verb)
    expect_identical(testedTerms(v), tested, info = verb)
    expect_identical(enrichmentInfo(v), info, info = verb)
  }
  # The verbs did change @result, so the tested table is not just a copy of it.
  expect_identical(nrow(views$filter[[1]]@result), nrow(plain) - 1L)
  expect_true(inherits(views$group_by[[1]]@result, "grouped_df"))
})

test_that("dplyr verbs on a GSEAResult change @result only and keep the analysis record", {
  skip_if_not_installed("fgsea")
  g <- .fgsea_run(minSize = 3, maxSize = 20, pvalue = 1, padj = NULL)
  expect_gt(nrow(g@result), 1)
  tested <- testedTerms(g)
  info <- enrichmentInfo(g)
  for (v in list(dplyr::filter(g, pval <= 0.5), dplyr::slice(g, 1),
                 dplyr::group_by(g, Term))) {
    expect_true(methods::is(v, "GSEAResult"))
    expect_identical(testedTerms(v), tested)
    expect_identical(enrichmentInfo(v), info)
  }
})
