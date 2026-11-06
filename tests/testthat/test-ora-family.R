# ORA tested family (SPEC-RICHR-002 D2, REQ-008 ... REQ-012).
# Expected values are computed here with stats::phyper and stats::p.adjust,
# never copied from richR output.

fam_sets <- list(
  S  = c("Q1", "X1"),                 # overlap 2, size 2 < minGSSize = 3
  T1 = c("Q1", "Q2", "Q3", "F1"),
  T2 = c("Q4", "Q5", "F2", "F3", "F4"),
  T4 = c("Q6", "X1", "F5", "F6"),
  Z1 = paste0("F", 8:15),             # size-eligible, zero overlap
  Z2 = paste0("F", 16:25),
  Z3 = paste0("F", 26:30)
)
fam_annot <- data.frame(
  GeneID = unlist(fam_sets, use.names = FALSE),
  Term = rep(names(fam_sets), lengths(fam_sets)),
  stringsAsFactors = FALSE
)
fam_query <- c(paste0("Q", 1:6), "X1")
# N and n are fixed before any term-size filtering: X1 belongs to the
# sub-minGSSize term S and to T4, and every annotated gene counts.
fam_N <- length(unique(fam_annot$GeneID))
fam_n <- length(fam_query)
fam_M <- lengths(fam_sets)
fam_k <- vapply(fam_sets, function(s) sum(s %in% fam_query), numeric(1))
fam_p <- ifelse(fam_k > 0,
                stats::phyper(fam_k - 1, fam_M, fam_N - fam_M, fam_n, lower.tail = FALSE),
                1)
fam_overlap <- c("T1", "T2", "T4")
fam_all <- c(fam_overlap, "Z1", "Z2", "Z3")
fam_bh_overlap <- stats::p.adjust(fam_p[fam_overlap], method = "BH")
fam_bh_all <- stats::p.adjust(fam_p[fam_all], method = "BH")

fam_expect_padj <- function(result, expected) {
  res <- result@result
  expect_true(nrow(res) > 0)
  expect_equal(unname(res$Padj), unname(expected[as.character(res$Annot)]),
               tolerance = 1e-12)
}

test_that("fixture discriminates the overlap family from all positive-overlap terms", {
  old_family <- stats::p.adjust(fam_p[c("S", fam_overlap)], method = "BH")
  expect_equal(fam_N, 36L)
  expect_equal(fam_n, 7L)
  expect_false(isTRUE(all.equal(unname(old_family[fam_overlap]), unname(fam_bh_overlap))))
  expect_false(isTRUE(all.equal(unname(fam_bh_all[fam_overlap]), unname(fam_bh_overlap))))
  expect_true(all(fam_bh_all[fam_overlap] < 1))
})

test_that("default family is size-bounded positive overlap with N and n fixed first", {
  res <- suppressMessages(richR::enrich(
    fam_query, fam_annot, minGSSize = 3, minSize = 1, pvalue = 1, padj = 1
  ))
  expect_setequal(res@result$Annot, fam_overlap)
  expect_false("S" %in% res@result$Annot)
  fam_expect_padj(res, fam_bh_overlap)
  ids <- as.character(res@result$Annot)
  k <- fam_k[ids]
  M <- fam_M[ids]
  expect_equal(unname(res@result$Pvalue), unname(fam_p[ids]), tolerance = 1e-12)
  expect_equal(unname(res@result$FoldEnrichment), unname(k / M * fam_N / fam_n),
               tolerance = 1e-12)
  mu <- M * fam_n / fam_N
  sigma <- mu * (fam_N - fam_n) * (fam_N - M) / fam_N / (fam_N - 1)
  expect_equal(unname(res@result$zscore), unname((k - mu) / sqrt(sigma)),
               tolerance = 1e-12)
  explicit <- suppressMessages(richR::enrich(
    fam_query, fam_annot, minGSSize = 3, minSize = 1, pvalue = 1, padj = 1,
    padj.family = "overlap"
  ))
  expect_equal(explicit@result, res@result)
})

test_that("maxGSSize bounds the family before adjustment", {
  res <- suppressMessages(richR::enrich(
    fam_query, fam_annot, minGSSize = 3, maxGSSize = 4, minSize = 1,
    pvalue = 1, padj = 1
  ))
  expected <- stats::p.adjust(fam_p[c("T1", "T4")], method = "BH")
  expect_setequal(res@result$Annot, c("T1", "T4"))
  fam_expect_padj(res, expected)
  expect_equal(richR::enrichmentInfo(res)$testing$tested_count, 2L)
})

test_that("reporting filters never change adjusted values", {
  variants <- list(
    list(minSize = 1, maxSize = 500, keepRich = TRUE, pvalue = 1, padj = 1),
    list(minSize = 1, maxSize = 500, keepRich = FALSE, pvalue = 1, padj = 1),
    list(minSize = 3, maxSize = 500, keepRich = TRUE, pvalue = 1, padj = 1),
    list(minSize = 1, maxSize = 2, keepRich = TRUE, pvalue = 1, padj = 1),
    list(minSize = 1, maxSize = 500, keepRich = TRUE, pvalue = 0.2, padj = NULL),
    list(minSize = 1, maxSize = 500, keepRich = TRUE, pvalue = 1, padj = 0.2)
  )
  seen <- character()
  for (v in variants) {
    res <- suppressMessages(do.call(richR::enrich, c(
      list(fam_query, fam_annot, minGSSize = 3), v
    )))
    fam_expect_padj(res, fam_bh_overlap)
    expect_equal(richR::enrichmentInfo(res)$testing$tested_count, 3L)
    seen <- union(seen, as.character(res@result$Annot))
  }
  expect_setequal(seen, fam_overlap)
})

test_that("padj.family = 'all' adjusts over every size-eligible term", {
  res <- suppressMessages(richR::enrich(
    fam_query, fam_annot, minGSSize = 3, minSize = 1, pvalue = 1, padj = 1,
    padj.family = "all"
  ))
  expect_setequal(res@result$Annot, fam_overlap)
  fam_expect_padj(res, fam_bh_all)
  expect_false(any(res@result$Significant == 0))
  info <- richR::enrichmentInfo(res)
  expect_equal(info$testing$family_policy, "size_bounded_all_eligible")
  expect_equal(info$testing$tested_count, length(fam_all))
  expect_equal(info$options$padj.family, "all")

  default <- suppressMessages(richR::enrich(
    fam_query, fam_annot, minGSSize = 3, minSize = 1, pvalue = 1, padj = 1
  ))
  info <- richR::enrichmentInfo(default)
  expect_equal(info$testing$family_policy, "size_bounded_positive_overlap")
  expect_equal(info$testing$tested_count, length(fam_overlap))
  expect_equal(info$options$padj.family, "overlap")
})

test_that("enrich(<Annot>) receives padj.family by name, as the data.frame method does", {
  # Catches: S4 dispatch handing padj.family positionally to another formal
  # of the Annot method, so padj.family = "all" on an Annot object silently
  # adjusted over the overlap family.
  annot <- methods::new("Annot", species = "x", anntype = "custom",
                        keytype = "SYMBOL", annot = fam_annot)
  for (family in c("overlap", "all")) {
    expected <- if (family == "all") fam_bh_all else fam_bh_overlap
    args <- list(minGSSize = 3, minSize = 1, pvalue = 1, padj = 1, padj.family = family)
    df <- suppressMessages(do.call(richR::enrich, c(list(fam_query, fam_annot), args)))
    an <- suppressMessages(do.call(richR::enrich, c(list(fam_query, annot), args)))
    fam_expect_padj(an, expected)
    expect_equal(richR::testedTerms(an)$Padj, richR::testedTerms(df)$Padj, tolerance = 1e-12)
    info <- richR::enrichmentInfo(an)
    expect_equal(info$options$padj.family, family, info = family)
    expect_equal(info$testing$family_policy,
                 if (family == "all") "size_bounded_all_eligible" else "size_bounded_positive_overlap",
                 info = family)
    expect_equal(info$testing$tested_count, length(expected), info = family)
  }
})

test_that("padj.family is validated", {
  expect_error(suppressMessages(richR::enrich(
    fam_query, fam_annot, minGSSize = 3, minSize = 1, padj.family = "bogus"
  )), "should be one of")
  expect_error(richR::richDAVID("1", padj.family = "all"), "unused argument")
  stats <- c(Q1 = 3, Q2 = 2, F1 = -1, F2 = -2)
  expect_error(richR::richGSEA(stats, fam_annot, minSize = 1, padj.family = "all"),
               "unused argument")
})

test_that("richGO family is ontology-specific and honours padj.family", {
  go_annot <- data.frame(
    GeneID = c(fam_annot$GeneID, "Q1", "Q2", "M1", "M2", "M3"),
    GOALL = c(fam_annot$Term, rep("MF1", 5)),
    Annot = c(paste("label", fam_annot$Term), rep("label MF1", 5)),
    ONTOLOGYALL = c(rep("BP", nrow(fam_annot)), rep("MF", 5)),
    stringsAsFactors = FALSE
  )
  for (family in c("overlap", "all")) {
    res <- suppressMessages(richR::richGO(
      fam_query, godata = go_annot, ontology = "BP", minGSSize = 3, minSize = 1,
      pvalue = 1, padj = 1, padj.family = family
    ))
    expected <- if (family == "all") fam_bh_all else fam_bh_overlap
    expect_setequal(res@result$Annot, fam_overlap)
    fam_expect_padj(res, expected)
    info <- richR::enrichmentInfo(res)
    expect_equal(info$testing$effective_N, fam_N)
    expect_equal(info$testing$tested_count, length(expected))
    expect_equal(info$options$padj.family, family)
  }
})

test_that("richKEGG and richLevel honour padj.family", {
  path_env <- new.env(parent = emptyenv())
  utils::data("path", package = "richR", envir = path_env)
  ko <- setNames(path_env$path$ko[seq_along(fam_sets)], names(fam_sets))
  level3 <- setNames(path_env$path$Level3[seq_along(fam_sets)], names(fam_sets))
  kegg_annot <- data.frame(GeneID = fam_annot$GeneID,
                           PATH = unname(ko[fam_annot$Term]),
                           stringsAsFactors = FALSE)
  for (family in c("overlap", "all")) {
    expected <- if (family == "all") fam_bh_all else fam_bh_overlap
    kegg <- suppressMessages(richR::richKEGG(
      fam_query, kodata = kegg_annot, minGSSize = 3, minSize = 1,
      pvalue = 1, padj = 1, padj.family = family
    ))
    expect_setequal(kegg@result$Annot, ko[fam_overlap])
    expect_equal(unname(kegg@result$Padj),
                 unname(expected[names(ko)[match(kegg@result$Annot, ko)]]), tolerance = 1e-12)
    expect_equal(richR::enrichmentInfo(kegg)$options$padj.family, family)

    level <- suppressMessages(richR::richLevel(
      fam_query, kodata = kegg_annot, level = "Level3", minGSSize = 3,
      minSize = 1, pvalue = 1, padj = 1, padj.family = family
    ))
    expect_setequal(level@result$Annot, level3[fam_overlap])
    expect_equal(unname(level@result$Padj),
                 unname(expected[names(level3)[match(level@result$Annot, level3)]]), tolerance = 1e-12)
    expect_equal(richR::enrichmentInfo(level)$testing$tested_count, length(expected))
  }
})

test_that("batchEnrich forwards padj.family", {
  batch <- suppressMessages(richR::batchEnrich(
    list(g = fam_query), fam_annot, verbose = FALSE, minGSSize = 3,
    minSize = 1, pvalue = 1, padj = 1, padj.family = "all"
  ))
  fam_expect_padj(batch$results$g, fam_bh_all)
  expect_equal(richR::enrichmentInfo(batch$results$g)$options$padj.family, "all")
})

test_that("GSEA metadata keeps the fgsea-returned family", {
  stats <- c(A = 2, B = 1, C = -1, D = -2)
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = c("T1", "T1", "T2", "T2"))
  result <- testthat::with_mocked_bindings(
    fgseaMultilevel = function(...) data.frame(
      pathway = c("T1", "T2"), pval = c(0.01, 0.04), NES = c(1, -1),
      size = c(2L, 2L), leadingEdge = I(list("A", "C"))
    ),
    richR::richGSEA(stats, annot, padj = 1, minSize = 1, table = FALSE),
    .package = "richR"
  )
  info <- richR::enrichmentInfo(result)
  expect_equal(info$testing$family_policy, "fgsea_returned_before_filtering")
  expect_null(info$options$padj.family)
})

test_that("keepRich exempts only minSize and never admits terms outside the family", {
  sets <- list(A = paste0("a", 1:3), B = "b1", C = paste0("c", 1:3),
               D = paste0("d", 1:6), Z = paste0("f", 1:12))
  annot <- data.frame(GeneID = unlist(sets, use.names = FALSE),
                      Term = rep(names(sets), lengths(sets)),
                      stringsAsFactors = FALSE)
  query <- c(paste0("a", 1:3), "b1", "c1", "c2", paste0("d", 1:4))
  N <- length(unique(annot$GeneID))
  n <- length(query)
  M <- lengths(sets)
  k <- vapply(sets, function(s) sum(s %in% query), numeric(1))
  p <- stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  expected <- stats::p.adjust(p[c("A", "C", "D")], method = "BH")
  run <- function(...) suppressMessages(richR::enrich(
    query, annot, minGSSize = 2, pvalue = 1, padj = 1, ...
  ))

  rich <- run(minSize = 4, keepRich = TRUE)
  plain <- run(minSize = 4, keepRich = FALSE)
  # maxSize = 2 needs minSize <= maxSize (validation), so minSize is 2 here.
  capped <- run(minSize = 2, maxSize = 2, keepRich = TRUE)

  expect_setequal(rich@result$Annot, c("A", "D"))
  expect_setequal(plain@result$Annot, "D")
  expect_false("A" %in% capped@result$Annot)
  for (res in list(rich, plain, capped)) {
    expect_false("B" %in% res@result$Annot)
    fam_expect_padj(res, expected)
  }
  expect_false("C" %in% rich@result$Annot)
  expect_false("C" %in% plain@result$Annot)

  default_sets <- list(L = paste0("l", 1:10), R = "r1", W = paste0("w", 1:15))
  default_annot <- data.frame(
    GeneID = unlist(default_sets, use.names = FALSE),
    Term = rep(names(default_sets), lengths(default_sets)),
    stringsAsFactors = FALSE
  )
  with_rich <- suppressMessages(richR::enrich(
    c("l1", "l2", "r1"), default_annot, keepRich = TRUE, pvalue = 1, padj = 1
  ))
  without_rich <- suppressMessages(richR::enrich(
    c("l1", "l2", "r1"), default_annot, keepRich = FALSE, pvalue = 1, padj = 1
  ))
  expect_equal(with_rich@result$Annot, without_rich@result$Annot)
  expect_equal(with_rich@result$Annot, "L")
})

test_that("default family reports a null singleton term (documented FDR caveat)", {
  universe_sets <- setNames(as.list(paste0("g", 1:100)), paste0("S", 1:100))
  annot <- data.frame(GeneID = unlist(universe_sets, use.names = FALSE),
                      Term = names(universe_sets), stringsAsFactors = FALSE)
  overlap <- suppressMessages(richR::enrich(
    "g1", annot, minGSSize = 1, minSize = 1, pvalue = 1, padj = 0.05
  ))
  expect_equal(overlap@result$Annot, "S1")
  expect_equal(overlap@result$Padj, 0.01, tolerance = 1e-12)
  all_family <- suppressMessages(richR::enrich(
    "g1", annot, minGSSize = 1, minSize = 1, pvalue = 1, padj = 0.05,
    padj.family = "all"
  ))
  expect_equal(nrow(all_family@result), 0L)
  expect_equal(richR::enrichmentInfo(all_family)$testing$tested_count, 100L)
})
