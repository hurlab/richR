# universe_policy (SPEC-RICHR-005 REQ-017): an opt-in ORA option counting
# tested-but-unannotated universe genes in N and n. The default ("annotated")
# must stay identical to the omitted argument on every entry point. Fixtures:
# F-AUD inline, F-ORA-C from the committed characterization fixtures.

fc_read <- function(f) read.csv(test_path("fixtures", f), colClasses = "character")
fc <- function() {
  oq <- fc_read("ora_query.csv")$GeneID
  oa <- fc_read("ora_annot.csv")
  ka <- fc_read("kegg_annot.csv")
  go <- data.frame(GeneID = c(oa$GeneID, "UN1"), GOALL = c(oa$Term, "MFTERM"),
                   Annot = c(oa$Annot, "mf term"),
                   ONTOLOGYALL = c(rep("BP", nrow(oa)), "MF"))
  uo <- unique(c(oa$GeneID, oq, "UN1", "UN2", "UN3"))
  list(oq = oq, oa = oa, ka = ka, go = go, uo = uo)
}
as_annot <- function(d) new("Annot", species = "x", anntype = "custom", keytype = "SYMBOL", annot = d)

# Every ORA entry point that accepts universe, with the genes its annotation tests.
fc_calls <- function(f) {
  args <- list(pvalue = 1, padj = NULL, minSize = 1, minGSSize = 1, universe = f$uo)
  bp_genes <- f$go$GeneID[f$go$ONTOLOGYALL == "BP"]
  run <- function(fun, x, annot, extra = list()) {
    function(...) suppressMessages(suppressWarnings(
      do.call(fun, c(list(f$oq, annot), args, extra, list(...)))))
  }
  batch <- function(annot, fun, extra = list()) {
    function(...) suppressMessages(suppressWarnings(
      do.call(batchEnrich, c(list(list(s = f$oq), annot, fun = fun, verbose = FALSE),
                             args, extra, list(...)))))$results$s
  }
  list(
    "enrich(df)" = list(call = run(enrich, f$oq, f$oa), genes = f$oa$GeneID),
    "enrich(Annot)" = list(call = run(enrich, f$oq, as_annot(f$oa)), genes = f$oa$GeneID),
    "richGO(df)" = list(call = run(richGO, f$oq, f$go, list(ontology = "BP")), genes = bp_genes),
    "richGO(Annot)" = list(call = run(richGO, f$oq, as_annot(f$go), list(ontology = "BP")),
                           genes = bp_genes),
    "richKEGG(df)" = list(call = run(richKEGG, f$oq, f$ka), genes = f$ka$GeneID),
    "richKEGG(Annot)" = list(call = run(richKEGG, f$oq, as_annot(f$ka)), genes = f$ka$GeneID),
    "batchEnrich(enrich)" = list(call = batch(f$oa, "enrich"), genes = f$oa$GeneID),
    "batchEnrich(richGO)" = list(call = batch(f$go, "richGO", list(ontology = "BP")),
                                 genes = bp_genes),
    "batchEnrich(richKEGG)" = list(call = batch(f$ka, "richKEGG"), genes = f$ka$GeneID)
  )
}

tested_sizes <- function(r) {
  t <- testedTerms(r)
  t[order(t$Annot), c("Annot", "Annotated", "Significant")]
}

test_that("AC-016 (a, b): the default is unchanged and \"measured\" reaches every entry point", {
  # Catches: the option silently changing the default (any entry point),
  # "annotated" computed from the whole GO frame instead of the selected
  # ontology, the option ignored by a wrapper (richGO, richKEGG, batchEnrich,
  # an Annot method), and the option altering term sizes or the tested family.
  f <- fc()
  expect_identical(length(f$uo), 31L)
  for (nm in names(calls <- fc_calls(f))) {
    cl <- calls[[nm]]
    omitted <- cl$call()
    annotated <- cl$call(universe_policy = "annotated")
    expect_identical(omitted, annotated, info = nm)
    expect_identical(enrichmentInfo(omitted), enrichmentInfo(annotated), info = nm)
    io <- enrichmentInfo(omitted)
    expect_identical(io$options$universe_policy, "annotated", info = nm)
    expect_identical(io$testing$effective_N, sum(f$uo %in% cl$genes), info = nm)

    measured <- cl$call(universe_policy = "measured")
    im <- enrichmentInfo(measured)
    expect_identical(im$options$universe_policy, "measured", info = nm)
    expect_identical(im$testing$effective_N, 31L, info = nm)
    expect_identical(im$testing$tested_count, io$testing$tested_count, info = nm)
    expect_identical(tested_sizes(measured), tested_sizes(omitted), info = nm)
  }
  go_bp <- calls[["richGO(df)"]]$call()
  expect_identical(enrichmentInfo(go_bp)$testing$effective_N, 27L)
  b <- mappingAudit(go_bp)
  expect_identical(b$decision[b$role == "background" & b$id == "UN1"], "not_annotated")
})

test_that("AC-016 (c): \"measured\" counts an unannotated universe gene in N and n", {
  # Catches: N / n changed while metadata or the audit still describe the
  # annotated sets, and p-values that do not follow the hypergeometric model
  # with the stated N and n.
  a <- data.frame(GeneID = paste0("G", 1:20), Term = rep(paste0("T", 1:4), each = 5))
  qu <- c("G1", "G2", "G3", "G6", "U1", "X9", "G19")
  ua <- c(paste0("G", 1:18), "U1", "U1", NA)
  run <- function(...) enrich(qu, a, pvalue = 1, padj = NULL, minSize = 1, minGSSize = 1,
                              universe = ua, ...)
  oracle <- function(r, N, n) {
    x <- r@result
    stats::phyper(x$Significant - 1, x$Annotated, N - x$Annotated, n, lower.tail = FALSE)
  }
  d <- run()
  id <- enrichmentInfo(d)
  expect_identical(c(id$testing$effective_N, id$testing$effective_n), c(18L, 4L))
  expect_identical(d@result$Annot, c("T1", "T2"))
  expect_equal(d@result$Pvalue, oracle(d, 18, 4), tolerance = 0)
  expect_equal(d@result$Pvalue, c(0.0441176, 0.76634), tolerance = 1e-5)

  m <- run(universe_policy = "measured")
  im <- enrichmentInfo(m)
  expect_identical(c(im$testing$effective_N, im$testing$effective_n), c(19L, 5L))
  expect_identical(m@result$Annot, c("T1", "T2"))
  expect_identical(m@result$Annotated, d@result$Annotated)
  expect_identical(m@result$Significant, d@result$Significant)
  expect_equal(m@result$Pvalue, oracle(m, 19, 5), tolerance = 0)
  expect_equal(m@result$Pvalue, c(0.0843653, 0.827829), tolerance = 1e-5)
  tt <- testedTerms(m)
  expect_equal(m@result$Padj, stats::p.adjust(tt$Pvalue, "BH")[match(m@result$Annot, tt$Annot)],
               tolerance = 0)
  expect_identical(length(im$background$effective_ids), 19L)
  expect_true("U1" %in% im$background$effective_ids)
  expect_identical(im$background$unmapped_ids, "U1")
  expect_identical(im$input$effective_ids, c("G1", "G2", "G3", "G6", "U1"))
  expect_identical(im$input$unmapped_ids, c("U1", "X9"))

  au <- mappingAudit(m)
  q <- au[au$role == "query", ]
  b <- au[au$role == "background", ]
  expect_identical(sum(b$decision == "effective"), 19L)
  expect_identical(q$decision[q$id %in% "U1"], "effective")
  expect_match(q$reason[q$id %in% "U1"], "not annotated")
  expect_identical(q$decision[q$id %in% c("X9", "G19")], rep("outside_background", 2))
  expect_identical(length(unique(q$id[q$decision == "effective"])), im$testing$effective_n)
  expect_identical(sum(q$occurrences, na.rm = TRUE), im$input$original_count)
  expect_identical(im$input$original_count, 7L)
})

test_that("AC-016 (d): \"measured\" without a universe, other values and richLevel error", {
  # Catches: a measured set invented when no universe was given, a measured run
  # with no testable overlap returned as a result, and an unsupported entry
  # point accepting and dropping the option.
  a <- data.frame(GeneID = paste0("G", 1:20), Term = rep(paste0("T", 1:4), each = 5))
  qu <- c("G1", "G2", "G3", "G6", "U1", "X9", "G19")
  ua <- c(paste0("G", 1:18), "U1", "U1", NA)
  run <- function(q, ...) enrich(q, a, pvalue = 1, padj = NULL, minSize = 1, minGSSize = 1, ...)
  expect_error(run(qu, universe = NULL, universe_policy = "measured"),
               "universe_policy.*universe")
  expect_error(run(qu, universe = ua, universe_policy = "all"), "universe_policy.*all")
  expect_error(run(c("U1", "X9"), universe = ua, universe_policy = "measured"),
               "no query genes remain in the effective universe")
  f <- fc()
  expect_error(richLevel(f$oq, f$ka, universe_policy = "measured"),
               "unused argument \\(universe_policy")
})
