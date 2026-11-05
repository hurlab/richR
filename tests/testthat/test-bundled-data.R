# Bundled KEGG snapshots (SPEC-RICHR-003 REQ-019, REQ-020): the shipped
# kegg.db, module and path objects say when they were retrieved and from which
# KEGG version, agree with what their documentation states, and are internally
# consistent, so that richLevel, richKEGG and richGSEA attach one hierarchy row
# to each pathway.

bundled_object <- function(dataset, object) {
  env <- new.env(parent = emptyenv())
  utils::data(list = dataset, package = "richR", envir = env)
  env[[object]]
}

bundled_rd_text <- function(topic) {
  rd <- tools::Rd_db("richR")[[paste0(topic, ".Rd")]]
  paste(unlist(lapply(rd, as.character)), collapse = "")
}

bundled_snapshots <- list(
  list(dataset = "kegg", object = "kegg.db", topic = "kegg.db"),
  list(dataset = "module", object = "module", topic = "module"),
  list(dataset = "path", object = "path", topic = "path")
)

test_that("each refreshed snapshot carries the provenance its documentation states", {
  for (snapshot in bundled_snapshots) {
    obj <- bundled_object(snapshot$dataset, snapshot$object)
    source <- attr(obj, "richR_source", exact = TRUE)
    expect_true(is.list(source), info = snapshot$object)
    expect_true(all(c("retrieved", "release", "url") %in% names(source)),
                info = snapshot$object)
    expect_match(source$retrieved, "^[0-9]{4}-[0-9]{2}-[0-9]{2}$", info = snapshot$object)
    expect_match(source$url, "^https://rest\\.kegg\\.jp/", info = snapshot$object)

    rd <- bundled_rd_text(snapshot$topic)
    stated <- regmatches(rd, regexec(paste0(
      "Retrieval date \\(", gsub(".", "\\.", snapshot$object, fixed = TRUE),
      "\\): ([^;]+); KEGG version: ([^;]+);"), rd))[[1]]
    expect_length(stated, 3L)
    expect_identical(stated[2], source$retrieved, info = snapshot$object)
    expect_identical(stated[3], if (is.na(source$release)) "NA" else source$release,
                     info = snapshot$object)

    rows <- regmatches(rd, regexec("A data\\.frame with ([0-9,]+) rows", rd))[[1]]
    expect_length(rows, 2L)
    expect_identical(as.integer(gsub(",", "", rows[2])), nrow(obj), info = snapshot$object)
  }
})

test_that("bundled KEGG labels report the snapshot's retrieval date", {
  retrieved <- attr(bundled_object("kegg", "kegg.db"), "richR_source", exact = TRUE)$retrieved
  module_retrieved <- attr(bundled_object("module", "module"), "richR_source",
                           exact = TRUE)$retrieved
  expect_false(is.null(retrieved))

  kodata <- data.frame(GeneID = c("1", "2", "3"), PATH = c("00010", "00010", "00020"))
  builtin <- suppressMessages(richR::richKEGG(c("1", "2"), kodata, pvalue = 1, padj = 1,
                                              minSize = 1, minGSSize = 1, builtin = TRUE))
  expect_identical(richR::enrichmentInfo(builtin)$provenance$labels$retrieved, retrieved)

  # A live request that falls back reports the date of the snapshot it used.
  expect_warning(
    fallback <- testthat::with_mocked_bindings(
      .getdbname = function(...) "org.Hs.eg.db",
      .getspeices = function(...) "hsa",
      keggLink = function(...) c("hsa:1" = "path:hsa00010", "hsa:2" = "path:hsa00020"),
      keggList = function(...) stop("simulated KEGG outage"),
      keggInfo = function(...) stop("simulated KEGG outage"),
      richR::buildAnnot(species = "human", keytype = "ENTREZID", anntype = "KEGG",
                        builtin = FALSE),
      .package = "richR"
    ),
    "bundled"
  )
  label_source <- richR::enrichmentInfo(fallback)$provenance$annotation$label_source
  expect_identical(label_source$source, "bundled")
  expect_identical(label_source$retrieved, retrieved)

  modules <- testthat::with_mocked_bindings(
    .getdbname = function(...) "org.Hs.eg.db",
    .getspeices = function(...) "hsa",
    keggLink = function(...) c("hsa:1" = "md:hsa:M00001", "hsa:2" = "md:hsa:M00002"),
    richR::buildAnnot(species = "human", keytype = "ENTREZID", anntype = "KEGGM"),
    .package = "richR"
  )
  expect_identical(
    richR::enrichmentInfo(modules)$provenance$annotation$label_source$retrieved,
    module_retrieved)
})

test_that("bundled KEGG identifiers are unique, well formed and labelled", {
  blank <- function(x) is.na(x) | !nzchar(trimws(x))

  kegg <- bundled_object("kegg", "kegg.db")
  expect_false(anyDuplicated(rownames(kegg)) > 0L)
  expect_true(all(grepl("^[0-9]{5}$", rownames(kegg))))
  expect_false(any(blank(kegg$annotation)))

  module <- bundled_object("module", "module")
  expect_false(anyDuplicated(rownames(module)) > 0L)
  expect_true(all(grepl("^M[0-9]{5}$", rownames(module))))
  expect_false(any(blank(module$annotation)))

  path <- bundled_object("path", "path")
  expect_identical(path$ko[duplicated(path$ko)], character(0))
  expect_true(all(grepl("^[0-9]{5}$", path$ko)))
  expect_identical(rownames(path), path$ko)
  for (level in c("Level3", "Level2", "Level1")) {
    expect_false(any(blank(path[[level]])), info = level)
  }
})

test_that("richLevel counts a pathway under the one hierarchy row richKEGG attaches", {
  # 00541 sat under two BRITE rows in the earlier snapshot, so richLevel
  # counted its genes twice while richKEGG and richGSEA used only the first row.
  genes <- sprintf("G%02d", 1:10)
  kodata <- data.frame(GeneID = genes, PATH = rep(c("00541", "00010"), each = 5))
  query <- genes[1:3]

  level3 <- suppressMessages(richR::richLevel(query, kodata, level = "Level3", pvalue = 1,
                                              padj = 1, minSize = 1, minGSSize = 1))
  kegg <- suppressMessages(richR::richKEGG(query, kodata, pvalue = 1, padj = 1,
                                           minSize = 1, minGSSize = 1))
  expect_identical(nrow(level3@result), 1L)
  expect_identical(level3@result$Term, kegg@result$Level3[kegg@result$Annot == "00541"])
})
