# writeRich() / readRich() (SPEC-RICHR-004 REQ-011..REQ-014, AC-005 I/O leg,
# AC-010..AC-014): a result object stored as a folder of TSV tables plus a JSON
# manifest. Fixtures are built with the public API (helper-tested.R) and every
# folder is written under tempdir().

io_dir <- function() tempfile("rich-io-")

# F-S22 (AC-010): 22 gene sets with five-digit identifiers holding leading
# zeros, 4 genes each in a 100-gene universe; the query overlaps 18 sets in 2
# genes and 4 sets in 1 gene, so enrich() reports 18 rows and records 4
# tested-but-not-reported rows.
fs22 <- function() {
  u <- sprintf("g%03d", 1:100)
  annot <- do.call(rbind, lapply(1:22, function(i) {
    data.frame(GeneID = u[(4 * i - 3):(4 * i)], Term = sprintf("%05d", 10 * i),
               Annot = paste("set", i), stringsAsFactors = FALSE)
  }))
  query <- c(unlist(lapply(1:18, function(i) u[(4 * i - 3):(4 * i - 2)])),
             u[4 * (19:22) - 3])
  suppressMessages(richR::enrich(query, annot, universe = u, pvalue = 1,
                                 padj = NULL, minSize = 2, minGSSize = 1))
}

# The 17 character values of acceptance.md L-P4.
lp4_values <- c(NA, "", "NA", "a\tb", "l1\nl2", "c\rr", "q\"uote", "it's",
                "back\\slash", "\\N", "\\\\N", "#hash", " lead", "trail ",
                "日本", "00010", "\\t")

# F-S22 edited by slot assignment only: special doubles, the L-P4 labels in
# @result and the same labels in @tested for the same identifiers, and a
# factor column in @detail.
fs22_edited <- function() {
  x <- fs22()
  x@result$zscore[1:6] <- c(NA, NaN, Inf, -Inf, 5e-324, 1 / 3)
  x@result$Term[1:17] <- lp4_values
  x@tested$Term[match(x@result$Annot[1:17], x@tested$Annot)] <- lp4_values
  x@detail$flag <- factor(rep_len(c("b", "a", NA), nrow(x@detail)),
                          levels = c("b", "a", "z"))
  x
}

kegg_result <- function() {
  suppressMessages(richKEGG(
    .char_read("ora_query.csv")$GeneID,
    .char_read("kegg_annot.csv", character = TRUE),
    pvalue = 0.05, padj = NULL, minSize = 2, minGSSize = 3
  ))
}

read_manifest <- function(p) {
  jsonlite::fromJSON(file.path(p, "manifest.json"), simplifyVector = FALSE)
}

dir_state <- function(p) {
  files <- sort(list.files(p, all.files = TRUE, no.. = TRUE))
  stats::setNames(unname(tools::md5sum(file.path(p, files))), files)
}

# ---- writeRich: folder content (AC-010, AC-012) -----------------------------

test_that("writeRich writes the manifest and identifiers verbatim", {
  # Catches identifiers converted to numbers (L-P2: "00010" read as 10) and a
  # manifest that does not say what the folder holds.
  x <- fs22()
  expect_identical(nrow(x@result), 18L)
  expect_identical(sum(!testedTerms(x)$reported), 4L)
  expect_true(all(grepl("^[0-9]{5}$", x@result$Annot)))
  p <- io_dir()
  expect_invisible(writeRich(x, p))

  expect_setequal(list.files(p), c("manifest.json", "result.tsv",
                                   "tested.tsv", "detail.tsv"))
  m <- read_manifest(p)
  expect_identical(m$format, "richR-result")
  expect_identical(m$format_version, 1L)
  expect_identical(m$class, "richResult")
  expect_true(m$tested_recorded)
  ids <- utils::read.delim(file.path(p, "result.tsv"), quote = "",
                           na.strings = character(0), comment.char = "",
                           colClasses = "character")$Annot
  expect_identical(ids, x@result$Annot)
})

test_that("writeRich records a zero-column table in the manifest only", {
  # Catches a header-only file for a zero-column table, which cannot be read
  # back (L-P5).
  x <- .fora_reference()
  x@detail <- x@detail[, 0, drop = FALSE]
  expect_gt(nrow(x@detail), 0)
  p <- io_dir()
  writeRich(x, p)
  expect_false(file.exists(file.path(p, "detail.tsv")))
  m <- read_manifest(p)
  detail <- Filter(function(t) identical(t$slot, "detail"), m$tables)[[1]]
  expect_null(detail$file)
  expect_identical(detail$ncol, 0L)
  expect_identical(detail$nrow, nrow(x@detail))
})

# ---- writeRich: refusals before any file operation (AC-005, AC-014) --------

test_that("writeRich refuses a grouped or incomplete view and creates nothing", {
  # Catches a grouped view silently stored ungrouped, and a view without its
  # identifier or label column stored as if complete (REQ-012).
  r <- .fora_reference()
  p <- io_dir()
  expect_error(writeRich(dplyr::group_by(r, Term), p), "@result.*grouped_df")
  expect_false(file.exists(p))
  expect_error(writeRich(dplyr::select(r, -Term), p), "Term")
  expect_false(file.exists(p))
  expect_error(writeRich(dplyr::select(r, -Annot), p), "Annot")
  expect_false(file.exists(p))
  g <- .fgsea_run(minSize = 3, maxSize = 9, pvalue = 0.05)
  expect_error(writeRich(dplyr::select(g, -pathway), p), "pathway")
  expect_false(file.exists(p))
})

test_that("writeRich never touches a folder it did not write", {
  # Catches user files deleted or mixed into a result folder, a silent
  # overwrite of an earlier result, and an unsupported column silently
  # coerced to text (REQ-012).
  x <- .fora_reference()
  d <- io_dir()
  dir.create(d)
  writeLines("user notes", file.path(d, "keep.txt"))
  before <- dir_state(d)
  expect_error(writeRich(x, d), "not empty")
  expect_identical(dir_state(d), before)
  expect_error(writeRich(x, d, overwrite = TRUE), "not written by writeRich")
  expect_identical(dir_state(d), before)

  w <- io_dir()
  writeRich(x, w)
  written <- dir_state(w)
  expect_error(writeRich(.fora_reference(pvalue = 1), w), "not empty")
  expect_identical(dir_state(w), written)

  y <- x
  y@result$when <- as.Date("2026-11-05")
  e <- io_dir()
  expect_error(writeRich(y, e), "when")
  expect_false(file.exists(e))
  y <- x
  y@tested$score <- 1:nrow(y@tested) + 0i
  expect_error(writeRich(y, e), "score")
  expect_false(file.exists(e))
  expect_error(writeRich(x@result, e), "richResult or GSEAResult")
  expect_false(file.exists(e))
})
