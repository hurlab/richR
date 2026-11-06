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

# Round-trip equality of acceptance.md § D.0: per table identical column
# names, classes, row names and columns (doubles bit-exact); a zero-column
# table identical as a whole; every other slot identical; the tested table
# unavailable exactly when it was; enrichmentInfo() identical.
expect_same_table <- function(b, a, what) {
  expect_identical(names(b), names(a), info = what)
  expect_identical(lapply(b, class), lapply(a, class), info = what)
  expect_identical(rownames(b), rownames(a), info = what)
  for (j in seq_along(a)) {
    expect_identical(b[[j]], a[[j]], info = paste(what, names(a)[j]))
  }
  if (!ncol(a)) expect_identical(b, a, info = what)
}

tested_available <- function(x) !is.null(suppressWarnings(testedTerms(x)))

expect_round_trip <- function(b, a) {
  expect_identical(class(b), class(a))
  tables <- if (methods::is(a, "richResult")) c("result", "detail") else "result"
  for (s in tables) expect_same_table(methods::slot(b, s), methods::slot(a, s), s)
  expect_identical(tested_available(b), tested_available(a))
  if (tested_available(a)) expect_same_table(b@tested, a@tested, "tested")
  for (s in setdiff(methods::slotNames(class(a)), c(tables, "tested"))) {
    expect_identical(methods::slot(b, s), methods::slot(a, s), info = s)
  }
  expect_identical(enrichmentInfo(b), enrichmentInfo(a))
}

round_trip <- function(x) {
  p <- io_dir()
  writeRich(x, p)
  readRich(p)
}

copy_folder <- function(src) {
  d <- io_dir()
  dir.create(d)
  file.copy(list.files(src, full.names = TRUE), d)
  d
}

write_manifest <- function(p, m) {
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, null = "null",
                              digits = I(17), pretty = TRUE),
             file.path(p, "manifest.json"))
}

edit_manifest <- function(p, f) write_manifest(p, f(read_manifest(p)))

# Edits the typed metadata the manifest embeds (jsonlite::serializeJSON form).
edit_metadata <- function(p, f) {
  edit_manifest(p, function(m) {
    meta <- jsonlite::unserializeJSON(jsonlite::toJSON(
      m$metadata, auto_unbox = TRUE, null = "null", digits = I(17)))
    m$metadata <- jsonlite::fromJSON(
      jsonlite::serializeJSON(f(meta), digits = I(17)), simplifyVector = FALSE)
    m
  })
}

edit_lines <- function(p, file, f) {
  path <- file.path(p, file)
  writeLines(f(readLines(path, encoding = "UTF-8")), path, useBytes = TRUE)
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

# ---- readRich: round trip (AC-005, AC-010..AC-012, AC-014) ------------------

test_that("readRich restores an ORA result exactly, statistics and labels included", {
  # Catches precision loss in stored statistics (L-P1), identifiers converted
  # to numbers (L-P2), NA / "" / "NA" conflated or rows lost (L-P4b), a TAB or
  # newline splitting a cell, factor levels lost, extra columns of a pipeline
  # result dropped, and metadata dropped or retyped on reload (L-P3).
  x <- fs22_edited()
  expect_no_error(methods::validObject(x))
  b <- round_trip(x)
  expect_round_trip(b, x)
  expect_identical(b@result$zscore[1:6], c(NA, NaN, Inf, -Inf, 5e-324, 1 / 3))
  expect_identical(b@result$Term[1:17], lp4_values)
  expect_identical(levels(b@detail$flag), c("b", "a", "z"))

  k <- kegg_result()
  expect_true(all(c("ko", "Level3", "Level2", "Level1") %in% names(k@result)))
  expect_round_trip(round_trip(k), k)
})

test_that("readRich restores GSEA results with their ranking, leading edges and provenance", {
  # Catches a GSEA ranking reordered or renamed on reload, a list column
  # flattened or split wrongly, an empty leading edge reloaded as list(), and
  # GSEA provenance dropped.
  skip_if_not_installed("fgsea")
  ranked <- .fgsea_ranked()
  ranked[["G05"]] <- ranked[["G04"]]
  run <- function(annot, ...) {
    out <- NULL
    utils::capture.output(out <- suppressWarnings(richR::richGSEA(
      ranked, annot, seed = 42, nproc = 1, minSize = 3, maxSize = 9,
      pvalue = 0.05, ...)))
    out
  }
  annot <- .char_read("gsea_annot.csv", character = TRUE)
  gt <- run(annot, table = TRUE)
  expect_gt(anyDuplicated(unname(gt@input)), 0)
  expect_true(any(gt@input < 0))
  gl <- run(annot, table = FALSE)
  expect_gt(nrow(gl@result), 0)
  expect_type(gl@result$leadingEdge, "list")
  gl@result$leadingEdge[[1]] <- character(0)
  ak <- annot
  ak$Term <- unname(c(T1 = "00010", T2 = "00061", T3 = "00190", T4 = "00230",
                      T5 = "00970")[ak$Term])
  gk <- run(ak, KEGG = TRUE)
  expect_true("Level1" %in% names(testedTerms(gk)))

  for (g in list(gt, gl, gk)) expect_round_trip(round_trip(g), g)
  bl <- round_trip(gl)
  expect_identical(bl@result$leadingEdge[[1]], character(0))
  expect_identical(names(bl@input), names(gl@input))
  info <- enrichmentInfo(bl)$provenance
  expect_false(is.null(info$fgsea_version))
  expect_false(is.null(info$rng))
})

test_that("readRich restores empty results, zero-column tables and unavailable records", {
  # Catches an empty result rejected or reloaded with logical columns, a
  # zero-column table that cannot be read back (L-P5) or loses its row count,
  # and a legacy object reloaded with a fabricated tested table or upgraded
  # metadata.
  none <- .fora_reference(pvalue = 1e-12)
  expect_identical(nrow(none@result), 0L)
  expect_gt(nrow(testedTerms(none)), 0)
  expect_round_trip(round_trip(none), none)

  emptied <- .fora_reference()
  emptied@tested <- emptied@tested[0, ]
  expect_round_trip(round_trip(emptied), emptied)

  zero <- .fora_reference()
  zero@detail <- zero@detail[, 0, drop = FALSE]
  expect_round_trip(round_trip(zero), zero)
  named <- .fora_reference()
  named@detail <- data.frame(a = 1:3, row.names = c("x", "y", "z"))[, 0, drop = FALSE]
  expect_round_trip(round_trip(named), named)

  as_obj <- methods::as(.fora_reference()@result, "richResult")
  expect_identical(dim(as_obj@detail), c(0L, 0L))
  p <- io_dir()
  writeRich(as_obj, p)
  expect_false(file.exists(file.path(p, "detail.tsv")))
  b <- readRich(p)
  expect_round_trip(b, as_obj)
  expect_false(enrichmentInfo(b)$available)

  for (name in c("legacy_richResult_fca5276.rds", "legacy_GSEAResult_fca5276.rds")) {
    legacy <- readRDS(test_path("fixtures", name))
    p <- io_dir()
    writeRich(legacy, p)
    expect_false(file.exists(file.path(p, "tested.tsv")))
    b <- readRich(p)
    expect_round_trip(b, legacy)
    expect_warning(tt <- testedTerms(b), "tested table unavailable")
    expect_null(tt)
    expect_identical(enrichmentInfo(b)$schema_version, 1L)
  }
})

test_that("readRich restores a GSEA result that reported nothing", {
  # Catches an empty GSEA result rejected or reloaded with retyped columns.
  skip_if_not_installed("fgsea")
  g <- .fgsea_run(minSize = 3, maxSize = 9, pvalue = 1e-12, padj = NULL)
  expect_identical(nrow(g@result), 0L)
  expect_gt(nrow(testedTerms(g)), 0)
  expect_round_trip(round_trip(g), g)
})

test_that("a filtered view round-trips with the original analysis record", {
  # Catches a derived view stored with a pruned tested table or rewritten
  # metadata (AC-005 I/O leg).
  r <- .fora_reference()
  f <- dplyr::filter(r, Annot != "T1")
  expect_lt(nrow(f@result), nrow(r@result))
  b <- round_trip(f)
  expect_round_trip(b, f)
  expect_identical(testedTerms(b), testedTerms(r))
  expect_identical(enrichmentInfo(b), enrichmentInfo(r))
})

test_that("overwrite = TRUE replaces a previous writeRich folder", {
  # Catches an overwrite that leaves the earlier result's files behind (a stale
  # tested.tsv would make the new folder unreadable).
  w <- io_dir()
  writeRich(.fora_reference(), w)
  x <- methods::as(.fora_reference(pvalue = 1)@result, "richResult")
  expect_no_error(writeRich(x, w, overwrite = TRUE))
  expect_false(file.exists(file.path(w, "tested.tsv")))
  expect_round_trip(readRich(w), x)
})

# ---- readRich: refusals (AC-013) ---------------------------------------------

test_that("readRich refuses a folder it cannot read exactly", {
  # Catches a newer, older or foreign folder read as if it were understood, a
  # truncated or hand-edited table returned as a valid result, and a corrupt
  # number silently read as NA.
  x <- .fora_reference(pvalue = 1)
  expect_identical(nrow(x@result), 2L)
  base <- io_dir()
  writeRich(x, base)
  # The edits below, not the rewriting itself, are what readRich refuses.
  same <- copy_folder(base)
  edit_manifest(same, identity)
  expect_round_trip(readRich(same), x)
  changed <- function(f) {
    p <- copy_folder(base)
    f(p)
    p
  }
  set_field <- function(field, value) {
    changed(function(p) edit_manifest(p, function(m) {
      m[field] <- list(value)
      if (is.null(value)) m[[field]] <- NULL
      m
    }))
  }
  set_schema <- function(value) {
    changed(function(p) edit_metadata(p, function(meta) {
      meta["schema_version"] <- list(value)
      if (is.null(value)) meta$schema_version <- NULL
      meta
    }))
  }

  expect_error(readRich(set_field("format_version", 2L)),
               "format version 2 .*newer.*supported: 1\\)")
  expect_error(readRich(set_field("format_version", 0L)),
               "format version 0 .*older.*supported: 1\\)")
  expect_error(readRich(set_field("format_version", "1")),
               "format version \"1\" .*not a whole number.*supported: 1\\)")
  expect_error(readRich(set_field("format_version", 1.5)),
               "format version 1.5 .*not a whole number.*supported: 1\\)")
  expect_error(readRich(set_field("format_version", NULL)),
               "no format version .*supported: 1\\)")
  expect_error(readRich(set_schema(3L)),
               "schema_version 3 .*newer.*supported: 1, 2\\)")
  expect_error(readRich(set_schema("2")),
               "schema_version \"2\" .*not a whole number.*supported: 1, 2\\)")
  expect_error(readRich(set_schema(NULL)),
               "no metadata schema_version .*supported: 1, 2\\)")
  expect_error(readRich(set_field("format", "other")),
               "format \"other\".*richR-result")
  expect_error(readRich(set_field("class", "Annot")),
               "class \"Annot\".*richResult, GSEAResult")

  expect_error(readRich(changed(function(p) unlink(file.path(p, "tested.tsv")))),
               "tested.tsv.*missing")
  untested <- io_dir()
  writeRich(methods::as(x@result, "richResult"), untested)
  extra <- copy_folder(untested)
  file.copy(file.path(base, "tested.tsv"), extra)
  expect_error(readRich(extra), "tested.tsv.*not listed")
  extra <- copy_folder(untested)
  file.copy(file.path(base, "detail.tsv"), extra)
  expect_error(readRich(extra), "detail.tsv.*not listed")

  expect_error(readRich(changed(function(p) {
    edit_lines(p, "result.tsv", function(l) l[-length(l)])
  })), "result.tsv.* 1 data rows.*records 2")
  expect_error(readRich(changed(function(p) {
    edit_lines(p, "result.tsv", function(l) {
      l[1] <- sub("Pvalue", "Pval", l[1], fixed = TRUE)
      l
    })
  })), "header of result.tsv.*'Pvalue'")
  set_cell <- function(row, column, value) {
    changed(function(p) edit_lines(p, "result.tsv", function(l) {
      cells <- strsplit(l[row + 1L], "\t", fixed = TRUE)[[1]]
      cells[match(column, names(x@result))] <- value
      l[row + 1L] <- paste(cells, collapse = "\t")
      l
    }))
  }
  expect_error(readRich(set_cell(2L, "Pvalue", "0.05x")),
               "result.tsv.*column 'Pvalue', row 2")
  expect_error(readRich(set_cell(1L, "Term", "a\\qb")),
               "result.tsv.*column 'Term', row 1.*unknown escape")
  expect_error(readRich(changed(function(p) unlink(file.path(p, "manifest.json")))),
               "manifest.json")
})
