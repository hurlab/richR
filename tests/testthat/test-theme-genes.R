# themeGenes(): per-group gene evidence of a theme (SPEC-RICHR-007 M6.2b).
# Fixtures and their role checks live in helper-themes.R.

.tg_stats <- function() {
  list(q1 = data.frame(gene = sprintf("G%02d", 1:7), log2FC = seq(-1.5, 1.5, length.out = 7)),
       q2 = data.frame(gene = sprintf("G%02d", 1:7), log2FC = (1:7) / 2))
}

.tg_themes <- function() {
  x <- .th_fth2()
  .th_fth2_roles(x$q1, x$q2)
  th <- richThemes(list(q1 = x$q1, q2 = x$q2, f = NULL))
  testthat::expect_identical(th$themes$representative[1], "K")
  th$themes$label[1] <- "Edited"
  th
}

.tg_row <- function(r, gene, group) r[r$gene == gene & r$group == group, ]

test_that("a theme's genes are shown per group with the members that list them and the group's status (catches one theme asserted for groups that different genes support, and a label edit that does not reach the gene view)", {
  th <- .tg_themes()
  r <- themeGenes(th, "theme_001", .tg_stats())
  expect_identical(names(r), c("theme_id", "label", "gene", "group", "evidence_members",
                               "n_evidence_members", "group_status", "log2FC", "stat_status"))
  expect_identical(nrow(r), 24L)
  expect_setequal(r$gene, sprintf("G%02d", 1:8))
  expect_identical(levels(r$group), c("q1", "q2", "f"))
  expect_true(all(r$theme_id == "theme_001"))
  expect_true(all(r$label == "Edited"))

  for (g in sprintf("G%02d", 1:6)) {
    expect_identical(.tg_row(r, g, "q1")$evidence_members, "K,M")
    expect_identical(.tg_row(r, g, "q1")$n_evidence_members, 2L)
    expect_identical(.tg_row(r, g, "q2")$evidence_members, "")
  }
  for (g in c("G07", "G08")) {
    expect_identical(.tg_row(r, g, "q2")$evidence_members, "K")
    expect_identical(.tg_row(r, g, "q1")$evidence_members, "")
    expect_identical(.tg_row(r, g, "q1")$n_evidence_members, 0L)
  }
  expect_true(all(r$evidence_members[r$group == "f"] == ""))
  expect_true(all(r$group_status[r$group == "q1"] == "reported"))
  expect_true(all(r$group_status[r$group == "q2"] == "tested_only"))
  expect_true(all(r$group_status[r$group == "f"] == "failed"))
})

test_that("a group where one member is reported and another only tested has status reported (catches the theme's status in a group taken from a member that is only tested)", {
  G <- sprintf("G%02d", 1:80)
  x <- .th_fth2()
  q3 <- .th_enrich(G[c(1, 2, 7, 8, 41, 42)])
  b <- bindRich(list(q1 = x$q1, q3 = q3))
  expect_identical(c(.th_state(b, "K", "q3"), .th_state(b, "M", "q3")),
                   c("reported", "tested_not_reported"))
  th <- richThemes(list(q1 = x$q1, q3 = q3))
  expect_identical(.th_partition(th), c("KM", "L"))
  r <- themeGenes(th, th$themes$theme_id[th$themes$representative == "K"])
  expect_true(all(r$group_status[r$group == "q3"] == "reported"))
  expect_identical(.tg_row(r, "G01", "q3")$evidence_members, "K,M")
  expect_identical(.tg_row(r, "G07", "q3")$evidence_members, "K")
})

test_that("statistics are joined per group and a missing one is NA with its reason, never 0 (catches missing statistics shown as 0, one group's statistics joined to another, and statistics shown for a failed group)", {
  th <- .tg_themes()
  st <- .tg_stats()
  r <- themeGenes(th, "theme_001", st)
  for (g in c("q1", "q2")) {
    for (gene in sprintf("G%02d", 1:7)) {
      row <- .tg_row(r, gene, g)
      expect_identical(row$log2FC, st[[g]]$log2FC[st[[g]]$gene == gene])
      expect_identical(row$stat_status, "joined")
    }
    row8 <- .tg_row(r, "G08", g)
    expect_identical(row8$log2FC, NA_real_)
    expect_identical(row8$stat_status, "gene_not_in_statistics")
  }
  expect_identical(.tg_row(r, "G04", "q1")$log2FC, 0)
  expect_identical(.tg_row(r, "G04", "q2")$log2FC, 2)
  expect_true(all(r$stat_status[r$group == "f"] == "group_failed"))
  expect_true(all(is.na(r$log2FC[r$group == "f"])))

  r1 <- themeGenes(th, "theme_001", st["q1"])
  expect_true(all(is.na(r1$log2FC[r1$group == "q2"])))
  expect_true(all(r1$stat_status[r1$group == "q2"] == "no_statistics_for_group"))
  expect_identical(r1$log2FC[r1$group == "q1"], r$log2FC[r$group == "q1"])

  rf <- themeGenes(th, "theme_001", c(st, list(f = data.frame(gene = "G01", log2FC = 9))))
  expect_true(all(rf$stat_status[rf$group == "f"] == "group_failed"))
  expect_true(all(is.na(rf$log2FC[rf$group == "f"])))
})

test_that("a statistic column whose class or factor levels differ between stats entries is refused, and equal ones join exactly (catches a statistic lost as NA while stat_status says joined, and numbers silently turned into text)", {
  th <- .tg_themes()
  genes <- sprintf("G%02d", 1:7)
  lev <- c("down", "up")
  dir1 <- rep(lev, length.out = 7)
  frame <- function(v) data.frame(gene = genes, direction = v)
  bad <- list(
    factor_levels = list(q1 = frame(factor(dir1, levels = lev)),
                         q2 = frame(factor(rep("c", 7), levels = c("c", "d")))),
    numeric_character = list(q1 = frame((1:7) / 2), q2 = frame(as.character(1:7))),
    character_numeric = list(q1 = frame(as.character(1:7)), q2 = frame((1:7) / 2)),
    # A Date and a double share typeof "double" and have no levels; only the
    # class tells them apart.
    date_double = list(q1 = frame(as.Date("2026-01-01") + 0:6), q2 = frame((1:7) / 2)),
    # The mismatch sits in a later statistic column, not the first one.
    second_column = list(q1 = data.frame(gene = genes, lfc = (1:7) / 2, direction = dir1),
                         q2 = data.frame(gene = genes, lfc = (1:7) / 4, direction = 1:7)))
  for (k in names(bad)) {
    e <- tryCatch({
      suppressWarnings(themeGenes(th, "theme_001", bad[[k]]))
      NA_character_
    }, error = conditionMessage)
    expect_false(is.na(e), info = k)
    expect_match(e, "'direction'", fixed = TRUE, info = k)
    expect_match(e, "'q1'", fixed = TRUE, info = k)
    expect_match(e, "'q2'", fixed = TRUE, info = k)
  }

  ok <- list(factor = list(q1 = frame(factor(dir1, levels = lev)),
                           q2 = frame(factor(rev(dir1), levels = lev))),
             character = list(q1 = frame(dir1), q2 = frame(rev(dir1))))
  for (k in names(ok)) {
    st <- ok[[k]]
    r <- themeGenes(th, "theme_001", st)
    expect_identical(class(r$direction), class(st$q1$direction), info = k)
    expect_identical(levels(r$direction), levels(st$q1$direction), info = k)
    for (g in c("q1", "q2")) {
      for (gene in genes) {
        row <- .tg_row(r, gene, g)
        expect_identical(row$direction, st[[g]]$direction[st[[g]]$gene == gene], info = k)
        expect_identical(row$stat_status, "joined", info = k)
      }
      expect_true(is.na(.tg_row(r, "G08", g)$direction), info = k)
    }
  }
})

test_that("without statistics ORA has no statistic column and GSEA uses the group's ranked input (catches a statistic invented for ORA and the GSEA drill-down to input statistics broken)", {
  th <- .tg_themes()
  r <- themeGenes(th)
  expect_identical(names(r), c("theme_id", "label", "gene", "group", "evidence_members",
                               "n_evidence_members", "group_status", "stat_status"))
  expect_setequal(unique(r$theme_id), th$themes$theme_id)
  expect_true(all(r$stat_status[r$group != "f"] == "no_statistics_for_group"))
  expect_true(all(r$stat_status[r$group == "f"] == "group_failed"))

  skip_if_not_installed("fgsea")
  t1 <- .th_fgsea(TRUE)
  tg <- themeGenes(richThemes(list(a = t1$a, b = t1$b)))
  expect_gt(nrow(tg), 0L)
  expect_true("rank_score" %in% names(tg))
  for (g in c("a", "b")) {
    rows <- tg[tg$group == g, ]
    expect_identical(rows$rank_score, unname(t1[[g]]@input[rows$gene]))
    expect_true(all(rows$stat_status == "joined"))
  }
})

test_that("invalid themes, theme ids and statistics stop with an error naming the offending element (catches an ambiguous join taking one of two rows, statistics for an unknown group ignored, and an unknown theme read as no evidence)", {
  th <- .tg_themes()
  fr <- function(...) data.frame(gene = c("G01", "G02"), ...)
  expect_error(themeGenes(list()), "themes", fixed = TRUE)
  expect_error(themeGenes(th, "theme_999"), "theme_999", fixed = TRUE)
  expect_error(themeGenes(th, stats = fr(v = 1:2)), "stats", fixed = TRUE)
  expect_error(themeGenes(th, stats = list(fr(v = 1:2))), "stats", fixed = TRUE)
  expect_error(themeGenes(th, stats = list(q3 = fr(v = 1:2))), "q3", fixed = TRUE)

  e <- tryCatch(themeGenes(th, stats = list(q1 = data.frame(g = "G01", v = 1))),
                error = conditionMessage)
  expect_match(e, "q1", fixed = TRUE); expect_match(e, "gene", fixed = TRUE)
  e <- tryCatch(themeGenes(th, stats = list(q1 = data.frame(gene = 1:2, v = 1:2))),
                error = conditionMessage)
  expect_match(e, "q1", fixed = TRUE); expect_match(e, "gene", fixed = TRUE)
  e <- tryCatch(themeGenes(th, stats = list(q2 = data.frame(gene = c("G01", "G01"), v = 1:2))),
                error = conditionMessage)
  expect_match(e, "q2", fixed = TRUE); expect_match(e, "G01", fixed = TRUE)
  expect_error(themeGenes(th, stats = list(q1 = data.frame(gene = "G01"))), "q1", fixed = TRUE)
  expect_error(themeGenes(th, stats = list(q1 = fr(group = c("a", "b")))), "group", fixed = TRUE)
  e <- tryCatch(themeGenes(th, stats = list(q1 = fr(log2FC = 1:2), q2 = fr(score = 1:2))),
                error = conditionMessage)
  expect_match(e, "log2FC", fixed = TRUE); expect_match(e, "score", fixed = TRUE)
})

test_that("an empty theme selection is refused, while an analysis without themes still gives an empty table (catches an empty request returning a table that reads as no evidence)", {
  th <- .tg_themes()
  expect_error(themeGenes(th, character(0)), "theme must name at least one theme_id",
               fixed = TRUE)
  tz <- themeGenes(richThemes(.th_fzero()))
  expect_identical(nrow(tz), 0L)
  expect_identical(names(tz), c("theme_id", "label", "gene", "group", "evidence_members",
                                "n_evidence_members", "group_status", "stat_status"))
})

# SPEC-RICHR-009 AC-005: an omitted zero-overlap term is not known absent.
test_that("zero overlap or untested members are undetermined, below positive evidence", {
  annot <- data.frame(GeneID = paste0("G", 1:60),
                       Annot = rep(paste0("T", 1:6), each = 10))
  annot$Term <- annot$Annot
  run <- function(q, ...) enrich(paste0("G", q), annot, pvalue = 1, padj = 1,
                                 minSize = 1, minGSSize = 1, ...)
  q1 <- run(c(1:8, 11:12))
  qa <- run(c(11:18, 21:24), padj.family = "all")
  qo <- run(c(11:18, 21:24))
  groups <- list(q1 = q1, qa = qa, qo = qo, qu = as(qo@result, "richResult"))
  b <- bindRich(groups)
  expect_identical(b$state[b$Annot == "T1"],
    c("reported", "untested_or_zero_overlap", "not_tested", "unavailable"))
  th <- richThemes(groups)
  tid <- th$members$theme_id[th$members$Annot == "T1"][1]
  r <- themeGenes(th, tid)
  expect_identical(unname(vapply(split(r$group_status, r$group), unique, "")),
                   c("reported", "undetermined", "absent", "unavailable"))

  annot <- rbind(data.frame(GeneID = paste0("G", 1:10), Annot = "K"),
    data.frame(GeneID = paste0("G", c(1:8, 41:44)), Annot = "M"),
    data.frame(GeneID = paste0("G", 21:30), Annot = "Z1"),
    data.frame(GeneID = paste0("G", 31:40), Annot = "Z2"),
    data.frame(GeneID = paste0("G", 45:60), Annot = "Z3"))
  annot$Term <- annot$Annot
  run <- function(q, ...) enrich(paste0("G", q), annot,
                                 minSize = 1, minGSSize = 1, ...)
  groups <- list(s1 = run(c(1:9, 41:43), pvalue = 1, padj = 1),
    sr = run(c(41:44, 21:23, 31:33), pvalue = 1, padj = 1, padj.family = "all"),
    st = run(c(41, 21:30, 31:35), pvalue = 0.001, padj = 1, padj.family = "all"))
  b <- bindRich(groups)
  expect_identical(b$state[b$Annot == "K"],
                   c("reported", "untested_or_zero_overlap", "untested_or_zero_overlap"))
  expect_identical(b$state[b$Annot == "M"],
                   c("reported", "reported", "tested_not_reported"))
  th <- richThemes(groups)
  expect_identical(th$members$theme_id[th$members$Annot == "K"][1],
                   th$members$theme_id[th$members$Annot == "M"][1])
  r <- themeGenes(th, th$members$theme_id[th$members$Annot == "K"][1])
  expect_identical(unname(vapply(split(r$group_status, r$group), unique, "")),
                   c("reported", "reported", "tested_only"))
})
