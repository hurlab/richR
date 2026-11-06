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
