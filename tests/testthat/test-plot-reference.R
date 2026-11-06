# AC-009 (SPEC-RICHR-006 REQ-006, REQ-008): every plot of set P (each method),
# richClusterDot and richUpset draw on the current build exactly what the M6.1
# base (fca5276) drew, on inputs whose rows all pass the default selection.
# The reference fixtures/plot_reference_m6base.rds was recorded by
# fixtures/make-plot-reference.R from an install of the base before any M6.1
# change; it is never regenerated to absorb a change.
#
# Changed by design, and therefore compared only outside the change:
#   - IC-2 / AC-005: the GSEA richLollipop, richGeneBar and richScatter draw NES
#     instead of a constant; those three calls are left out here (AC-005 checks
#     what they draw).
#   - IC-5 / AC-008: richHeatmap draws an absent cell as NA instead of p = 1 and
#     names its fill after the column drawn. Everything else (row order, labels,
#     every recorded p-value, every other fill) must equal the reference.
# Nothing else may differ.

.m61_ref <- function() readRDS(testthat::test_path("fixtures", "plot_reference_m6base.rds"))

# By-design exclusions (IC-2 / AC-005).
.m61_ic2_calls <- c("gsea_richLollipop", "gsea_richGeneBar", "gsea_richScatter")

# Compare one heat-map call outside the IC-5 change: absent cells (NA now,
# 1 on the base) and the drawn column's name / fill title.
.m61_expect_heatmap_equal <- function(now, ref, name) {
  if (!is.null(ref$table)) {
    tab <- now$table
    absent <- is.na(as.matrix(tab[, -1, drop = FALSE]))
    expect_true(all(as.matrix(ref$table[, -1, drop = FALSE])[absent] == 1), info = name)
    tab[, -1][absent] <- 1
    expect_identical(tab, ref$table, info = name)
    return(invisible())
  }
  d <- now$data
  absent <- is.na(d[[3]])
  expect_true(all(ref$data[[3]][absent] == 1), info = name)
  d[[3]][absent] <- 1
  names(d)[3] <- names(ref$data)[3]
  expect_identical(d, ref$data, info = name)
  expect_equal(length(now$layers), length(ref$layers), info = name)
  for (i in seq_along(ref$layers)) {
    ln <- now$layers[[i]]; lr <- ref$layers[[i]]
    keep <- !absent
    expect_identical(ln$fill[keep], lr$fill[keep], info = paste(name, "layer", i, "fill"))
    ln$fill <- NULL; lr$fill <- NULL
    expect_identical(ln, lr, info = paste(name, "layer", i))
  }
  labs_now <- now$labs; labs_ref <- ref$labs
  labs_now$fill <- NULL; labs_ref$fill <- NULL
  expect_identical(labs_now, labs_ref, info = name)
}

.m61_expect_reference <- function(now, ref) {
  expect_gt(length(now), 0)
  for (name in names(now)) {
    expect_true(name %in% names(ref), info = name)
    if (name %in% .m61_ic2_calls) next
    if (startsWith(name, "richHeatmap")) {
      .m61_expect_heatmap_equal(now[[name]], ref[[name]], name)
    } else {
      expect_identical(now[[name]], ref[[name]], info = name)
    }
  }
}

test_that("ORA plots draw what the M6.1 base drew (AC-009)", {
  # Catches: any unintended change in what an existing plot call draws.
  ref <- .m61_ref()
  calls <- .m61_reference_calls(testthat::test_path("fixtures"), gsea = FALSE)
  now <- .m61_record(calls)
  expect_identical(sort(names(now)), sort(grep("^gsea_", names(ref$plots), value = TRUE, invert = TRUE)))
  .m61_expect_reference(now, ref$plots)
})

test_that("GSEA plots draw what the M6.1 base drew (AC-009)", {
  # Catches: any unintended change in what an existing GSEA plot call draws.
  .char_skip_unless_recorded_fgsea()
  ref <- .m61_ref()
  expect_identical(ref$meta$fgsea, as.character(utils::packageVersion("fgsea")))
  calls <- .m61_reference_calls(testthat::test_path("fixtures"), gsea = TRUE)
  calls <- calls[startsWith(names(calls), "gsea_")]
  now <- .m61_record(calls)
  expect_identical(sort(names(now)), sort(grep("^gsea_", names(ref$plots), value = TRUE)))
  .m61_expect_reference(now, ref$plots)
})
