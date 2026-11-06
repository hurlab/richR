# Regenerates plot_reference_m6base.rds (SPEC-RICHR-006 AC-009, REQ-006).
#
# The reference holds the plot data, layer data and titles of every plot of
# set P (richBar, richDot, richNetplot, richNetwork, richNetmap, richLollipop,
# richCircle, richNES, richScatter, richECDF, richTermSim, richGeneDot,
# richGeneHeat, richGeneBar, richVolcano, richHeatmap, richCompareDot; every
# method) plus richClusterDot and richUpset, on inputs whose rows all pass the
# default selection. It was recorded from an install of the M6.1 base
# (fca5276) BEFORE any M6.1 change and must never be regenerated to absorb a
# later change: test-plot-reference.R compares the current build against it.
#
# The calls and inputs live in ../helper-plot-fixtures.R (shared with the test).
#
# Usage, from any directory, with the richR under test first on the library path:
#   R_LIBS=<lib with richR>:<provider lib> Rscript --vanilla make-plot-reference.R
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
suppressPackageStartupMessages(library(richR))
source(file.path(here, "..", "helper-plot-fixtures.R"))
cat("richR", as.character(utils::packageVersion("richR")), "from", find.package("richR"), "\n")
fgsea_version <- if (requireNamespace("fgsea", quietly = TRUE)) {
  as.character(utils::packageVersion("fgsea"))
} else NA_character_
calls <- .m61_reference_calls(here, gsea = !is.na(fgsea_version))
ref <- list(
  plots = .m61_record(calls),
  meta = list(fgsea = fgsea_version,
              ggplot2 = as.character(utils::packageVersion("ggplot2")),
              R = paste(R.version$major, R.version$minor, sep = "."))
)
out <- file.path(here, "plot_reference_m6base.rds")
saveRDS(ref, out, version = 3)
cat("recorded", length(ref$plots), "plot calls into", out, "\n")
