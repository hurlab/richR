# Regenerates the small deterministic fixtures used by the characterization
# tests (test-characterization-*.R). No random numbers, no network.
#
# Run from tests/testthat/fixtures/:  Rscript make-fixtures.R
# The committed CSV files are the source of truth for the tests; this script
# documents how they were derived so a reviewer can reproduce them.

write_fixture <- function(x, file) {
  utils::write.csv(x, file, row.names = FALSE, quote = TRUE)
}

## ---- ORA: 30 annotated genes, 6 terms, 11 query genes -----------------------
genes <- sprintf("G%02d", 1:30)
mk <- function(term, name, idx) {
  data.frame(GeneID = genes[idx], Term = term, Annot = name, stringsAsFactors = FALSE)
}
ora_annot <- rbind(
  mk("GO:0000001", "alpha process", 1:8),
  mk("GO:0000002", "beta process", 5:12),
  mk("GO:0000003", "gamma process", c(1, 2, 13:15)),
  mk("GO:0000004", "delta process", 16:24),
  mk("GO:0000005", "epsilon process", c(3, 4, 9, 25:27)),
  mk("GO:0000006", "zeta process", 1:3)
)
write_fixture(ora_annot, "ora_annot.csv")

ora_query <- genes[c(1:6, 9, 10, 13, 14, 28)]
write_fixture(data.frame(GeneID = ora_query, stringsAsFactors = FALSE), "ora_query.csv")

# Differential-expression style input (data.frame x, genes as the GeneID column).
ora_degs <- data.frame(
  GeneID = ora_query,
  lfc = c(2.5, 1.75, 1.5, 1.25, 1.0, 0.75, -0.5, -1.0, -1.5, -2.0, 0.25),
  stringsAsFactors = FALSE
)
write_fixture(ora_degs, "ora_degs.csv")

# KEGG-style membership: terms are bundled KEGG map ids with distinct Level2
# groups (so richLevel has several groups to test),
# so richKEGG/richLevel resolve names and levels from data/kegg.rdata and
# data/path.rda without any network access.
kegg_idx <- list(
  "00010" = 1:8,  # Carbohydrate metabolism
  "00061" = 5:12,  # Lipid metabolism
  "00190" = c(1, 2, 13:15),  # Energy metabolism
  "00230" = 16:24,  # Nucleotide metabolism
  "00970" = c(3, 4, 9, 25:27)  # Translation (Level1 differs)
)
kegg_annot <- do.call(rbind, lapply(names(kegg_idx), function(k) {
  data.frame(GeneID = genes[kegg_idx[[k]]], Term = k, stringsAsFactors = FALSE)
}))
write_fixture(kegg_annot, "kegg_annot.csv")

## ---- GSEA: 40 ranked genes, 5 gene sets ------------------------------------
ggenes <- sprintf("G%02d", 1:40)
mkg <- function(term, name, idx) {
  data.frame(GeneID = ggenes[idx], Term = term, Annot = name, stringsAsFactors = FALSE)
}
gsea_annot <- rbind(
  mkg("T1", "up pathway", c(1:6, 8, 10)),
  mkg("T2", "down pathway", 31:38),
  mkg("T3", "mixed pathway", c(2, 9, 15, 22, 28, 33, 39)),
  mkg("T4", "top pathway", 1:5),
  mkg("T5", "middle pathway", 16:25)
)
write_fixture(gsea_annot, "gsea_annot.csv")

# Strictly decreasing-with-jitter scores (no ties), rounded to 6 decimals.
score_a <- rev(seq(-2, 2, length.out = 40)) + 0.01 * ((1:40) %% 7)
score_b <- score_a + 0.02 * ((1:40) %% 5 - 2)
write_fixture(data.frame(GeneID = ggenes, score = round(score_a, 6)), "gsea_ranked_a.csv")
write_fixture(data.frame(GeneID = ggenes, score = round(score_b, 6)), "gsea_ranked_b.csv")

## ---- richCluster: five terms, two overlapping neighbourhoods ----------------
cluster_terms <- data.frame(
  Term = LETTERS[1:5], Annot = letters[1:5],
  GeneID = c(
    paste0("g", 1:5, collapse = ","),
    paste0("g", 2:6, collapse = ","),
    paste0("g", 11:15, collapse = ","),
    paste0("g", 12:16, collapse = ","),
    paste0("g", 17:40, collapse = ",")
  ),
  Pvalue = c(1e-8, 1e-7, 1e-6, 1e-5, 0.5),
  Padj = c(1e-8, 1e-7, 1e-6, 1e-5, 0.5),
  stringsAsFactors = FALSE
)
write_fixture(cluster_terms, "cluster_terms.csv")
