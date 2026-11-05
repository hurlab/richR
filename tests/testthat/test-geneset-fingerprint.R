# Gene-set fingerprints (REQ-015 .. REQ-018), canonical form richR-geneset-v1.
# Catches: a fingerprint that is missing, R-version-dependent, computed over the
# wrong columns, unstable under row order, blind to membership changes, or copied
# from a stale Annot into an analysis result. Offline; providers are mocked.

sha <- function(text) digest::digest(text, algo = "sha256", serialize = FALSE)

membership_kav <- "60d6ebd80be3bfd0698c6b4f7855aa4a7b948e16968211ebc2a00a5a15f8192f"
full_labels_kav <- "c7fa39ee64ad724d49b4b9a3a4f46d9cdc1fd2368aa13063d1fe29c0a62734cc"
full_setname_kav <- "dbea0a4fe1ec8024d693bc8cbbe326ed3405f9d289e2f51755b32e02c0596126"

fingerprint_of <- function(object) {
  richR::enrichmentInfo(object)$provenance$annotation$fingerprint
}

quiet <- function(expr) {
  utils::capture.output(value <- suppressMessages(expr))
  value
}

test_that("every builder records both known-answer fingerprints and only stated database versions", {
  msig_rows <- function(collection) {
    if (identical(collection, "H")) {
      data.frame(gs_cat = "H", gs_subcat = NA_character_,
                 gene_symbol = c("G1", "G2", "G2"),
                 gs_exact_source = c("x1", "x1", "x2"),
                 gs_name = c("HALLMARK_T1", "HALLMARK_T1", "HALLMARK_T2"),
                 db_version = "2025.1")
    } else {
      data.frame(gs_cat = "C5", gs_subcat = "GO:BP",
                 gene_symbol = c("G1", "G2", "G2"),
                 gs_exact_source = c("T1", "T1", "T2"),
                 gs_name = c("GOBP_alpha", "GOBP_alpha", "GOBP_beta"),
                 db_version = "2025.1")
    }
  }
  build_msig <- function(anntype, collection) {
    quiet(testthat::with_mocked_bindings(
      .getmsig = function(...) "Homo sapiens",
      msigdbr = function(...) msig_rows(collection),
      richR::buildMSIGDB(species = "human", keytype = "SYMBOL", anntype = anntype),
      .package = "richR"
    ))
  }
  gmt_file <- tempfile(fileext = ".gmt")
  writeLines(c("T1\tdescription\tG1\tG2", "T2\tdescription\tG2"), gmt_file)

  set_name <- list(
    list = quiet(richR::buildAnnotFromList(list(T1 = c("G1", "G2"), T2 = "G2"))),
    gmt = quiet(richR::readGMT(gmt_file)),
    msigdb_hallmark = build_msig("HALLMARK", "H")
  )
  separate <- list(
    buildAnnot_go = quiet(testthat::with_mocked_bindings(
      .makeGOdata = function(...) data.frame(
        GeneID = c("G1", "G2", "G2"), GOALL = c("T1", "T1", "T2"),
        ONTOLOGYALL = "BP", Annot = c("alpha", "alpha", "beta")),
      .getdbname = function(...) "stats",
      richR::buildAnnot(species = "human", keytype = "SYMBOL", anntype = "GO"),
      .package = "richR"
    )),
    msigdb_bp = build_msig("BP", "C5")
  )

  for (name in names(set_name)) {
    fp <- fingerprint_of(set_name[[name]])
    expect_equal(fp$membership, membership_kav, info = name)
    expect_equal(fp$full, full_setname_kav, info = name)
  }
  for (name in names(separate)) {
    fp <- fingerprint_of(separate[[name]])
    expect_equal(fp$membership, membership_kav, info = name)
    expect_equal(fp$full, full_labels_kav, info = name)
  }
  for (object in c(set_name, separate)) {
    fp <- fingerprint_of(object)
    expect_equal(fp$algorithm, "sha256")
    expect_equal(fp$canonicalization, "richR-geneset-v1")
    expect_equal(c(fp$n_pairs, fp$n_sets, fp$n_genes), c(3L, 2L, 2L))
  }

  # REQ-018: a database version only where the source states one; the provider
  # package version is a separate field.
  for (name in c("msigdb_hallmark", "msigdb_bp")) {
    annotation <- richR::enrichmentInfo(c(set_name, separate)[[name]])$provenance$annotation
    expect_equal(annotation$source_version, "2025.1", info = name)
    expect_equal(annotation$provider_version,
                 as.character(utils::packageVersion("msigdbr")), info = name)
  }
  for (name in setdiff(names(c(set_name, separate)), c("msigdb_hallmark", "msigdb_bp"))) {
    annotation <- richR::enrichmentInfo(c(set_name, separate)[[name]])$provenance$annotation
    expect_true(is.na(annotation$source_version), info = name)
  }
  reactome <- quiet(testthat::with_mocked_bindings(
    .reactome_map = function(name) list(
      reactomePATHID2EXTID = list("R-HSA-1" = c("G1", "G2"), "R-HSA-2" = "G2"),
      reactomePATHID2NAME = list("R-HSA-1" = "Homo sapiens: alpha",
                                 "R-HSA-2" = "Homo sapiens: beta")
    )[[name]],
    richR::buildAnnot(species = "human", keytype = "ENTREZID", anntype = "Reactome"),
    .package = "richR"
  ))
  reactome_annotation <- richR::enrichmentInfo(reactome)$provenance$annotation
  expect_true(is.na(reactome_annotation$source_version))
  expect_equal(reactome_annotation$fingerprint$membership,
               sha("richR-geneset-v1\nR-HSA-1\tG1\nR-HSA-1\tG2\nR-HSA-2\tG2\n"))
  expect_equal(reactome_annotation$fingerprint$full,
               sha("richR-geneset-v1/full\nR-HSA-1\tG1\talpha\nR-HSA-1\tG2\talpha\nR-HSA-2\tG2\tbeta\n"))
})

test_that("buildOwn records both known-answer fingerprints", {
  skip_if_not_installed("org.Hs.eg.db")
  own <- quiet(testthat::with_mocked_bindings(
    keys = function(x, keytype) c("G1", "G2"),
    getann = function(...) data.frame(annotation = c("alpha", "beta"),
                                      row.names = c("T1", "T2")),
    .package = "richR",
    testthat::with_mocked_bindings(
      select = function(x, keys, keytype, columns) {
        data.frame(SYMBOL = c("G1", "G2", "G2"), GOALL = c("T1", "T1", "T2"),
                   ONTOLOGYALL = "BP")
      },
      richR::buildOwn(dbname = "org.Hs.eg.db", anntype = "GO", keytype = "SYMBOL"),
      .package = "AnnotationDbi"
    )
  ))
  fp <- fingerprint_of(own)
  expect_equal(c(fp$membership, fp$full), c(membership_kav, full_labels_kav))
  expect_equal(c(fp$n_pairs, fp$n_sets, fp$n_genes), c(3L, 2L, 2L))
  expect_true(is.na(richR::enrichmentInfo(own)$provenance$annotation$source_version))
})

test_that("fingerprints ignore layout noise and follow every membership and label change", {
  fingerprint <- getFromNamespace(".geneset_fingerprint", "richR")
  base <- data.frame(GeneID = c("G1", "G2", "G2"), Term = c("T1", "T1", "T2"),
                     Annot = c("alpha", "alpha", "beta"), stringsAsFactors = FALSE)
  ref <- fingerprint(base)
  expect_equal(c(ref$membership, ref$full), c(membership_kav, full_labels_kav))

  same <- list(
    permuted = base[c(3, 1, 2), ],
    duplicated = base[c(1, 2, 2, 3), ],
    row_names = `rownames<-`(base, c("r1", "r2", "r3")),
    extra_column = cbind(base, Extra = c("x", "y", "z")),
    factors = data.frame(GeneID = factor(base$GeneID), Term = factor(base$Term),
                         Annot = factor(base$Annot))
  )
  for (name in names(same)) {
    fp <- fingerprint(same[[name]])
    expect_equal(fp$membership, ref$membership, info = name)
    expect_equal(fp$full, ref$full, info = name)
  }

  relabelled <- fingerprint(transform(base, Annot = c("alpha", "alpha", "gamma")))
  expect_equal(relabelled$membership, ref$membership)
  expect_false(identical(relabelled$full, ref$full))

  changed <- list(
    added = fingerprint(rbind(base, data.frame(GeneID = "G3", Term = "T2", Annot = "beta"))),
    removed = fingerprint(base[-1, ]),
    changed = fingerprint(transform(base, GeneID = c("G1", "G3", "G2")))
  )
  for (name in names(changed)) {
    expect_false(identical(changed[[name]]$membership, ref$membership), info = name)
    expect_false(identical(changed[[name]]$full, ref$full), info = name)
  }
  memberships <- vapply(changed, `[[`, character(1), "membership")
  expect_equal(length(unique(memberships)), 3)

  tabbed <- fingerprint(data.frame(GeneID = "G1", Term = "T\t1", Annot = "lab"))
  expect_equal(tabbed$membership, sha("richR-geneset-v1\nT\\t1\tG1\n"))
  expect_equal(tabbed$full, sha("richR-geneset-v1/full\nT\\t1\tG1\tlab\n"))
})

test_that("analyses fingerprint the annotation actually passed, not the stored one", {
  x <- quiet(richR::buildAnnotFromList(list(T1 = c("G1", "G2"), T2 = "G2")))
  stored <- fingerprint_of(x)
  x@annot <- rbind(x@annot, data.frame(GeneID = "G3", Term = "T2", Annot = "T2"))
  edited_membership <- sha("richR-geneset-v1\nT1\tG1\nT1\tG2\nT2\tG2\nT2\tG3\n")
  edited_full <- sha("richR-geneset-v1/full\nT1\tG1\tT1\nT1\tG2\tT1\nT2\tG2\tT2\nT2\tG3\tT2\n")
  expect_false(identical(edited_membership, stored$membership))

  ora <- quiet(richR::enrich(c("G1", "G3"), x, pvalue = 1, padj = 1, minSize = 1,
                             minGSSize = 1))
  expect_equal(fingerprint_of(ora)$membership, edited_membership)
  expect_equal(fingerprint_of(ora)$full, edited_full)

  # richGO narrows to one ontology; the fingerprint covers the frame passed.
  go_frame <- data.frame(GeneID = c("G1", "G2", "G2", "G4"),
                         GOALL = c("T1", "T1", "T2", "T3"),
                         ONTOLOGYALL = c("BP", "BP", "BP", "MF"),
                         Annot = c("alpha", "alpha", "beta", "delta"))
  go <- quiet(richR::richGO(c("G1", "G2"), go_frame, ontology = "BP", pvalue = 1,
                            padj = 1, minSize = 1, minGSSize = 1))
  expect_equal(fingerprint_of(go)$membership,
               sha("richR-geneset-v1\nT1\tG1\nT1\tG2\nT2\tG2\nT3\tG4\n"))
  expect_equal(fingerprint_of(go)$full,
               sha("richR-geneset-v1/full\nT1\tG1\talpha\nT1\tG2\talpha\nT2\tG2\tbeta\nT3\tG4\tdelta\n"))

  skip_if_not_installed("fgsea")
  scores <- c(G1 = 3, G2 = 1, G3 = -1, G4 = -2)
  gsea <- quiet(richR::richGSEA(scores, object = x, minSize = 1, maxSize = 20,
                                pvalue = 1, padj = 1, nproc = 1, seed = 42))
  expect_equal(fingerprint_of(gsea)$membership, edited_membership)
  expect_equal(fingerprint_of(gsea)$full, edited_full)
})
