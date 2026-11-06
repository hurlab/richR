# SPEC-RICHR-014: all fixture rows have overlapping genes except the explicit control.
.review014_inputs <- function() {
  a <- .m61_fp3()
  b <- a; b@result$Pvalue <- c(0.02, 0.03, 0.04)
  b@result$Padj <- c(0.025, 0.035, 0.045)
  b@result$Significant <- c(12, 13, 14)
  # Repeated identifiers within a source must also remain separate.
  b@result$Annot[2] <- "A1"
  c <- a; c@result <- c@result[1:2, , drop = FALSE]
  c@result$Annot <- c("A1 [result 1]", "isolated")
  c@result$Term <- c("protected", "isolated")
  c@result$GeneID[2] <- "g30"
  list(a, b, c)
}

test_that("netmap qualifies colliding source IDs without pooling statistics", {
  inputs <- .review014_inputs(); before <- serialize(inputs, NULL)
  ids <- unlist(lapply(inputs, function(x) x@result$Annot), use.names = FALSE)
  sources <- rep(seq_along(inputs), vapply(inputs, function(x) nrow(x@result), 1L))
  dup <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
  ids[dup] <- tail(make.unique(c(ids, paste0(ids[dup], " [result ", sources[dup], "]"))), sum(dup))
  object <- do.call(rbind, lapply(inputs, function(x) x@result))
  object$Annot <- ids
  object <- object[order(object$Pvalue), ]
  rownames(object) <- object$Annot
  graphfile <- tempfile("richr014-graph-")
  on.exit(unlink(paste0(graphfile, ".graphml")), add = TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  for (use in c(TRUE, FALSE)) {
    set.seed(1)
    p <- richNetmap(inputs, top = 20, pvalue = 1, padj = NULL, useTerm = use, weightcut = 0,
                     writeCyt = TRUE, cytoscapeFile = graphfile)
    expect_identical(nrow(p$data), 7L)
    set.seed(1)
    ref <- richNetwork(object, gene = inputs[[1]]@gene, top = 20,
                       pvalue = 1, padj = NULL, useTerm = use, weightcut = 0)
    expect_equal(p$data, ref$data)
    expect_false("isolated" %in% p$data$label)
    if (!use) {
      expect_true("A1 [result 1]" %in% p$data$label)
      graph <- igraph::read_graph(paste0(graphfile, ".graphml"), format = "graphml")
      vertices <- igraph::as_data_frame(graph, what = "vertices")
      idx <- match(vertices$name, object$Annot)
      expect_false(anyNA(idx))
      expect_equal(vertices$size, log10(object$Significant[idx]) * 10)
      # Per-source adjusted p-values must survive, including differing colliding rows.
      evidence <- -log10(object$Padj[idx])
      colors <- richR:::.color_scale("red", "orange")
      expected <- colors[vapply(evidence, richR:::.getIdx, numeric(1), min(evidence), max(evidence))]
      expect_equal(vertices$color, unname(expected))
    }
  }
  expect_identical(serialize(inputs, NULL), before)
})

test_that("public network wrappers honor HTML saving only for interactive output", {
  dir <- tempfile("richr014-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  x <- .m61_fp3()
  # Native GSEA frame exercises its already-supported network conversion too.
  g <- methods::new("GSEAResult")
  g@result <- data.frame(pathway = x@result$Annot, Term = x@result$Term,
                        pval = x@result$Pvalue, padj = x@result$Padj,
                        size = x@result$Significant)
  g@result$leadingEdge <- strsplit(x@result$GeneID, ",", fixed = TRUE)
  g@gene <- x@gene; g@sep <- ","
  calls <- list(function(...) richNetwork(x, ...),
                function(...) richNetwork(x@result, gene = x@gene, ...),
                function(...) richNetwork(g, gene = x@gene, ...),
                function(...) richNetmap(list(x), ...))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_along(calls)) for (save in c(FALSE, TRUE)) for (interactive in c(TRUE, FALSE)) {
    path <- file.path(dir, paste(i, save, interactive, sep = "-"))
    p <- calls[[i]](pvalue = 1, padj = NULL, weightcut = 0, visNet = interactive,
                    savehtml = save, filename = path)
    expect_identical(file.exists(paste0(path, ".html")), save && interactive)
    if (interactive) expect_s3_class(p, "visNetwork") else expect_s3_class(p, "ggplot")
    if (save && interactive && file.exists(paste0(path, ".html")))
      expect_match(paste(readLines(paste0(path, ".html"), warn = FALSE), collapse = ""), "visNetwork")
  }
})
