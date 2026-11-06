#' Read enrichment diagnostics and provenance
#'
#' Metadata is recorded on newly created Annot, richResult, and GSEAResult
#' objects. Legacy objects return an explicit unavailable record.
#'
#' The \code{testing} element names the family that \code{padj.method}
#' adjusted over (\code{family_policy}) and its size (\code{tested_count}):
#' \describe{
#'   \item{\code{"size_bounded_positive_overlap"}}{ORA with
#'     \code{padj.family = "overlap"} (default): terms whose annotated size in
#'     the background lies in \code{[minGSSize, maxGSSize]} and that share at
#'     least one gene with the query.}
#'   \item{\code{"size_bounded_all_eligible"}}{ORA with
#'     \code{padj.family = "all"}: every term whose annotated size lies in
#'     \code{[minGSSize, maxGSSize]}, zero-overlap terms contributing p = 1.}
#'   \item{\code{"fgsea_returned_before_filtering"}}{GSEA: every pathway
#'     returned by fgsea, before the significance cutoffs.}
#' }
#' For ORA results, \code{options$padj.family} records the family argument.
#'
#' For results, \code{testing$zero_overlap_count} is the number of tested
#' terms with zero overlap that the tested table (\code{\link{testedTerms}})
#' does not list: non-zero only for ORA with \code{padj.family = "all"}, and
#' \code{0} otherwise. The tested table's row count plus
#' \code{zero_overlap_count} equals \code{tested_count}.
#'
#' \code{input$occurrences} (ORA: the number of times each identifier of
#' \code{input$ids} occurred in the valid input; GSEA: 1 per ranked gene) and
#' \code{input$mapping} are recorded by richR 0.2.0. \code{input$mapping} is
#' the \code{\link{mapGenes}()} record carried by the analysed vector when its
#' targets (and scores) match it, \code{list(available = FALSE, reason =
#' "mapping record does not match the submitted input")} when they do not, and
#' absent (\code{NULL}) when the input carried no record; read it with
#' \code{\link{mappingAudit}()}. For ORA results,
#' \code{options$universe_policy} records how \code{N} and \code{n} were
#' counted (\code{"annotated"} or \code{"measured"}, see \code{\link{enrich}});
#' metadata recorded without it describes the \code{"annotated"} behaviour.
#'
#' \code{schema_version} is \code{2} for \code{richResult} and
#' \code{GSEAResult} metadata recorded by richR 0.2.0, and \code{1} for
#' \code{Annot} metadata, for results created by earlier versions (which have
#' no \code{zero_overlap_count} and no tested table) and for the unavailable
#' record.
#' @param x Annot, richResult, or GSEAResult object
#' @return A named list containing diagnostics and provenance
#' @export
enrichmentInfo <- function(x) {
  supported <- methods::is(x, "Annot") || methods::is(x, "richResult") ||
    methods::is(x, "GSEAResult")
  if (!supported) {
    stop("enrichmentInfo() supports Annot, richResult, or GSEAResult objects.",
         call. = FALSE)
  }
  metadata <- attr(x, "richR_metadata", exact = TRUE)
  if (is.null(metadata)) {
    return(list(schema_version = 1L, available = FALSE,
                reason = "metadata not recorded", diagnostics = NULL,
                provenance = NULL))
  }
  metadata
}

.richr_version <- function() {
  if ("richR" %in% loadedNamespaces()) {
    return(as.character(utils::packageVersion("richR")))
  }
  NA_character_
}

.package_version_or_na <- function(package) {
  if (length(package) == 1L && !is.na(package) && nzchar(package) &&
      requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}

.annotation_provenance <- function(object = NULL, constructor = NULL,
                                   species = NULL, keytype = NULL,
                                   anntype = NULL, source_type = NULL,
                                   source_name = NULL, source_version = NA_character_,
                                   provider_package = NULL,
                                   provider_version = NULL,
                                   constructor_normalization = NULL,
                                   label_source = NULL) {
  if (methods::is(object, "Annot")) {
    previous <- attr(object, "richR_metadata", exact = TRUE)
    if (is.list(previous) && is.list(previous$provenance$annotation)) {
      annotation <- previous$provenance$annotation
      annotation$species <- object@species
      annotation$keytype <- object@keytype
      annotation$anntype <- object@anntype
      return(annotation)
    }
    species <- object@species
    keytype <- object@keytype
    anntype <- object@anntype
    source_type <- source_type %||% "unknown"
    source_name <- source_name %||% NA_character_
  }
  provenance <- list(
    constructor = constructor,
    species = species,
    keytype = keytype,
    anntype = anntype,
    source_type = source_type %||% "custom",
    source_name = source_name %||% NA_character_,
    source_version = source_version,
    provider_package = provider_package,
    provider_version = provider_version,
    constructor_normalization = constructor_normalization
  )
  ## Where the labels came from (KEGG / KEGGM), kept apart from the membership
  ## provider above so a bundled fallback is never read as a live retrieval.
  if (!is.null(label_source)) provenance$label_source <- label_source
  provenance
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.annotation_metadata <- function(normalization, provenance, options = list()) {
  list(schema_version = 1L, available = TRUE, analysis = "annotation",
       annotation = normalization$counts,
       provenance = list(richR_version = .richr_version(),
                         R_version = as.character(getRversion()),
                         fgsea_version = NULL,
                         annotation = provenance, rng = NULL),
       options = options)
}

#' Gene-set fingerprints of an annotation frame (canonical form richR-geneset-v1)
#'
#' Pairs are (identifier, gene) under the shared identity rule, restricted to
#' valid rows, de-duplicated and sorted in C-locale byte order; values are
#' UTF-8 and escaped (backslash, TAB, LF, CR). \code{membership} hashes
#' \code{id<TAB>gene} lines; \code{full} adds each identifier's label (first
#' non-blank \code{Annot} in row order, else the identifier). Both are SHA-256
#' of the text itself, never of a serialized R object.
#' @param annot annotation data.frame or Annot object
#' @return list(algorithm, canonicalization, membership, full, n_pairs, n_sets, n_genes)
#' @importFrom digest digest
#' @noRd
.geneset_fingerprint <- function(annot) {
  frame <- .annotation_identity_frame(annot)
  ids <- as.character(frame$id)
  genes <- as.character(frame$GeneID)
  labels <- as.character(frame$label)
  ## Blank = only spaces, TAB, CR or LF (what trimws() strips).
  blank <- function(x) grepl("^[ \t\r\n]*$", x, perl = TRUE)
  valid <- !is.na(ids) & !is.na(genes) & !blank(ids) & !blank(genes)
  ids <- enc2utf8(ids[valid])
  genes <- enc2utf8(genes[valid])
  labels <- labels[valid]
  ## Label of each identifier: its first non-blank label in row order (REQ-003).
  labelled <- which(!is.na(labels) & !blank(labels))
  first <- labelled[!duplicated(ids[labelled])]
  ## Unique pairs in C-locale byte order: sort, then drop repeats of the previous pair.
  o <- order(ids, genes, method = "radix")
  id <- ids[o]
  gene <- genes[o]
  n <- length(id)
  keep <- if (n > 1L) c(TRUE, id[-1L] != id[-n] | gene[-1L] != gene[-n]) else rep(TRUE, n)
  id <- id[keep]
  gene <- gene[keep]
  label <- labels[first][match(id, ids[first])]
  label[is.na(label)] <- id[is.na(label)]
  escape <- function(x) {
    if (!any(grepl("[\\\\\t\n\r]", x, perl = TRUE))) return(x)
    x <- gsub("\\", "\\\\", x, fixed = TRUE)
    x <- gsub("\t", "\\t", x, fixed = TRUE)
    x <- gsub("\n", "\\n", x, fixed = TRUE)
    gsub("\r", "\\r", x, fixed = TRUE)
  }
  id_text <- escape(id)
  gene_text <- escape(gene)
  label <- escape(enc2utf8(label))
  ## Hash the canonical text's bytes, written field by field into memory (much
  ## faster than paste(collapse =) on large annotations; same bytes).
  hash <- function(header, fields) {
    con <- rawConnection(raw(0), "wb")
    on.exit(close(con))
    writeLines(header, con, sep = "\n", useBytes = TRUE)
    if (length(id)) {
      writeLines(as.vector(do.call(rbind, fields)), con, sep = "", useBytes = TRUE)
    }
    digest(rawConnectionValue(con), algo = "sha256", serialize = FALSE)
  }
  m <- length(id)
  list(
    algorithm = "sha256", canonicalization = "richR-geneset-v1",
    membership = hash("richR-geneset-v1", list(id_text, "\t", gene_text, "\n")),
    full = hash("richR-geneset-v1/full",
                list(id_text, "\t", gene_text, "\t", label, "\n")),
    n_pairs = m,
    n_sets = if (m) 1L + sum(id[-1L] != id[-m]) else 0L,
    n_genes = length(unique(gene))
  )
}

.attach_annotation_metadata <- function(object, normalization, provenance,
                                        options = list()) {
  provenance$fingerprint <- .geneset_fingerprint(object@annot)
  attr(object, "richR_metadata") <- .annotation_metadata(
    normalization, provenance, options
  )
  object
}

.analysis_metadata <- function(analysis, input, annotation, background,
                               testing, options, annotation_provenance,
                               fgsea_version = NULL, rng = NULL,
                               labels = NULL) {
  ## Result metadata is schema 2 (testing$zero_overlap_count; the object carries
  ## a tested slot). Annot metadata and the unavailable record stay schema 1.
  metadata <- list(schema_version = 2L, available = TRUE, analysis = analysis,
       input = input, annotation = annotation, background = background,
       testing = testing, options = options,
       provenance = list(richR_version = .richr_version(),
                         R_version = as.character(getRversion()),
                         fgsea_version = fgsea_version,
                         annotation = annotation_provenance,
                         rng = rng))
  ## Label source used at analysis time (richKEGG), when the analysis looked labels up.
  if (!is.null(labels)) metadata$provenance$labels <- labels
  metadata
}
