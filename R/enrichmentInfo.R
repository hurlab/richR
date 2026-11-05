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

.attach_annotation_metadata <- function(object, normalization, provenance,
                                        options = list()) {
  attr(object, "richR_metadata") <- .annotation_metadata(
    normalization, provenance, options
  )
  object
}

.analysis_metadata <- function(analysis, input, annotation, background,
                               testing, options, annotation_provenance,
                               fgsea_version = NULL, rng = NULL,
                               labels = NULL) {
  metadata <- list(schema_version = 1L, available = TRUE, analysis = analysis,
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
