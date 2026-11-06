.validate_gsea_inputs <- function(x, pvalue, padj, minSize, maxSize,
                                 padj.method,
                                 eps, scoreType, nproc, BPPARAM, seed) {
  scalar_size <- function(value, name, allow_infinite = FALSE) {
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        value < 1 || (!is.finite(value) && !allow_infinite) ||
        (is.finite(value) && value != floor(value))) {
      stop(name, " must be one positive whole-number scalar.", call. = FALSE)
    }
  }
  if (!is.atomic(x) || !is.numeric(x) || !is.null(dim(x)) || length(x) == 0L ||
      is.null(names(x)) || anyNA(x) || any(!is.finite(x))) {
    stop("GSEA stats must be a non-empty named numeric vector with finite numbers.",
         call. = FALSE)
  }
  ids <- names(x)
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("GSEA stats must have unique, non-missing, non-blank gene names.",
         call. = FALSE)
  }
  .check_cutoff(pvalue, "pvalue")
  .check_cutoff(padj, "padj", allow_null = TRUE)
  scalar_size(minSize, "minSize")
  scalar_size(maxSize, "maxSize", allow_infinite = TRUE)
  if (minSize > maxSize) {
    stop("Minimum GSEA size must not exceed its maximum size.", call. = FALSE)
  }
  if (length(padj.method) != 1L || is.na(padj.method) ||
      !padj.method %in% stats::p.adjust.methods) {
    stop("padj.method must be one of: ", paste(stats::p.adjust.methods, collapse = ", "),
         call. = FALSE)
  }
  if (!is.null(eps) && (length(eps) != 1L || !is.numeric(eps) ||
                        !is.finite(eps) || eps < 0)) {
    stop("eps must be one finite numeric scalar greater than or equal to 0.",
         call. = FALSE)
  }
  if (!is.null(scoreType) && (!is.character(scoreType) || length(scoreType) != 1L ||
                              is.na(scoreType) ||
                              !scoreType %in% c("std", "pos", "neg"))) {
    stop("scoreType must be one of: std, pos, neg.", call. = FALSE)
  }
  if (!is.null(nproc) && (length(nproc) != 1L || !is.numeric(nproc) ||
                          !is.finite(nproc) || nproc < 0 ||
                          nproc != floor(nproc) || nproc > .Machine$integer.max)) {
    stop("nproc must be one whole-number scalar between 0 and .Machine$integer.max.",
         call. = FALSE)
  }
  if (!is.null(BPPARAM) && !methods::is(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from BiocParallelParam.", call. = FALSE)
  }
  if (!is.null(nproc) && !is.null(BPPARAM)) {
    stop("Specify only one of nproc or BPPARAM.", call. = FALSE)
  }
  if (!is.null(seed) && (length(seed) != 1L || !is.numeric(seed) ||
                         !is.finite(seed) || seed < 0 ||
                         seed != floor(seed) || seed > .Machine$integer.max)) {
    stop("seed must be one whole-number scalar between 0 and .Machine$integer.max.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# ---- GSEA gene-set identity ---------------------------------------------------
# Gene sets are keyed by identifier, never by label: distinct identifiers that
# share a label stay distinct gene sets. The identifier and label columns follow
# the shared annotation identity rule (.annotation_identity_frame, R/utils.R).

#' Gene-set memberships keyed by identifier
#' @param object annotation data.frame or Annot object
#' @return named list of gene vectors, one per identifier
#' @keywords internal
#' @noRd
.gsea_memberships <- function(object) {
  frame <- .annotation_identity_frame(object)
  split(as.character(frame$GeneID), as.character(frame$id))
}

#' Gene-set memberships for the pathways of a GSEA result
#'
#' Keyed by identifier; a result computed before gene sets were keyed by
#' identifier stores the label in `pathway`, so a pathway absent from the
#' identifier-keyed sets takes the label-keyed membership instead.
#' @param object annotation data.frame or Annot object
#' @param ids the `pathway` column of the GSEA result
#' @return named list of gene vectors
#' @keywords internal
#' @noRd
.gsea_result_memberships <- function(object, ids) {
  sets <- .gsea_memberships(object)
  missing <- setdiff(as.character(ids), names(sets))
  if (length(missing)) {
    frame <- .annotation_identity_frame(object)
    by_label <- split(as.character(frame$GeneID), as.character(frame$label))
    sets <- c(sets, by_label[intersect(missing, names(by_label))])
  }
  sets
}

#' Labels for identifiers (the identifier itself where no label exists)
#' @keywords internal
#' @noRd
.gsea_terms <- function(ids, map) {
  ids <- as.character(ids)
  terms <- unname(map[ids])
  missing <- is.na(terms)
  terms[missing] <- ids[missing]
  terms
}

#' Labels of a GSEA result table: the `Term` column, or `pathway` for a legacy
#' result without `Term`
#' @keywords internal
#' @noRd
.gsea_result_terms <- function(result) {
  ids <- as.character(result$pathway)
  if (!"Term" %in% colnames(result)) {
    return(ids)
  }
  terms <- as.character(result$Term)
  missing <- is.na(terms) | !nzchar(trimws(terms))
  terms[missing] <- ids[missing]
  terms
}

#' Display labels that stay unique per identifier
#'
#' A label shared by several of the given identifiers gets the identifier
#' appended, so facets, legends and network nodes never merge distinct gene sets.
#' @param ids identifiers to label
#' @param id_pool,term_pool identifiers and their labels
#' @keywords internal
#' @noRd
.gsea_display_labels <- function(ids, id_pool, term_pool) {
  ids <- as.character(ids)
  terms <- as.character(term_pool)[match(ids, as.character(id_pool))]
  terms[is.na(terms)] <- ids[is.na(terms)]
  u_ids <- unique(ids)
  u_terms <- terms[match(u_ids, ids)]
  shared <- u_terms %in% u_terms[duplicated(u_terms)]
  u_terms[shared] <- paste0(u_terms[shared], " [", u_ids[shared], "]")
  u_terms[match(ids, u_ids)]
}

#' Canonical KEGG pathway identifier (5 digits), NA for a non-KEGG identifier
#' @keywords internal
#' @noRd
.kegg_canonical_id <- function(x) {
  x <- as.character(x)
  is_kegg <- !is.na(x) & grepl("^([a-z]{2,4})?[0-9]{5}$", x)
  out <- rep(NA_character_, length(x))
  out[is_kegg] <- sub("^[a-z]{2,4}", "", x[is_kegg])
  out
}

#' Resolve pathway selectors to identifiers
#'
#' Each selector resolves, in order, by exact identifier, by KEGG identifier
#' with the organism prefix stripped, then by a unique label. A label carried by
#' several identifiers, or a selector matching nothing, is an error.
#' @param selectors character vector of selectors
#' @param ids candidate identifiers
#' @param terms labels of `ids`
#' @param what selector name used in messages
#' @return character vector of resolved identifiers (unique, selector order)
#' @keywords internal
#' @noRd
.gsea_resolve_selectors <- function(selectors, ids, terms, what = "pathway") {
  ids <- as.character(ids)
  terms <- as.character(terms)
  canon_ids <- .kegg_canonical_id(ids)
  resolved <- character()
  for (selector in as.character(selectors)) {
    hit <- unique(ids[ids == selector])
    if (!length(hit)) {
      canon <- .kegg_canonical_id(selector)
      if (!is.na(canon)) hit <- unique(ids[!is.na(canon_ids) & canon_ids == canon])
    }
    if (!length(hit)) {
      hit <- unique(ids[!is.na(terms) & terms == selector])
    }
    if (length(hit) > 1L) {
      stop(what, " '", selector, "' matches several pathways: ",
           paste(hit, collapse = ", "), ". Use a pathway identifier.", call. = FALSE)
    }
    if (!length(hit)) {
      stop(what, " '", selector, "' was not found in the GSEA results.", call. = FALSE)
    }
    resolved <- c(resolved, hit)
  }
  unique(resolved)
}

#' KEGG level columns for GSEA pathways
#'
#' KEGG identifiers (5 digits, optionally organism-prefixed) join the bundled
#' `path` table on `ko`; any other value joins on the `Level3` name.
#' @keywords internal
#' @noRd
.gsea_kegg_levels <- function(pathway, path) {
  canon <- .kegg_canonical_id(pathway)
  idx <- ifelse(is.na(canon), match(as.character(pathway), path$Level3),
                match(canon, path$ko))
  levels <- path[idx, , drop = FALSE]
  rownames(levels) <- NULL
  levels
}

#' ORA-shaped network frame from a GSEA result, built by column name
#' @param result GSEA result table
#' @param sep gene separator for a list-column leadingEdge
#' @return data.frame with Annot, Term, Pvalue, Padj, Significant, GeneID
#' @keywords internal
#' @noRd
.gsea_network_frame <- function(result, sep = ",") {
  le <- result$leadingEdge
  if (is.list(le)) {
    le <- vapply(le, function(x) paste(x, collapse = sep), character(1))
  }
  ids <- as.character(result$pathway)
  data.frame(
    Annot = ids,
    Term = .gsea_display_labels(ids, ids, .gsea_result_terms(result)),
    Pvalue = result$pval,
    Padj = result$padj,
    Significant = result$size,
    GeneID = as.character(le),
    stringsAsFactors = FALSE
  )
}

#' Enrichment analysis for any type of annotation data
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param keytype gene ID type
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param KEGG a logical evaluating to TRUE or FALSE indicating whether KEGG GSEA were peformed or not.
#' @param minSize Minimal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea). Smaller
#'   gene sets are not tested.
#' @param maxSize Maximal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea). Larger
#'   gene sets are not tested.
#' @param padj.method p value adjust method
#' @param organism organism name
#' @param ontology ontology type
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @param eps optional fgsea convergence tolerance
#' @param scoreType optional fgsea score type (std, pos, or neg)
#' @param nproc optional fgsea worker count
#' @param BPPARAM optional BiocParallel backend
#' @param seed optional seed scoped to this GSEA call
#' @param annotation_provenance optional annotation provenance captured from an Annot object
#' @importFrom fgsea fgseaMultilevel
#' @author Kai Guo
richGSEA_internal<-function(x,object,keytype="",pvalue=0.05,padj=0.05,KEGG=FALSE,minSize=15,maxSize=500,
                            padj.method="BH",organism=NULL,ontology=NULL,table=TRUE,sep=",",
                            eps=NULL,scoreType=NULL,nproc=NULL,BPPARAM=NULL,seed=NULL,
                            annotation_provenance=NULL){
  .validate_gsea_inputs(x, pvalue, padj, minSize, maxSize, padj.method,
                        eps, scoreType, nproc, BPPARAM, seed)
  rng_summary <- list(seed = NULL, RNGkind = RNGkind(), backend = NULL)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    old_kind <- RNGkind()
    on.exit({
      kind_args <- as.list(old_kind)
      names(kind_args) <- c("kind", "normal.kind", "sample.kind")[seq_along(kind_args)]
      do.call(RNGkind, kind_args)
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
    rng_summary <- list(seed = as.integer(seed), RNGkind = RNGkind(), backend = NULL)
  }
  metadata_padj <- padj
  original_input_ids <- names(x)
  # Read before sort() drops it (SPEC-RICHR-005 REQ-006).
  mapping_record <- .mapping_for_input(x)
  x<-sort(x)
  # Memberships come from gene-identifier pairs; labels are attached afterwards.
  identity_frame <- .annotation_identity_frame(object)
  fingerprint <- .geneset_fingerprint(object)
  normalization <- .normalize_annotation(identity_frame, func_name = "richGSEA")
  label_map <- .annotation_label_map(identity_frame, func_name = "richGSEA")
  object <- normalization$data[, 1:2, drop = FALSE]
  normalization$counts$eligible_rows <- as.integer(nrow(object))
  normalization$counts$effective_rows <- normalization$counts$eligible_rows
  annod<-sf(object);
  backend_args <- list(pathways = annod, stats = x,
                       minSize = minSize, maxSize = maxSize)
  if (!is.null(eps)) backend_args$eps <- eps
  if (!is.null(scoreType)) backend_args$scoreType <- scoreType
  if (!is.null(nproc)) backend_args$nproc <- nproc
  if (!is.null(BPPARAM)) backend_args$BPPARAM <- BPPARAM
  res <- as.data.frame(do.call(fgseaMultilevel, backend_args))
  rng_summary$backend <- if (!is.null(BPPARAM)) class(BPPARAM) else if (!is.null(nproc)) {
    list(nproc = nproc)
  } else NULL
  tested_count <- nrow(res)
  res$padj <- stats::p.adjust(res$pval, method = padj.method)
  if(isTRUE(table)){
    res$leadingEdge<-unlist(lapply(res$leadingEdge, function(x)paste(gsub(" ","",x),collapse = sep,sep="")))
  }
  # Term (the label) goes directly after leadingEdge, so fgsea's columns keep
  # their positions and KEGG level columns follow it.
  res$Term <- .gsea_terms(res$pathway, label_map)
  le_pos <- match("leadingEdge", colnames(res))
  if (!is.na(le_pos)) {
    cols <- colnames(res)
    res <- res[, c(cols[seq_len(le_pos)], "Term",
                   setdiff(cols[-seq_len(le_pos)], "Term")), drop = FALSE]
  }
  # Tested record (REQ-002): every pathway fgsea returned, captured before the
  # significance predicate; it takes the same ordering and KEGG level columns
  # as @result below, and nothing is recomputed.
  tested <- res
  # richR's own padj (above) is computed over the full returned family before
  # the significance predicate runs.
  res <- res[.sig_keep(res, pvalue, padj, p_col = "pval", padj_col = "padj"), ,
             drop = FALSE]
  if (is.null(padj)) padj <- numeric()
  res<-res[order(res$pval),]
  res <- res[!is.na(res$pathway),]
  tested <- tested[order(tested$pval),]
  tested <- tested[!is.na(tested$pathway),]
  if(isTRUE(KEGG)){
    path_env <- new.env(parent = emptyenv())
    utils::data("path", package = "richR", envir = path_env)
    path <- path_env$path
    res<-cbind(res,.gsea_kegg_levels(res$pathway, path))
    tested <- cbind(tested, .gsea_kegg_levels(tested$pathway, path))
  }
  tested$reported <- tested$pathway %in% res$pathway
  if(is.null(organism)){
    organism=character()
  }
  if(is.null(ontology)){
    ontology=character()
  }
  result<-new("GSEAResult",
              result=res,
              pvalueCutoff   = pvalue,
              pAdjustMethod  = padj.method,
              padjCutoff   = padj,
              genenumber    = length(x),
              organism       = organism,
              gene           = names(x),
              input = x,
              ontology = ontology,
              keytype        = keytype,
              sep = sep
  )
  result@tested <- tested
  annotation_genes <- unique(as.character(object[, 1]))
  ranked_ids <- unique(as.character(original_input_ids[!is.na(original_input_ids) &
                                                        nzchar(trimws(original_input_ids))]))
  input_info <- list(
    original_count = as.integer(length(original_input_ids)),
    valid_count = as.integer(sum(!is.na(original_input_ids) & nzchar(trimws(original_input_ids)))),
    unique_count = as.integer(length(ranked_ids)),
    duplicate_count = as.integer(sum(!is.na(original_input_ids) & nzchar(trimws(original_input_ids))) -
                                   length(ranked_ids)),
    invalid_count = as.integer(sum(is.na(original_input_ids) | !nzchar(trimws(original_input_ids)))),
    ids = ranked_ids, mapped_ids = ranked_ids[ranked_ids %in% annotation_genes],
    unmapped_ids = ranked_ids[!ranked_ids %in% annotation_genes],
    mapped_count = as.integer(sum(ranked_ids %in% annotation_genes)),
    unmapped_count = as.integer(sum(!ranked_ids %in% annotation_genes)),
    effective_ids = ranked_ids[ranked_ids %in% annotation_genes],
    effective_count = as.integer(sum(ranked_ids %in% annotation_genes)),
    excluded_by_universe_ids = character(),
    occurrences = rep(1L, length(ranked_ids))
  )
  input_info$mapping <- mapping_record
  background <- list(mode = "ranked_input", original_count = NA_integer_,
                     valid_count = NA_integer_, unique_count = NA_integer_,
                     invalid_count = NA_integer_, duplicate_count = NA_integer_,
                     unmapped_count = NA_integer_, unmapped_ids = character(),
                     effective_count = as.integer(length(ranked_ids)),
                     effective_ids = ranked_ids)
  testing <- list(family_policy = "fgsea_returned_before_filtering",
                  tested_count = as.integer(tested_count),
                  zero_overlap_count = 0L,
                  effective_N = NA_integer_, effective_n = NA_integer_)
  options <- list(pvalue = pvalue, padj = metadata_padj, minSize = minSize,
                  maxSize = maxSize, padj.method = padj.method,
                  ontology = ontology, keytype = keytype, KEGG = KEGG,
                  table = table, sep = sep, eps = eps, scoreType = scoreType,
                  nproc = nproc,
                  BPPARAM = if (is.null(BPPARAM)) NULL else class(BPPARAM),
                  seed = seed)
  if (is.null(annotation_provenance)) {
    annotation_provenance <- .annotation_provenance(
      species = organism, keytype = keytype, anntype = ontology,
      source_type = "custom", source_name = NA_character_
    )
  }
  # Computed from the frame passed, replacing any copy stored on the Annot (REQ-017).
  annotation_provenance$fingerprint <- fingerprint
  attr(result, "richR_metadata") <- .analysis_metadata(
    "GSEA", input_info, normalization$counts, background, testing, options,
    annotation_provenance, fgsea_version = .package_version_or_na("fgsea"),
    rng = rng_summary
  )
  result
}

#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object data.frame annotation
#' @param keytype gene ID type
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off
#' @param KEGG logical, whether KEGG GSEA were performed
#' @param minSize Minimal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea)
#' @param ontology ontology type
#' @param maxSize Maximal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea)
#' @param padj.method pvalue adjust method
#' @param organism organism name
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes
#' @param eps optional fgsea convergence tolerance
#' @param scoreType optional fgsea score type (std, pos, or neg)
#' @param nproc optional fgsea worker count
#' @param BPPARAM optional BiocParallel backend
#' @param seed optional seed scoped to this GSEA call
#' @details Gene sets are keyed by identifier, so distinct identifiers that
#'   share a label are tested as distinct gene sets. The identifier column is
#'   the first of \code{GOALL}, \code{PATH}, \code{Module} or \code{Term}
#'   present in \code{object}; an annotation carrying none of these names uses
#'   its second column. The label is the column named \code{Annot}.
#'   \code{result$pathway} holds the identifier and \code{result$Term}
#'   (directly after \code{leadingEdge}) the label: the identifier itself when
#'   no label exists, and the first non-missing label in annotation row order
#'   when an identifier carries several labels (with a warning naming those
#'   identifiers). With \code{KEGG = TRUE}, KEGG identifiers (five digits,
#'   optionally organism-prefixed) are joined to the KEGG levels by identifier,
#'   other values by the \code{Level3} name.
#' @section Identifier convention:
#' \code{Term} (like \code{GOALL}, \code{PATH} and \code{Module}) is the
#' identifier column in an annotation frame and \code{Annot} holds the label;
#' in a GSEA result, \code{pathway} is the identifier and \code{Term} the
#' label (in an ORA result, \code{Annot} and \code{Term}).
#' @export
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsako<-as.data.frame(hsako)
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' }
#' @author Kai Guo
setMethod("richGSEA", signature(object = "data.frame"),definition = function(x,object,keytype="",pvalue=0.05,padj=0.05,KEGG=FALSE,minSize=15,ontology="",
                                                                             maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",",
                                                                             eps=NULL,scoreType=NULL,nproc=NULL,BPPARAM=NULL,seed=NULL) {
  richGSEA_internal(x,object,keytype=keytype,pvalue=pvalue,padj=padj,KEGG=KEGG,minSize=minSize,ontology=ontology,
                    maxSize=maxSize,padj.method=padj.method,organism=organism,table=table,sep=sep,
                    eps=eps,scoreType=scoreType,nproc=nproc,BPPARAM=BPPARAM,seed=seed)
})
#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object Annot object
#' @param keytype gene ID type
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off
#' @param KEGG logical, whether KEGG GSEA were performed
#' @param minSize Minimal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea)
#' @param ontology ontology type
#' @param maxSize Maximal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea)
#' @param padj.method pvalue adjust method
#' @param organism organism name
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' }
#' @section Identifier convention:
#' \code{Term} (like \code{GOALL}, \code{PATH} and \code{Module}) is the
#' identifier column in an annotation frame and \code{Annot} holds the label;
#' in a GSEA result, \code{pathway} is the identifier and \code{Term} the
#' label (in an ORA result, \code{Annot} and \code{Term}).
#' @export
#' @author Kai Guo
#' @param eps optional fgsea convergence tolerance
#' @param scoreType optional fgsea score type (std, pos, or neg)
#' @param nproc optional fgsea worker count
#' @param BPPARAM optional BiocParallel backend
#' @param seed optional seed scoped to this GSEA call
#' @details Gene sets are keyed by identifier and labelled afterwards; see the
#'   \code{data.frame} method for the identifier and label rules.
setMethod("richGSEA", signature(object = "Annot"),definition = function(x,object,keytype="",pvalue=0.05,padj=0.05,KEGG=FALSE,minSize=15,ontology="",
                                                                        maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",",
                                                                        eps=NULL,scoreType=NULL,nproc=NULL,BPPARAM=NULL,seed=NULL) {
  provenance <- .annotation_provenance(object)
  richGSEA_internal(x,object@annot,keytype=object@keytype,pvalue=pvalue,padj=padj,KEGG=KEGG,minSize=minSize,ontology=object@anntype,
                    maxSize=maxSize,padj.method=padj.method,organism=object@species,table=table,sep=sep,
                    annotation_provenance=provenance,eps=eps,scoreType=scoreType,
                    nproc=nproc,BPPARAM=BPPARAM,seed=seed)
})


#' @name richGSEAplot
#' @title plot the gsea result
#' @param object Annot object
#' @param term the significant term
#' @param gseaRes GSEAResult object
#' @param top number of top pathways to display
#' @param default whether to use default table plot or individual enrichment plots
#' @importFrom fgsea plotEnrichment plotGseaTable
#' @importFrom ggplot2 ggtitle
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' richGSEAplot(term = res$pathway, object = hsako, gseaRes = res, default = FALSE)
#' }
#' @export
#' @author Kai Guo
richGSEAplot <- function(term, object, gseaRes, top = 10, default = TRUE) {
  # Input validation
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  # Extract fold change vector from GSEAResult object
  x <- gseaRes@input
  x <- sort(x)

  gseaRes_result <- gseaRes@result

  if(nrow(gseaRes_result) == 0) {
    stop("No pathways found in GSEA results")
  }

  if(nrow(gseaRes_result) < top) {
    top <- nrow(gseaRes_result)
  }

  # Gene sets keyed by identifier (never by label); titles show the Term label.
  # `term` resolves like every other GSEA selector: identifier, KEGG identifier
  # with or without organism prefix, or a unique label.
  annod <- .gsea_result_memberships(object, gseaRes_result$pathway)
  titles <- .gsea_result_terms(gseaRes_result)
  term <- .gsea_resolve_selectors(term, gseaRes_result$pathway, titles, what = "term")

  if(length(term) > 1 & !is.null(gseaRes_result)) {
    if(isTRUE(default)) {
      plotGseaTable(annod[term], stats = x, gseaRes_result, gseaParam = 0.5)
    } else {
      res <- lapply(gseaRes_result$pathway[seq_len(top)], function(y) {
        plotEnrichment(annod[[y]], stats = x) +
          ggtitle(titles[match(y, gseaRes_result$pathway)])
      })
      if (!requireNamespace("cowplot", quietly = TRUE))
        stop("Package 'cowplot' is required for grid GSEA plots. Install it with install.packages('cowplot').")
      cowplot::plot_grid(plotlist = res, ncol = 5)
    }
  } else {
    plotEnrichment(annod[[term]], stats = x)
  }
}

#' plot multiple significant pathways with enhanced customization and directional coloring
#' @importFrom ggplot2 ggplot geom_hline aes geom_point geom_segment geom_line theme_bw element_blank theme scale_color_manual ggtitle element_text annotate ggsave guides guide_legend
#' @param object Annot object
#' @param gseaRes GSEAResult object
#' @param mycol a vector indicate the colors used for the figure (used when show_direction = FALSE)
#' @param up_color color for up-regulated pathways (positive NES)
#' @param down_color color for down-regulated pathways (negative NES)
#' @param show_direction logical, whether to use directional coloring based on NES
#' @param pathways pathways to display (if NULL, use top/pvalue/padj filtering):
#'   pathway identifiers, KEGG identifiers with or without organism prefix, or
#'   unique \code{Term} labels; a label shared by several pathways, or a value
#'   matching nothing, is an error
#' @param pathway_pattern regex pattern matched against pathway identifiers and \code{Term} labels
#' @param category pathway category filter
#' @param min_genes minimum gene set size filter
#' @param max_genes maximum gene set size filter
#' @param top number of terms you want to display (used when pathways is NULL)
#' @param pvalue raw p-value cutoff, used only when \code{pathways} and
#'   \code{pathway_pattern} are NULL: a pathway is drawn when
#'   \code{pval <= pvalue} and, unless \code{padj} is NULL,
#'   \code{padj <= padj}. The explicit \code{pathways} and
#'   \code{pathway_pattern} selections are not significance-filtered by design;
#'   they choose among the rows the \code{GSEAResult} already holds.
#' @param padj adjusted p-value cutoff (used only when \code{pathways} and
#'   \code{pathway_pattern} are NULL); NULL filters on the raw p-value only
#' @param gseaParam GSEA parameter for enrichment score calculation
#' @param ticksSize size of the tick marks
#' @param plot_title title for the plot
#' @param show_pathway_names whether to show pathway names on plot
#' @param pathway_name_size size of pathway name text
#' @param line_alpha transparency of lines
#' @param point_alpha transparency of points
#' @param legend_position position of legend ("top", "bottom", "left", "right", "none")
#' @param legend_ncol number of columns in legend
#' @param return_data whether to return plot data along with plot
#' @param interactive whether to create interactive plot using plotly
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' # Plot with directional coloring (default)
#' richGSEAcurve(object = hsako, gseaRes = res, show_direction = TRUE)
#' # Plot specific pathways
#' richGSEAcurve(object = hsako, gseaRes = res,
#'          pathways = c("MAPK signaling pathway", "PI3K-Akt signaling pathway"))
#' # Plot pathways matching pattern
#' richGSEAcurve(object = hsako, gseaRes = res, pathway_pattern = "signaling")
#' # Traditional multi-color approach
#' richGSEAcurve(object = hsako, gseaRes = res, show_direction = FALSE)
#' }
#' @export
#' @author Kai Guo
richGSEAcurve <- function(object, gseaRes, mycol = NULL,
                     up_color = "#E31A1C", down_color = "#1F78B4", show_direction = TRUE,
                     pathways = NULL, pathway_pattern = NULL, category = NULL,
                     min_genes = NULL, max_genes = NULL, top = 10,
                     pvalue = 0.05, padj = 0.05, gseaParam = 1, ticksSize = 0.2,
                     plot_title = "GSEA Enrichment Plot", show_pathway_names = FALSE,
                     pathway_name_size = 10, line_alpha = 0.8, point_alpha = 0.6,
                     legend_position = "bottom", legend_ncol = 2,
                     return_data = FALSE, interactive = FALSE) {

  # Input validation
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  # Extract fold change vector from GSEAResult object
  x <- gseaRes@input

  # Extract results from GSEAResult object
  gseaRes_result <- gseaRes@result

  if(nrow(gseaRes_result) == 0) {
    stop("No pathways found in GSEA results")
  }

  if(!is.null(pathways) && !is.character(pathways)) {
    stop("pathways must be a character vector")
  }

  # Cutoffs are validated whichever selection branch runs (explicit
  # selections are not significance-filtered, but the arguments are checked).
  .check_cutoff(pvalue, "pvalue")
  .check_cutoff(padj, "padj", allow_null = TRUE)

  result_terms <- .gsea_result_terms(gseaRes_result)

  # Select pathways to display (identifiers)
  if(!is.null(pathways)) {
    # Resolve user selectors: identifier, prefix-stripped KEGG identifier, unique label
    sigpathway <- .gsea_resolve_selectors(pathways, gseaRes_result$pathway,
                                          result_terms, what = "pathways")
    # Get corresponding NES values
    pathway_indices <- match(sigpathway, gseaRes_result$pathway)
    sig_nes <- gseaRes_result$NES[pathway_indices]
  } else if(!is.null(pathway_pattern)) {
    # Use pattern-based selection on identifiers and labels
    pattern_indices <- grepl(pathway_pattern, gseaRes_result$pathway, ignore.case = TRUE) |
      grepl(pathway_pattern, result_terms, ignore.case = TRUE)
    sigpathway <- gseaRes_result$pathway[pattern_indices]
    sig_nes <- gseaRes_result$NES[pattern_indices]
    if(length(sigpathway) == 0) {
      stop("No pathways match the specified pattern")
    }
  } else {
    # Use filtering criteria: the significance predicate (raw AND adjusted)
    sig_indices <- .sig_keep(gseaRes_result, pvalue, padj,
                             p_col = "pval", padj_col = "padj")

    sigpathway <- gseaRes_result$pathway[sig_indices]
    sig_nes <- gseaRes_result$NES[sig_indices]
    if(length(sigpathway) == 0) {
      stop("richGSEAcurve: no pathway passes the cutoffs (pvalue = ", pvalue,
           ", padj = ", if (is.null(padj)) "NULL" else padj, ").", call. = FALSE)
    }

    if(top > length(sigpathway)) {
      top <- length(sigpathway)
    }
    sigpathway <- sigpathway[seq_len(top)]
    sig_nes <- sig_nes[seq_len(top)]
  }

  # Apply gene set size filtering
  if(!is.null(min_genes) || !is.null(max_genes)) {
    pathway_indices <- match(sigpathway, gseaRes_result$pathway)
    valid_mask <- rep(TRUE, length(pathway_indices))

    if(!is.null(min_genes)) {
      valid_mask <- valid_mask & (gseaRes_result$size[pathway_indices] >= min_genes)
    }
    if(!is.null(max_genes)) {
      valid_mask <- valid_mask & (gseaRes_result$size[pathway_indices] <= max_genes)
    }

    sigpathway <- sigpathway[valid_mask]
    sig_nes <- sig_nes[valid_mask]
  }

  if(length(sigpathway) == 0) {
    stop("No pathways remaining after filtering")
  }

  # From here on `sig_ids` holds the identifiers and `sigpathway` the display
  # labels (Term, unique per identifier) used for facets, colours and legends.
  sig_ids <- sigpathway
  sigpathway <- .gsea_display_labels(sig_ids, gseaRes_result$pathway, result_terms)

  # Set up colors based on direction preference
  if(show_direction) {
    # Create direction-based colors for individual pathways
    pathway_direction <- ifelse(sig_nes > 0, "Up-regulated", "Down-regulated")

    # Generate different shades/variants for each pathway within direction
    up_indices <- which(sig_nes > 0)
    down_indices <- which(sig_nes < 0)

    # Create color variations for up-regulated pathways
    if(length(up_indices) > 0) {
      up_colors <- colorRampPalette(c(up_color,
                                      adjustcolor(up_color, alpha.f = 0.7)))(length(up_indices))
    } else {
      up_colors <- character(0)
    }

    # Create color variations for down-regulated pathways
    if(length(down_indices) > 0) {
      down_colors <- colorRampPalette(c(down_color,
                                        adjustcolor(down_color, alpha.f = 0.7)))(length(down_indices))
    } else {
      down_colors <- character(0)
    }

    # Assign colors to pathways
    pathway_colors <- character(length(sigpathway))
    pathway_colors[up_indices] <- up_colors
    pathway_colors[down_indices] <- down_colors
    names(pathway_colors) <- sigpathway

    legend_title <- "Pathways (by Direction)"

  } else {
    # Use traditional multi-color approach
    if(is.null(mycol)) {
      mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9",
                 "#7A142C", "#5D90BA", "#431A3D", "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D",
                 "#7CC767")
    }
    pathway_colors <- rep(mycol, length.out = length(sigpathway))
    names(pathway_colors) <- sigpathway
    legend_title <- "Pathway"
    pathway_direction <- rep("Individual", length(sigpathway))
  }

  # Calculate GSEA curves
  fc <- x
  res <- lapply(seq_along(sig_ids), function(i) {
    curve <- .calGSEA(object, sig_ids[i], fc, gseaParam = gseaParam, ticksSize = ticksSize)
    curve$toPlot$pathway <- sig_ids[i]
    curve$toPlot$Group <- sigpathway[i]
    curve$pathway$pathway <- sig_ids[i]
    curve$pathway$Group <- sigpathway[i]
    curve
  })

  toPlot <- do.call(rbind, lapply(res, '[[', 'toPlot'))
  pathway_segments <- do.call(rbind, lapply(res, '[[', 'pathway'))
  tops <- unlist(lapply(res, '[[', 'tops'))
  bottoms <- unlist(lapply(res, '[[', 'bottoms'))

  # Keep pathways separate - each pathway keeps its individual identity
  # Add direction information for potential grouping in legend
  if(show_direction) {
    direction_mapping <- data.frame(
      Group = sigpathway,
      Direction = pathway_direction,
      stringsAsFactors = FALSE
    )
    toPlot <- merge(toPlot, direction_mapping, by = "Group", all.x = TRUE)
    pathway_segments <- merge(pathway_segments, direction_mapping, by = "Group", all.x = TRUE)
  }

  diff <- (max(tops) - min(bottoms))/8

  # Create the plot - each pathway maintains its individual line
  p <- ggplot(toPlot, aes(x = x, y = y, color = Group)) +
    geom_point(size = 0.1, alpha = point_alpha) +
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") +
    geom_line(alpha = line_alpha) +
    theme_bw() +
    geom_segment(data = pathway_segments,
                 mapping = aes(x = x, y = -diff/4, xend = x, yend = diff/4, color = Group),
                 linewidth = ticksSize) +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = legend_position,
          plot.title = element_text(hjust = 0.5, size = 14)) +
    labs(x = "Rank", y = "Enrichment Score") +
    ggtitle(plot_title)

  # Customize legend based on approach
  if(legend_position != "none") {
    if(show_direction) {
      # Create pathway labels with NES values and group by direction
      pathway_indices <- match(sig_ids, gseaRes_result$pathway)
      nes_values <- round(gseaRes_result$NES[pathway_indices], 2)

      # Create labels with pathway names and NES values
      pathway_labels_with_nes <- paste0(sigpathway, " (", nes_values, ")")

      # Truncate long pathway names but keep NES values
      max_pathway_length <- 40  # Adjust as needed
      truncated_pathways <- ifelse(nchar(sigpathway) > max_pathway_length,
                                   paste0(substr(sigpathway, 1, max_pathway_length-3), "..."),
                                   sigpathway)
      pathway_labels_with_nes <- paste0(truncated_pathways, " (", nes_values, ")")

      # Group pathways by direction
      up_indices <- which(sig_nes > 0)
      down_indices <- which(sig_nes < 0)

      # Create ordered labels (up-regulated first, then down-regulated)
      ordered_pathways <- character(0)
      ordered_colors <- character(0)
      ordered_labels <- character(0)

      if(length(up_indices) > 0) {
        # Sort up-regulated by NES (highest first)
        up_order <- up_indices[order(sig_nes[up_indices], decreasing = TRUE)]
        ordered_pathways <- c(ordered_pathways, sigpathway[up_order])
        ordered_colors <- c(ordered_colors, pathway_colors[sigpathway[up_order]])
        ordered_labels <- c(ordered_labels, pathway_labels_with_nes[up_order])
      }

      if(length(down_indices) > 0) {
        # Sort down-regulated by NES (least negative first, i.e., highest to lowest)
        down_order <- down_indices[order(sig_nes[down_indices], decreasing = TRUE)]
        ordered_pathways <- c(ordered_pathways, sigpathway[down_order])
        ordered_colors <- c(ordered_colors, pathway_colors[sigpathway[down_order]])
        ordered_labels <- c(ordered_labels, pathway_labels_with_nes[down_order])
      }

      # Add group separators in legend
      if(length(up_indices) > 0 && length(down_indices) > 0) {
        legend_title_with_groups <- paste0(legend_title, "\n",
                                           "Up-regulated (", length(up_indices), ") | ",
                                           "Down-regulated (", length(down_indices), ")")
      } else if(length(up_indices) > 0) {
        legend_title_with_groups <- paste0(legend_title, "\nUp-regulated (", length(up_indices), ")")
      } else {
        legend_title_with_groups <- paste0(legend_title, "\nDown-regulated (", length(down_indices), ")")
      }

      # Apply the colors and labels - single scale_color_manual call
      p <- p + scale_color_manual(values = setNames(ordered_colors, ordered_pathways),
                                  breaks = ordered_pathways,
                                  labels = ordered_labels,
                                  name = legend_title_with_groups) +
        guides(color = guide_legend(
          override.aes = list(size = 2, alpha = 1, linewidth = 1),
          ncol = legend_ncol,
          title.theme = element_text(size = 10, hjust = 0),
          label.theme = element_text(size = 8),
          keywidth = unit(1, "cm"),
          keyheight = unit(0.5, "cm")
        ))

    } else {
      # Traditional approach - single scale_color_manual call
      legend_labels <- names(pathway_colors)
      if(max(nchar(legend_labels)) > 50) {
        legend_labels <- ifelse(nchar(legend_labels) > 50,
                                paste0(substr(legend_labels, 1, 47), "..."),
                                legend_labels)
        names(pathway_colors) <- legend_labels
      }

      p <- p + scale_color_manual(values = pathway_colors, name = legend_title) +
        guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                    ncol = legend_ncol,
                                    title = legend_title))
    }
  } else {
    # No legend case - still need to set colors
    p <- p + scale_color_manual(values = pathway_colors, guide = "none")
  }

  # Add pathway names as annotations if requested
  if(show_pathway_names && length(sigpathway) <= 10) {
    annotation_y <- seq(max(tops) * 0.8, min(bottoms) * 0.8, length.out = length(sigpathway))
    p <- p + annotate("text", x = max(toPlot$x) * 0.7, y = annotation_y,
                      label = sigpathway, size = pathway_name_size/3, hjust = 0, vjust = 0.5)
  }

  # Create interactive plot if requested
  if(interactive) {
    if(!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package not available, returning static plot")
    } else {
      p <- plotly::ggplotly(p)
    }
  }

  # Return data if requested
  if(return_data) {
    return(list(
      plot = p,
      data = toPlot,
      pathways = stats::setNames(sig_ids, sigpathway),
      directions = if(show_direction) pathway_direction else NULL,
      nes_values = sig_nes,
      pathway_colors = pathway_colors,
      stats = list(tops = tops, bottoms = bottoms),
      pathway_data = pathway_segments
    ))
  } else {
    return(p)
  }
}

#' Keep GSEA result rows passing the significance predicate
#'
#' \code{pval_cutoff} and \code{padj_cutoff} are validated with the shared rule;
#' \code{padj_cutoff = NULL} filters on the raw p-value only, and
#' \code{pval_cutoff = NULL} applies no raw cutoff.
#' @param result GSEA result data.frame (\code{pval}, \code{padj} columns)
#' @param pval_cutoff raw p-value cutoff or NULL
#' @param padj_cutoff adjusted p-value cutoff or NULL
#' @keywords internal
#' @noRd
.gsea_sig_filter <- function(result, pval_cutoff, padj_cutoff) {
  .check_cutoff(pval_cutoff, "pval_cutoff", allow_null = TRUE)
  .check_cutoff(padj_cutoff, "padj_cutoff", allow_null = TRUE)
  if (is.null(pval_cutoff)) pval_cutoff <- 1
  result[.sig_keep(result, pval_cutoff, padj_cutoff, p_col = "pval", padj_col = "padj"), ,
         drop = FALSE]
}

#' Get available pathways from GSEA results
#' @param gseaRes GSEAResult object
#' @param padj_cutoff adjusted p-value cutoff in (0, 1]: a pathway is kept when
#'   \code{padj <= padj_cutoff}; NULL applies no adjusted cutoff
#' @param pval_cutoff raw p-value cutoff in (0, 1], applied together with
#'   \code{padj_cutoff} (\code{pval <= pval_cutoff}); NULL applies no raw cutoff
#' @param min_size minimum pathway size
#' @param max_size maximum pathway size
#' @export
#' @examples
#' \dontrun{
#' pathways <- getPathways(gseaRes, padj_cutoff = 0.05)
#' }
#' @author Kai Guo
getPathways <- function(gseaRes, padj_cutoff = 0.05, pval_cutoff = 0.05, min_size = NULL, max_size = NULL) {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  result <- gseaRes@result

  # Apply significance filtering
  result <- .gsea_sig_filter(result, pval_cutoff, padj_cutoff)

  # Apply size filtering
  if(!is.null(min_size)) {
    result <- result[result$size >= min_size, ]
  }
  if(!is.null(max_size)) {
    result <- result[result$size <= max_size, ]
  }

  return(result$pathway)
}

#' Search pathways by keyword
#' @param gseaRes GSEAResult object
#' @param keyword search term (can be regex pattern)
#' @param ignore_case whether to ignore case when searching
#' @return pathway identifiers whose identifier or \code{Term} label matches,
#'   named by their \code{Term} label
#' @export
#' @examples
#' \dontrun{
#' signaling_pathways <- searchPathways(gseaRes, "signaling")
#' mapk_pathways <- searchPathways(gseaRes, "MAPK")
#' }
#' @author Kai Guo
searchPathways <- function(gseaRes, keyword, ignore_case = TRUE) {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  pathways <- as.character(gseaRes@result$pathway)
  terms <- .gsea_result_terms(gseaRes@result)
  hit <- grepl(keyword, pathways, ignore.case = ignore_case) |
    grepl(keyword, terms, ignore.case = ignore_case)
  # identifiers, named by their Term label
  matched <- stats::setNames(pathways[hit], terms[hit])

  if(length(matched) == 0) {
    message("No pathways found matching the keyword: ", keyword)
  }

  return(matched)
}

#' Save GSEA plot with optimal dimensions
#' @param plot_obj ggplot object from richGSEAcurve
#' @param filename output filename
#' @param width plot width in inches
#' @param height plot height in inches
#' @param dpi resolution in dots per inch
#' @param device output device (auto-detected from filename if NULL)
#' @importFrom ggplot2 ggsave
#' @export
#' @examples
#' \dontrun{
#' p <- richGSEAcurve(object = hsako, gseaRes = res)
#' saveGSEAplot(p, "gsea_plot.pdf")
#' saveGSEAplot(p, "gsea_plot.png", dpi = 600)
#' }
#' @author Kai Guo
saveGSEAplot <- function(plot_obj, filename, width = 12, height = 8, dpi = 300, device = NULL) {
  if(!inherits(plot_obj, c("gg", "ggplot"))) {
    stop("plot_obj must be a ggplot object")
  }

  ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, device = device)
  message("Plot saved to: ", filename)
}

#' Create multiple GSEA plots for different comparisons
#' @param gsea_list named list of GSEAResult objects
#' @param object Annot object
#' @param output_dir directory to save plots
#' @param width plot width
#' @param height plot height
#' @param dpi resolution
#' @param format output format (pdf, png, etc.)
#' @param ... additional arguments passed to richGSEAcurve
#' @export
#' @examples
#' \dontrun{
#' gsea_results <- list(
#'   "D7_WT_vs_TFAM" = res1,
#'   "D14_WT_vs_D7_WT" = res2
#' )
#' plots <- batchGSEAplot(gsea_results, object = hsako)
#' }
#' @author Kai Guo
batchGSEAplot <- function(gsea_list, object, output_dir = "gsea_plots",
                          width = 12, height = 8, dpi = 300, format = "pdf", ...) {
  if(!is.list(gsea_list) || is.null(names(gsea_list))) {
    stop("gsea_list must be a named list")
  }

  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  plots <- lapply(names(gsea_list), function(name) {
    message("Processing: ", name)

    tryCatch({
      p <- richGSEAcurve(object = object, gseaRes = gsea_list[[name]],
                    plot_title = paste("GSEA:", name), ...)
      filename <- file.path(output_dir, paste0(name, "_gsea.", format))
      ggsave(filename, p, width = width, height = height, dpi = dpi)
      return(p)
    }, error = function(e) {
      warning("Failed to create plot for ", name, ": ", e$message)
      return(NULL)
    })
  })

  names(plots) <- names(gsea_list)
  plots <- plots[!sapply(plots, is.null)]  # Remove failed plots

  message("Successfully created ", length(plots), " plots in ", output_dir)
  return(plots)
}

#' Get pathway statistics summary
#' @param gseaRes GSEAResult object
#' @param top number of top pathways to summarize
#' @export
#' @examples
#' \dontrun{
#' summary <- getPathwayStats(gseaRes, top = 20)
#' print(summary)
#' }
#' @author Kai Guo
getPathwayStats <- function(gseaRes, top = 10) {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  result <- gseaRes@result

  if(nrow(result) == 0) {
    return(data.frame())
  }

  # Get top pathways
  top_pathways <- head(result, top)

  # Create summary
  summary_stats <- data.frame(
    pathway = top_pathways$pathway,
    Term = .gsea_result_terms(top_pathways),
    NES = round(top_pathways$NES, 3),
    pvalue = formatC(top_pathways$pval, format = "e", digits = 2),
    padj = formatC(top_pathways$padj, format = "e", digits = 2),
    size = top_pathways$size,
    direction = ifelse(top_pathways$NES > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )

  return(summary_stats)
}

#' Filter pathways by multiple criteria
#' @param gseaRes GSEAResult object
#' @param nes_cutoff absolute NES cutoff
#' @param padj_cutoff adjusted p-value cutoff in (0, 1]: a pathway is kept when
#'   \code{padj <= padj_cutoff}; NULL applies no adjusted cutoff
#' @param pval_cutoff raw p-value cutoff in (0, 1], applied together with
#'   \code{padj_cutoff} (\code{pval <= pval_cutoff}); NULL applies no raw cutoff
#' @param min_size minimum pathway size
#' @param max_size maximum pathway size
#' @param direction filter by direction ("up", "down", or "both")
#' @export
#' @examples
#' \dontrun{
#' filtered <- filterPathways(gseaRes, nes_cutoff = 1.5, padj_cutoff = 0.01, direction = "up")
#' }
#' @author Kai Guo
filterPathways <- function(gseaRes, nes_cutoff = NULL, padj_cutoff = 0.05, pval_cutoff = 0.05,
                           min_size = NULL, max_size = NULL, direction = "both") {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  result <- gseaRes@result

  # Apply NES filtering
  if(!is.null(nes_cutoff)) {
    result <- result[abs(result$NES) >= nes_cutoff, ]
  }

  # Apply significance filtering
  result <- .gsea_sig_filter(result, pval_cutoff, padj_cutoff)

  # Apply size filtering
  if(!is.null(min_size)) {
    result <- result[result$size >= min_size, ]
  }
  if(!is.null(max_size)) {
    result <- result[result$size <= max_size, ]
  }

  # Apply direction filtering
  if(direction == "up") {
    result <- result[result$NES > 0, ]
  } else if(direction == "down") {
    result <- result[result$NES < 0, ]
  }

  return(result)
}

#' GSEA Enrichment analysis function for data frame with log2FC and name specifc
#' @param x data.frame include the gene name and log2FC
#' @param object annotation file for all genes
#' @param log2FC the column name for log2FoldChange
#' @param gene.col the gene name column if rownames won't be used as gene name
#' @param simplify boolean value to indicate if the simple table should be returned
#' @param keytype gene ID type
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off
#' @param KEGG logical, whether KEGG GSEA were performed
#' @param minSize Minimal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea)
#' @param ontology ontology type
#' @param maxSize Maximal size of a gene set to test, counted as the gene-set
#'   membership intersected with the ranked genes (passed to fgsea)
#' @param padj.method p value adjust method
#' @param organism organism name
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @param eps optional fgsea convergence tolerance
#' @param scoreType optional fgsea score type (std, pos, or neg)
#' @param nproc optional fgsea worker count
#' @param BPPARAM optional BiocParallel backend
#' @param seed optional seed scoped to this GSEA call
#' @export
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsako<-as.data.frame(hsako)
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' d<-data.frame(Gene=name,log2FoldChange=gene)
#' rownames(d)<-d$Gene
#' res<-parGSEA(d,object = hsako)
#' }
#' @author Kai Guo
parGSEA<-function(x,object,log2FC="log2FoldChange",gene.col=NULL,
                  simplify=FALSE,keytype="",
                  pvalue=0.05,padj=0.05,KEGG=FALSE,minSize=15,ontology="",
                  maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",",
                  eps=NULL,scoreType=NULL,nproc=NULL,BPPARAM=NULL,seed=NULL){
  fc<-x[,log2FC]
  if(is.null(gene.col)){
    names(fc)<-rownames(x)
  }else{
    names(fc)<-x[,gene.col]
  }
  tmp<-richGSEA(fc,object,keytype=keytype,pvalue=pvalue,padj=padj,KEGG=KEGG,minSize=minSize,ontology=ontology,
                maxSize=maxSize,padj.method=padj.method,organism=organism,table=table,sep=sep,
                eps=eps,scoreType=scoreType,nproc=nproc,BPPARAM=BPPARAM,seed=seed)
  if(isTRUE(simplify)){
    tmp<-result(tmp)
  }
  return(tmp)
}

###
.calGSEA<-function (obj,sigpathway, fc, gseaParam = 1, ticksSize = 0.2){
  frame <- .annotation_identity_frame(obj)
  pathway <- split(as.character(frame$GeneID), as.character(frame$id))[[sigpathway]]
  if (is.null(pathway)) {
    # A result computed before gene sets were keyed by identifier stores the
    # label in `pathway`; its gene set was the label-keyed membership.
    pathway <- split(as.character(frame$GeneID), as.character(frame$label))[[sigpathway]]
  }
  stats <- fc
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                 returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  pathway<-data.frame(x = pathway,Group=sigpathway)
  toPlot$Group<-sigpathway
  return(list(toPlot=toPlot,pathway=pathway,tops=max(tops),bottoms=min(bottoms)))
}
