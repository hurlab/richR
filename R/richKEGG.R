#' KEGG Pathway Enrichment analysis function
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata GO annotation data
#' @param ontology KEGG
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimum annotated term size for testing: the number of the term's genes in the background. \code{[minGSSize, maxGSSize]} (inclusive) bounds the tested family before p-value adjustment.
#' @param maxGSSize maximum annotated term size for testing (inclusive upper bound of the tested family; see \code{minGSSize}).
#' @param keepRich keep terms with rich factor value equal 1 even when their overlap is below \code{minSize} (default: TRUE). It exempts only the \code{minSize} reporting filter and never admits a term outside the tested family.
#' @param padj.family family of terms that \code{padj.method} adjusts over:
#'   \code{"overlap"} (default) holds the size-eligible terms (see \code{minGSSize})
#'   that share at least one gene with the query; \code{"all"} holds every
#'   size-eligible term, zero-overlap terms contributing p = 1 (they are never
#'   reported). \code{minSize}, \code{maxSize}, \code{keepRich} and the cutoffs are
#'   applied afterwards and never change an adjusted value. The default family
#'   matches clusterProfiler/DOSE but is a convention, not an FDR guarantee: with a
#'   universe of 100 single-gene terms and one selected gene under the null (with
#'   \code{minGSSize = 1, minSize = 1}), the only overlapping term has p = 0.01 and
#'   BH over that one term gives 0.01, so a null term is reported every time; BH
#'   over all 100 gives 1. Use \code{padj.family = "all"} for adjustment over every
#'   size-eligible term.
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want to use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @param universe optional character vector defining the ORA background
#' @param annotation_provenance optional annotation provenance captured from an Annot object
#' @author Kai Guo
richKEGG_internal<-function(x,kodata,pvalue=0.05,padj=0.05,ontology="KEGG",
                            organism=NULL,keytype="SYMBOL",minSize=2,maxSize=500,
                            minGSSize = 10, maxGSSize = 500,
                            keepRich=TRUE, filename=NULL,padj.method="BH",builtin=TRUE,sep=",",universe=NULL,
                            annotation_provenance=NULL,padj.family=c("overlap","all")){
    .validateParams(pvalue=pvalue, minSize=minSize, maxSize=maxSize,
                    minGSSize=minGSSize, maxGSSize=maxGSSize, func_name="richKEGG")
    .check_cutoff(padj, "padj", allow_null = TRUE)
    normalization <- .normalize_annotation(kodata, func_name = "richKEGG")
  ## KEGG-specific: term name lookup + pathway level annotation
  if (ontology == "KEGGM") {
    all_ko <- .get_kgm.data()
    retrieved <- attr(all_ko, "richR_source", exact = TRUE)$retrieved
    label_provenance <- list(source = "bundled", release = NA_character_,
      retrieved = if (is.null(retrieved)) NA_character_ else as.character(retrieved))
  } else {
    ## The release is not queried here; it is recorded by buildAnnot.
    labels <- .kegg_labels(builtin = builtin)
    all_ko <- labels$data
    label_provenance <- labels[c("source", "release", "retrieved")]
  }
  term_names <- setNames(all_ko[, 1], rownames(all_ko))
  ## Labels: KEGG term list first, then the frame's Annot, then the identifier.
  covered <- names(term_names)[!is.na(term_names) & nzchar(trimws(term_names))]
  term_names <- .fill_labels(
    term_names, .frame_fallback_labels(normalization$data, covered, "richKEGG"))
  path_env <- new.env(parent = emptyenv())
  utils::data("path", package = "richR", envir = path_env)
    .run_ora(x, annot = normalization, term_names = term_names, extra_cols = path_env$path,
           pvalue = pvalue, padj = padj, padj.method = padj.method,
           minSize = minSize, maxSize = maxSize,
           minGSSize = minGSSize, maxGSSize = maxGSSize,
           keepRich = keepRich, organism = organism, ontology = ontology,
           keytype = keytype, filename = filename, sep = sep, universe = universe,
           annotation_provenance = annotation_provenance,
           metadata_ontology = ontology,
           padj.family = padj.family, label_provenance = label_provenance)
}
#' KEGG Pathway Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param ontology KEGG
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimum annotated term size for testing: the number of the term's genes in the background. \code{[minGSSize, maxGSSize]} (inclusive) bounds the tested family before p-value adjustment.
#' @param maxGSSize maximum annotated term size for testing (inclusive upper bound of the tested family; see \code{minGSSize}).
#' @param keepRich keep terms with rich factor value equal 1 even when their overlap is below \code{minSize} (default: TRUE). It exempts only the \code{minSize} reporting filter and never admits a term outside the tested family.
#' @param padj.family family of terms that \code{padj.method} adjusts over:
#'   \code{"overlap"} (default) holds the size-eligible terms (see \code{minGSSize})
#'   that share at least one gene with the query; \code{"all"} holds every
#'   size-eligible term, zero-overlap terms contributing p = 1 (they are never
#'   reported). \code{minSize}, \code{maxSize}, \code{keepRich} and the cutoffs are
#'   applied afterwards and never change an adjusted value. The default family
#'   matches clusterProfiler/DOSE but is a convention, not an FDR guarantee: with a
#'   universe of 100 single-gene terms and one selected gene under the null (with
#'   \code{minGSSize = 1, minSize = 1}), the only overlapping term has p = 0.01 and
#'   BH over that one term gives 0.01, so a null term is reported every time; BH
#'   over all 100 gives 1. Use \code{padj.family = "all"} for adjustment over every
#'   size-eligible term.
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want to use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#' }
#' @export
#' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("richKEGG", signature(kodata = "data.frame"),definition = function(x,kodata,pvalue=0.05,padj=0.05,organism=NULL,ontology="KEGG",
                                                                             keytype=NULL,minSize=2,maxSize=500,minGSSize = 10, maxGSSize = 500,
                                                                             keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",",universe=NULL,padj.family=c("overlap","all")) {
  richKEGG_internal(x,kodata=kodata,pvalue=pvalue,padj=padj,
                    organism=organism,ontology=ontology,keytype=keytype,minSize=minSize,
                    maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,builtin=builtin,sep=sep,universe=universe,padj.family=padj.family)
})

#' KEGG Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata Annot object with KEGG annotation
#' @param ontology KEGG
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimum annotated term size for testing: the number of the term's genes in the background. \code{[minGSSize, maxGSSize]} (inclusive) bounds the tested family before p-value adjustment.
#' @param maxGSSize maximum annotated term size for testing (inclusive upper bound of the tested family; see \code{minGSSize}).
#' @param keepRich keep terms with rich factor value equal 1 even when their overlap is below \code{minSize} (default: TRUE). It exempts only the \code{minSize} reporting filter and never admits a term outside the tested family.
#' @param padj.family family of terms that \code{padj.method} adjusts over:
#'   \code{"overlap"} (default) holds the size-eligible terms (see \code{minGSSize})
#'   that share at least one gene with the query; \code{"all"} holds every
#'   size-eligible term, zero-overlap terms contributing p = 1 (they are never
#'   reported). \code{minSize}, \code{maxSize}, \code{keepRich} and the cutoffs are
#'   applied afterwards and never change an adjusted value. The default family
#'   matches clusterProfiler/DOSE but is a convention, not an FDR guarantee: with a
#'   universe of 100 single-gene terms and one selected gene under the null (with
#'   \code{minGSSize = 1, minSize = 1}), the only overlapping term has p = 0.01 and
#'   BH over that one term gives 0.01, so a null term is reported every time; BH
#'   over all 100 gives 1. Use \code{padj.family = "all"} for adjustment over every
#'   size-eligible term.
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want to use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#' }
#' @export
#' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("richKEGG", signature(kodata = "Annot"),definition = function(x,kodata,pvalue=0.05,padj=0.05,organism=NULL,ontology="KEGG",
                                                                        keytype=NULL,minSize=2,maxSize=500,minGSSize = 10, maxGSSize = 500,
                                                                        keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",",universe=NULL,padj.family=c("overlap","all")) {
  provenance <- .annotation_provenance(kodata)
  richKEGG_internal(x,kodata@annot,pvalue=pvalue,padj=padj,
                  organism=kodata@species,ontology=kodata@anntype,keytype=kodata@keytype,minSize=minSize,
                  maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                  padj.method=padj.method,builtin=builtin,sep=sep,universe=universe,padj.family=padj.family,
                  annotation_provenance=provenance)
})



