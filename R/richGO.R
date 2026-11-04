#' GO Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param godata GO annotation data
#' @param ontology BP,MF or CC
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimum annotated term size for testing: the number of the term's genes in the background. \code{[minGSSize, maxGSSize]} (inclusive) bounds the tested family before p-value adjustment.
#' @param maxGSSize maximum annotated term size for testing (inclusive upper bound of the tested family; see \code{minGSSize}).
#' @param keepRich keep terms with rich factor value equal 1 even when their overlap is below \code{minSize} (default: FALSE). It exempts only the \code{minSize} reporting filter and never admits a term outside the tested family.
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
#'   The annotation is narrowed to \code{ontology} before the background is built, so the tested family (and the default background) is ontology-specific.
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @param universe optional character vector defining the ORA background
#' @param annotation_provenance optional annotation provenance captured from an Annot object
#' @author Kai Guo
richGO_internal<-function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,
                 organism=NULL,keytype="SYMBOL",minSize=2,maxSize=500,
                 minGSSize = 10, maxGSSize = 500,
                 keepRich=FALSE, filename=NULL,padj.method="BH",sep=",",universe=NULL,
                 annotation_provenance=NULL,padj.family=c("overlap","all")){
  .validateParams(pvalue=pvalue, padj=padj, minSize=minSize, maxSize=maxSize,
                  minGSSize=minGSSize, maxGSSize=maxGSSize, func_name="richGO")
  normalization <- .normalize_annotation(godata, func_name = "richGO")
  godata <- normalization$data
  ## GO-specific: filter annotation to terms in this ontology
  if ("ONTOLOGYALL" %in% colnames(godata)) {
    godata_filt <- godata[which(godata$ONTOLOGYALL == ontology), ]
  } else {
    all_go <- .get_go_dat(ont = ontology)
    godata_filt <- godata[godata[, 2] %in% rownames(all_go), ]
  }
  if (nrow(godata_filt) == 0) {
    stop("richGO: no GO terms found for ontology '", ontology,
         "'. Check that your annotation contains this ontology.", call. = FALSE)
  }
  ## term_names: named vector  GO_ID -> readable term
  if ("Annot" %in% colnames(godata_filt) && "GOALL" %in% colnames(godata_filt)) {
    ## Build from existing annotation columns (avoids requiring GO.db)
    tm <- unique(godata_filt[, c("GOALL", "Annot")])
    term_names <- setNames(tm$Annot, tm$GOALL)
    annot_2col <- godata_filt[, c("GeneID", "GOALL"), drop = FALSE]
  } else {
    all_go <- .get_go_dat(ont = ontology)
    term_names <- setNames(all_go[, 1], rownames(all_go))
    annot_2col <- godata_filt[, 1:2, drop = FALSE]
  }
  normalization$counts$eligible_rows <- as.integer(nrow(annot_2col))
  normalization$data <- annot_2col
  .run_ora(x, annot = normalization, term_names = term_names,
           pvalue = pvalue, padj = padj, padj.method = padj.method,
           minSize = minSize, maxSize = maxSize,
           minGSSize = minGSSize, maxGSSize = maxGSSize,
           keepRich = keepRich, organism = organism, ontology = ontology,
           keytype = keytype, filename = filename, sep = sep, universe = universe,
           annotation_provenance = annotation_provenance,
           metadata_ontology = ontology,
           padj.family = padj.family)
}
#' GO Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param godata GO annotation data
#' @param ontology BP,MF or CC
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
#'   The annotation is narrowed to \code{ontology} before the background is built, so the tested family (and the default background) is ontology-specific.
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   hsago <- as.data.frame(hsago)
#'   gene <- sample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#' }
#' @export
#' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("richGO", signature(godata = "data.frame"),definition = function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,
                                                                           organism=NULL,keytype=NULL,minSize=2,maxSize=500,
                                                                           minGSSize = 10, maxGSSize = 500,
                                                                           keepRich=TRUE, filename=NULL,padj.method="BH",sep=",",universe=NULL,padj.family=c("overlap","all")) {
 # godata<-as(godata,"Annot")
  richGO_internal(x,godata,ontology=ontology,pvalue=pvalue,padj=padj,
                  organism=organism,keytype=keytype,minSize=minSize,maxSize=maxSize,
                  minGSSize = minGSSize, maxGSSize = maxGSSize,
                  keepRich=keepRich,filename=filename,padj.method=padj.method,sep=sep,universe=universe,padj.family=padj.family)
})

#' GO Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param godata GO annotation data
#' @param ontology BP,MF or CC
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
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
#'   The annotation is narrowed to \code{ontology} before the background is built, so the tested family (and the default background) is ontology-specific.
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   gene=sample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#' }
#' @export
#' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("richGO", signature(godata = "Annot"),definition = function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,minSize=2,maxSize=500,
                                                                      minGSSize = 10, maxGSSize = 500,
                                                                      keepRich=TRUE,filename=NULL,padj.method="BH",sep=",",universe=NULL,padj.family=c("overlap","all")) {
  provenance <- .annotation_provenance(godata)
  richGO_internal(x,godata@annot,ontology=ontology,pvalue=pvalue,padj=padj,
                  organism=godata@species,keytype=godata@keytype,minSize=minSize,maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,
                  keepRich=keepRich,filename=filename,padj.method=padj.method,sep=sep,universe=universe,padj.family=padj.family,
                  annotation_provenance=provenance)
})
