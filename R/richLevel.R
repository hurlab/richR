#' Pathway Enrichment analysis for different level function
#' @importFrom dplyr filter left_join
#' @importFrom rlang sym
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param level pathway level ("Level1", "Level2", or "Level3")
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param ontology ontology type
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
#'   Terms are the entries of the requested \code{level}, so each call adjusts over that hierarchy level only (one family per level).
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @author Kai Guo
richLevel_internal<-function(x,kodata,level="Level2",pvalue =0.05, padj=0.05,
                             organism=NULL,keytype="SYMBOL",ontology="",minSize=2,maxSize=500,
                             minGSSize = 10, maxGSSize = 500,
                             keepRich=TRUE, filename=NULL,padj.method="BH",sep=",",padj.family=c("overlap","all")){
  path_env <- new.env(parent = emptyenv())
  utils::data("path", package = "richR", envir = path_env)
  annot<-left_join(kodata,path_env$path,by=c('PATH'="ko"))
  annot<-annot[,c("GeneID",level)]
  annot<-annot[!is.na(annot[,2]),]
  result<-enrich(x = x,object = annot,pvalue=pvalue,padj=padj,organism=organism,minSize=minSize,
              maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,keytype=keytype,ontology=ontology,filename=filename,
              padj.method=padj.method,sep = sep,padj.family=padj.family)
  return(result)
}
#' Pathway Enrichment analysis for different level function
#' @importFrom dplyr filter left_join
#' @importFrom rlang sym
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param level pathway level ("Level1", "Level2", or "Level3")
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param ontology ontology type
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
#'   Terms are the entries of the requested \code{level}, so each call adjusts over that hierarchy level only (one family per level).
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richLevel(gene,kodata = hsako,level="Level2")
#' }
#' @export
#' @author Kai Guo
setMethod("richLevel", signature(kodata = "data.frame"),definition = function(x,kodata,level="Level2",pvalue =0.05, padj=0.05,
                                                                              organism=NULL,keytype="SYMBOL",ontology="KEGG",minSize=2,maxSize=500,
                                                                              minGSSize = 10, maxGSSize = 500,
                                                                              keepRich=TRUE, filename=NULL,padj.method="BH",sep=",",padj.family=c("overlap","all")) {
  richLevel_internal(x,kodata=kodata,level=level, pvalue=pvalue,padj=padj,
                    organism=organism,keytype=keytype,ontology=ontology,minSize=minSize,
                    maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,sep=sep,padj.family=padj.family)
})

#' Pathway Enrichment analysis for different level function
#' @importFrom dplyr filter left_join
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param level pathway level ("Level1", "Level2", or "Level3")
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param ontology ontology type
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
#'   Terms are the entries of the requested \code{level}, so each call adjusts over that hierarchy level only (one family per level).
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richLevel(gene,kodata = hsako,level="Level2")
#' }
#' @export
#' @author Kai Guo
setMethod("richLevel", signature(kodata = "Annot"),definition = function(x,kodata,level="Level2",pvalue =0.05, padj=0.05,
                                                                         organism=NULL,keytype="SYMBOL",ontology="",minSize=2,maxSize=500,
                                                                         minGSSize = 10, maxGSSize = 500,
                                                                         keepRich=TRUE, filename=NULL,padj.method="BH",sep=",",padj.family=c("overlap","all")) {
  richLevel_internal(x,kodata@annot,level=level,pvalue=pvalue,padj=padj,
                    organism=kodata@species,keytype=kodata@keytype,ontology=kodata@anntype,minSize=minSize,
                    maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,sep=sep,padj.family=padj.family)
})
