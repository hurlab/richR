#' Enrichment analysis for any type of annotation data
#' @param x vector contains gene names or dataframe with DEGs information
#' @param object annotation data
#' @param ontology ontology type
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
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @param universe optional character vector defining the ORA background
#' @param annotation_provenance optional annotation provenance captured from an Annot object
#' @author Kai Guo
enrich_internal<-function(x,object,ontology= "",pvalue=0.05,padj=NULL,organism=NULL,minSize=1,maxSize=500,
                          minGSSize = 10, maxGSSize = 500,
                          keepRich=TRUE,keytype="",filename=NULL,padj.method="BH",sep = ",",universe=NULL,
                          annotation_provenance=NULL,padj.family=c("overlap","all")){
  .validateParams(pvalue=pvalue, padj=padj, minSize=minSize, maxSize=maxSize,
                  minGSSize=minGSSize, maxGSSize=maxGSSize, func_name="enrich")
  normalization <- .normalize_annotation(object, func_name = "enrich")
  .run_ora(x, annot = normalization, term_names = NULL, extra_cols = NULL,
           pvalue = pvalue, padj = padj, padj.method = padj.method,
           minSize = minSize, maxSize = maxSize,
           minGSSize = minGSSize, maxGSSize = maxGSSize,
           keepRich = keepRich, organism = organism, ontology = "",
           keytype = keytype, filename = filename, sep = sep, universe = universe,
           annotation_provenance = annotation_provenance,
           metadata_ontology = ontology,
           padj.family = padj.family)
}
#' Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param object annotation data
#' @param ontology ontology type
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
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   hsago <- as.data.frame(hsago)
#'   gene <- sample(unique(hsago$GeneID),1000)
#'   res<-enrich(gene,godata = hsago)
#' }
#' @export
#' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("enrich", signature(object = "data.frame"),definition = function(x,object,ontology="",pvalue=0.05,padj=NULL,organism=NULL,
                                                                             keytype="",filename=NULL,minSize=2,maxSize=500,
                                                                           minGSSize = 10, maxGSSize = 500,
                                                                             keepRich=TRUE,padj.method="BH",sep=",",universe=NULL,padj.family=c("overlap","all")) {
  enrich_internal(x,object=object,ontology=ontology,pvalue=pvalue,padj=padj,
                    organism=organism,keytype=keytype,minSize=minSize,maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,
                    filename=filename,padj.method=padj.method,sep=sep,universe=universe,padj.family=padj.family)
})

#' Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param object annotation data
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
#' @export
#' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("enrich", signature(object = "Annot"),definition = function(x,object,pvalue=0.05,padj=NULL,organism=NULL,
                                                                        keytype="",filename=NULL,minSize=2,maxSize=500,
                                                                      minGSSize = 10, maxGSSize = 500,
                                                                        keepRich=TRUE,padj.method="BH",sep=",",universe=NULL,builtin=TRUE,padj.family=c("overlap","all")) {
  provenance <- .annotation_provenance(object)
  enrich_internal(x=x,object=object@annot,ontology=object@anntype,pvalue=pvalue,padj=padj,
                    organism=object@species,keytype=object@keytype,minSize=minSize,maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,
                    filename=filename,padj.method=padj.method,sep=sep,universe=universe,padj.family=padj.family,
                    annotation_provenance=provenance)
})


#' Functional enrichment analysis with DAVID
#'
#' @param gene vector contains gene names
#' @param keytype key type export
#' @param species species name
#' @param anntype GOTERM_BP_FAT, KEGG_PATHWAY,GOTERM_CC_FAT
#' @param universe background gene universe
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method (default: "BH")
#' @param sep character string used to separate the genes when concatenating
#' @param david.user richR@und.edu
#' @return Annot object
#' @examples
#' \dontrun{
#' hsako <- buildAnnot(species="human",keytype="ENTREZID",anntype = "KEGG")
#' hsako <- as.data.frame(hsako)
#' gene <- sample(unique(hsako$GeneID),1000)
#' res <- richDAVID(gene,keytype="ENTREZID",species="human")
#' }
#' @export
#' @author Kai Guo
richDAVID <- function(gene,keytype="ENTREZ_GENE_ID",species="human",anntype="GOTERM_BP_FAT",
                      universe,pvalue=0.05,padj=NULL,minSize=2,maxSize=500,
                      keepRich=TRUE, filename=NULL,padj.method="BH",sep=",",
                      david.user="richR@und.edu"){
  pkg <- "RDAVIDWebService"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for DAVID analysis. Install it with:\n",
         "  BiocManager::install('", pkg, "')", call. = FALSE)
  }
  #require(pkg, character.only=TRUE)
  idtype=keytype
  if(keytype=="SYMBOL"){
    gene<-idconvert(species=species,keys=gene,fkeytype = "SYMBOL",tkeytype = "ENTREZID")
    keytype="ENTREZ_GENE_ID"
  }else if(keytype=="ENTREZID"){
    keytype="ENTREZ_GENE_ID"
  }else if(keytype=="ENSEMBL"){
    keytype="ENSEMBL_GENE_ID"
  }else if(keytype=="REFSEQ"){
    keytype="REFSEQ_MRNA"
  }else{
    keytype=keytype
  }
  david <- DAVIDWebService$new(email=david.user,
                               url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  keytype <- match.arg(keytype, getIdTypes(david))
  rh <- addList(david, gene, idType=keytype,
                listName="richR",
                listType="Gene")
  if (rh$inDavid == 0) {
    stop("Please check 'keytype' parameter...")
  }
  if (!missing(universe)) {
    rh <- addList(david, universe, idType=keytype,
                  listName="universe",
                  listType="Background")
  }
  #####
  setAnnotationCategories(david, anntype)
  res <- getFunctionalAnnotationChart(david, threshold=1, count=minSize)
  if (length(res@.Data) == 0) {
    warning("No significant enrichment found...")
    return(NULL)
  }
  terms <- res$Term
  if (length(grep("~", terms[1])) == 0) {
    sep <- ":"
  } else {
    sep <- "~"
  }
  terml <- sapply(terms, function(y) strsplit(y, split=sep))
  term <- do.call("rbind", terml)
  Annot <- term[,1]
  Term <- term[,2]
  Annotated <- res$Pop.Hits
  Significant <- res$Count
  RichFactor <- Significant / Annotated
  resultFis <- data.frame(Annot          = Annot,
                          Term = Term,
                          Annotated   = Annotated,
                          Significant     = Significant,
                          RichFactor = RichFactor,
                          Pvalue      = res$PValue,
                          stringsAsFactors = FALSE)
  row.names(resultFis) <- Annot
  if (padj.method == "bonferroni") {
    resultFis$Padj <- res$Bonferroni
  } else if(padj.method == "FDR") {
    resultFis$Padj <- res$FDR
  }else{
    resultFis$Padj <- res$Benjamini
  }
  resultFis$GeneID <- gsub(' ','',res$Genes)
  resultFis<-resultFis[order(resultFis$Pvalue),]
  if(is.null(padj)){
    resultFis<-resultFis[resultFis$Pvalue<pvalue,]
    padj=numeric()
  }else{
    resultFis<-resultFis[resultFis$Padj<padj,]
  }
  resultFis<-subset(resultFis, Significant<=maxSize)
  if(keepRich==FALSE){
    resultFis<-subset(resultFis, Significant>=minSize)
  }else{
    resultFis<-subset(resultFis, Significant>=minSize|(Annotated/Significant)==1)
  }
  rownames(resultFis)<-resultFis$Annot
  if(idtype=="SYMBOL"){
    xx<-suppressMessages(lapply(strsplit(resultFis$GeneID,","),function(x)idconvert(species=species,keys=x,fkeytype = "ENTREZID",tkeytype = "SYMBOL")))
    xx<-unlist(lapply(xx,function(x)paste(x,sep="",collapse = ",")))
    resultFis$GeneID <- xx
  }
  if(!is.null(filename)){
    write.table(resultFis,file=paste(filename,anntype,"res.txt",sep="_"),sep="\t",quote=F,row.names=F)
  }
  genes<-strsplit(as.vector(resultFis$GeneID),split=sep)
  names(genes)<-resultFis$Annot
  gened<-data.frame("TERM"=rep(names(genes),times=unlist(lapply(genes,length))),
                    "Annot"=rep(resultFis$Term,times=unlist(lapply(genes,length))),
                    "GeneID"=unlist(genes),row.names=NULL,
                    "Pvalue"=rep(resultFis$Pvalue,times=unlist(lapply(genes,length))),
                    "Padj"=rep(resultFis$Padj,times=unlist(lapply(genes,length)))
  )
  gened$GeneID<-as.character(gened$GeneID)
  detail<-gened
  species <- getSpecieNames(david)
  species <- gsub("\\(.*\\)", "", species)
  new("richResult",
      result=resultFis,
      detail=detail,
      pvalueCutoff   = pvalue,
      pAdjustMethod  = padj.method,
      padjCutoff   = padj,
      genenumber    = length(gene),
      organism       = species,
      ontology       = anntype,
      gene           = gene,
      keytype        = keytype,
      sep=",")
}

