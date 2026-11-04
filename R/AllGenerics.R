##' richGO
##'
#' @name richGO
#' @title GO Enrichment analysis function
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
#' @param ... additional arguments
#' @return richResult
##' @examples
#' \dontrun{
#' hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#' hsago <- as.data.frame(hsago)
#' gene <- sample(unique(hsago$GeneID),1000)
#' res <- richGO(gene,godata = hsago,ontology ="BP")
#' }
##' @export
##' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setGeneric("richGO", function(x,godata,ontology="BP",pvalue=0.05,padj=0.05,organism=NULL,keytype=NULL,minSize=2,maxSize=500,
                              minGSSize = 10, maxGSSize = 500,
                              keepRich=TRUE,filename=NULL,padj.method="BH",sep=",",universe=NULL,padj.family=c("overlap","all"),...)
  standardGeneric("richGO"))

##' richKEGG
#' @title KEGG Enrichment analysis function
#' @name richKEGG
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param pvalue cutoff pvalue
#' @param ontology KEGG
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
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @param ... additional arguments
##' @return richResult
#' @examples
#' \dontrun{
#' hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' gene <- sample(unique(hsako$GeneID),1000)
#' res <- richKEGG(gene,kodata = hsako)
#' }
##' @export
##' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setGeneric("richKEGG", function(x,kodata,pvalue=0.05,padj=0.05,organism=NULL,ontology="KEGG",
                                keytype=NULL,minSize=2,maxSize=500,minGSSize = 10, maxGSSize = 500,
                                keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",",universe=NULL,padj.family=c("overlap","all"),...)
  standardGeneric("richKEGG"))

##' richGSEA
##'
##' @name richGSEA
##' @rdname richGSEA-methods
##' @title richGSEA method
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
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
#' @param scoreType optional fgsea score type
#' @param nproc optional fgsea worker count
#' @param BPPARAM optional BiocParallel backend
#' @param seed optional seed scoped to this GSEA call
#' @param ... additional arguments
##' @return GSEAResult
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name <- sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene) <- name
#' res <- richGSEA(gene,object = hsako)
#' }
##' @export
##' @author Kai Guo
setGeneric("richGSEA", function(x,object,keytype="",pvalue=0.05,padj=0.05,KEGG=FALSE,minSize=15,ontology="",
                                maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",",
                                eps=NULL,scoreType=NULL,nproc=NULL,BPPARAM=NULL,seed=NULL,...)
  standardGeneric("richGSEA"))
##' enrich
##'
#' @name enrich
#' @title enrich method
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
#' @param ... additional arguments
##' @return richResult
##' @examples
##' \dontrun{
##' library(bioAnno)
##' fromKEGG(species="ath")
##' athgo <- buildOwn(dbname="org.ath.eg.db",anntype="GO")
##' athgo <- as.data.frame(athgo)
##' athgo <- na.omit(athgo)
##' gene <- sample(unique(athgo$GeneID),1000)
##' res <- enrich(gene,athgo)
##' }
##' @export
##' @author Kai Guo
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setGeneric("enrich", function(x,object,ontology="",pvalue=0.05,padj=0.05,organism=NULL,
                              keytype="",filename=NULL,minSize=2,maxSize=500,
                              minGSSize = 10, maxGSSize = 500,
                              keepRich=TRUE,padj.method="BH",sep=",",universe=NULL,padj.family=c("overlap","all"),...)
  standardGeneric("enrich"))


##' richLevel
##'
##' @name richLevel
##' @rdname richLevel-methods
##' @title richLevel method
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
#' @param ... additional arguments
#' @examples
#' \dontrun{
#'   hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene <- sample(unique(hsako$GeneID),1000)
#'   res <- richLevel(gene,kodata = hsako,level="Level2")
#' }
#' @export
#' @author Kai Guo
setGeneric("richLevel", function(x,kodata,level="Level2",pvalue =0.05, padj=0.05,
                                 organism=NULL,keytype="SYMBOL",ontology="",minSize=2,maxSize=500,
                                 minGSSize = 10, maxGSSize = 500,
                                 keepRich=TRUE, filename=NULL,padj.method="BH",sep=",",padj.family=c("overlap","all"),...)
  standardGeneric("richLevel"))

##' richDot
##'
##' @name richDot
##' @rdname richDot-methods
##' @title richDot method
##' @param object enrichment result or data.frame
##' @param top number of terms you want to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param order order by Term or richFactor
##' @param low low color
##' @param high high color
##' @param alpha transparency alpha
##' @param font.x font of x axis
##' @param font.y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param orderp order by p value (adjusted p value)
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param ... additional arguments
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene <- sample(unique(hsako$GeneID),1000)
#'   res <- richKEGG(gene,kodata = hsako)
#'   richDot(res)
#' }
##' @export
##' @author Kai Guo
setGeneric("richDot", function(object,top=50,pvalue=0.05,order=FALSE,
                             low="lightpink",high="red",alpha=0.7,
                             font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                             short=FALSE,
                             padj=0.05,usePadj=TRUE,orderp=FALSE,filename=NULL,width=10,height=8,...)
    standardGeneric("richDot"))

##' richBar
##'
##' @name richBar
##' @rdname richBar-methods
##' @title ggbar method
##' @param object enrichment result or data.frame
##' @param top number of terms you want to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param order order by Term or richFactor
##' @param usePadj use p adjust value as color or not
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param fontsize.text font size for bar text labels
##' @param angle angle of x axis text
##' @param orderp order by p value (adjusted p value)
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param horiz horizontal bar plot or not
##' @param ... additional arguments
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   richBar(res)
#' }
##' @export
##' @author Kai Guo
setGeneric("richBar", function(object,top=50,pvalue=0.05,padj=0.05,order=FALSE,
                             usePadj=TRUE,fontsize.x=10,fontsize.y=10,short=FALSE,fontsize.text=3,angle=0,orderp=FALSE,
                             filename=NULL,
                             width=10,height=8,horiz=TRUE,...)
    standardGeneric("richBar"))
##' richNetplot
##' @rdname richNetplot-method
##' @title ggnetplot method
##' @param object richResult or dataframe
##' @param top number of terms to show (default: 50)
##' @param pvalue cutoff p value for enrichment result
##' @param padj cutoff p adjust value for enrichment result
##' @param usePadj use adjust p value for the color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param useTerm use terms for nodes (default: TRUE)
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param label.color label color
##' @param label.size label size
##' @param node.shape vector of shape and names of the vector should be the terms (default: 20)
##' @param layout layout method
##' @param savefig save figures or not
##' @param filename output figure name
##' @param width width for output figure
##' @param height height for output figure
##' @param node.alpha alpha-transparency scales
##' @param node.shape shape of the node
##' @param repel use ggrepel text function or not
##' @param segment.size segment size for ggrepel text
##' @param sep character string used to separate the genes when concatenating
##' @param ... additional arguments
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene<-sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   richNetplot(res)
#' }
##' @export
##' @author Kai Guo
setGeneric("richNetplot",function(object,top=50, pvalue=0.05, padj=0.05,
                                usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                                writeCyt=FALSE, cytoscapeFile=NULL,
                                label.color = "black", label.size = 2, node.shape=NULL,
                                layout = "fruchtermanreingold",savefig=FALSE,filename=NULL,
                                width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2,sep=",",...)
  standardGeneric("richNetplot"))


##' richNetwork
##'
##' @name richNetwork
##' @rdname richNetwork-methods
##' @title ggnetwork method
##' @param object richResult,GSEAResult object or dataframe
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param low color used for small value
##' @param high color used for large value
##' @param weightcut cutoff value for edge weight
##' @param useTerm use Term name instead of Annotation ID
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat output file format
##' @param label.color label color
##' @param label.size label size
##' @param node.shape shape of the nodes
##' @param layout layout method
##' @param savefig save the figure or not
##' @param visNet use visNetwork for interactive plot
##' @param smooth smooth edges in visNetwork
##' @param nodeselect enable node selection in visNetwork
##' @param edit enable editing in visNetwork
##' @param savehtml save visNetwork as HTML file
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param sep character string used to separate the genes when concatenating
##' @param ... additional arguments
##' @return ggplot2 object
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene<-sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   richNetwork(res)
#' }
##' @export
setGeneric("richNetwork", function(object,gene,top = 50, pvalue = 0.05, padj = 0.05,usePadj=TRUE,low = "orange",high = "red",
                                 weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                                 label.color = "black", label.size = 2,node.shape=NULL, layout = "fruchtermanreingold",savefig=FALSE,
                                 visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename=NULL,
                                 width=7,height=7,segment.size=0.2,node.alpha=0.7,sep=",",...)
  standardGeneric("richNetwork"))


##' result generic
##' @param x richResult object
##' @return result return dataframe and print summary
##' @export
result<-function(x){
  UseMethod("result",x)
}
##' detail generic
##' @param x richResult object
##' @return detail return detail for these significant genes
##' @export
detail<-function(x){
  UseMethod("detail",x)
}
#' kappa cluster analysis
##' @param x richResult object or dataframe
##' @param gene (Optional).a vector of gene list
##' @param useTerm to use the term or not (TRUE/FALSE)
##' @param cutoff kappa score threshold for significant dispersion results
##' @param overlap cutoff value of the overlap between two Terms
##' @param minSize minimal number of terms in the cluster
##' @param escore kappa enrichment score cutoff value (default: 3)
#' @examples
#' \dontrun{
#'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   gene<-ample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#'   resc<-richCluster(res)
#' }
#' @export
#' @author Kai Guo
setGeneric("richCluster", function(x,gene=NULL,useTerm=FALSE,cutoff=0.5,overlap=0.5,minSize=5,escore=3)
  standardGeneric("richCluster"))


if (!isGeneric("head")) {
  setGeneric("head", function(x, n = 6L, ...) standardGeneric("head"))
}

if (!isGeneric("tail")) {
  setGeneric("tail", function(x, n = 6L, ...) standardGeneric("tail"))
}

if (!isGeneric("dim")) {
  setGeneric("dim", function(x) standardGeneric("dim"))
}
if(!isGeneric("summary")){
  setGeneric("summary", function(object, ...) standardGeneric("summary"))
}

#' S4 generic for extracting the result slot
#'
#' @param x An S4 object (e.g., \code{richResult} or \code{GSEAResult})
#' @return Typically a data frame (implementation-specific)
#' @export
#' @docType methods
#' @name result
#'
#' @examples
#' \dontrun{
#'   # usage: result(x)
#' }
NULL

if (!isGeneric("result")) {
  setGeneric("result", function(x) standardGeneric("result"))
}




