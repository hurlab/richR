##' @method as.data.frame Annot
##' @export
as.data.frame.Annot<-function(x,...){
  as.data.frame(x@annot)
}
##' @method as.data.frame richResult
##' @export
as.data.frame.richResult <- function(x, ...) {
  as.data.frame(x@result, ...)
}
##' @method as.data.frame GSEAResult
##' @export
as.data.frame.GSEAResult <- function(x, ...) {
  as.data.frame(x@result, ...)
}

# S4 row- and column-related methods (docs in accessors.R)

# --- row.names() methods ---

#' @noRd
setMethod(
  "row.names",
  signature(x = "Annot"),
  function(x) {
    row.names(x@annot)
  }
)

#' @noRd
setMethod(
  "row.names",
  signature(x = "richResult"),
  function(x) {
    row.names(x@result)
  }
)

#' @noRd
setMethod(
  "row.names",
  signature(x = "GSEAResult"),
  function(x) {
    # if GSEAResult extends richResult
    row.names(x@result)
  }
)

#' @noRd
setMethod(
  "rownames",
  signature(x = "Annot"),
  function(x) {
    rownames(x@annot)
  }
)

#' @noRd
setMethod(
  "rownames",
  signature(x = "richResult"),
  function(x) {
    rownames(x@result)
  }
)

#' @noRd
setMethod(
  "rownames",
  signature(x = "GSEAResult"),
  function(x) {
    rownames(x@result)
  }
)

# --- names() methods ---

#' @noRd
setMethod(
  "names",
  signature(x = "Annot"),
  function(x) {
    names(x@annot)
  }
)

#' @noRd
setMethod(
  "names",
  signature(x = "richResult"),
  function(x) {
    names(x@result)
  }
)

#' @noRd
setMethod(
  "names",
  signature(x = "GSEAResult"),
  function(x) {
    names(x@result)
  }
)

# --- colnames() methods ---

#' @noRd
setMethod(
  "colnames",
  signature(x = "Annot"),
  function(x) {
    colnames(x@annot)
  }
)

#' @noRd
setMethod(
  "colnames",
  signature(x = "richResult"),
  function(x) {
    colnames(x@result)
  }
)

#' @noRd
setMethod(
  "colnames",
  signature(x = "GSEAResult"),
  function(x) {
    colnames(x@result)
  }
)

#' S4 methods for head, tail, dim on Annot, richResult, GSEAResult
#'
#' Defines how to call \code{head}, \code{tail}, and \code{dim} for these S4 classes.
#'
#' @param x An object of class \code{Annot}, \code{richResult}, or \code{GSEAResult}
#' @param n Number of rows to display (for \code{head} or \code{tail})
#' @param ... Further arguments passed to the underlying data.frame method
#'
#' @name utilities-methods
#' @rdname utilities-methods
NULL

# ---------------------------
#  HEAD METHODS
# ---------------------------
#' @rdname utilities-methods
#' @aliases head,Annot-method
setMethod(
  "head",
  signature = signature(x = "Annot"),
  function(x, n = 6L, ...) {
    cat(
      "=== species is:", x@species,
      "and Annotation is", x@anntype,
      "keytype is", x@keytype, "===\n"
    )
    utils::head(x@annot, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases head,richResult-method
setMethod(
  "head",
  signature = signature(x = "richResult"),
  function(x, n = 6L, ...) {
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::head(x@result, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases head,GSEAResult-method
setMethod(
  "head",
  signature = signature(x = "GSEAResult"),
  function(x, n = 6L, ...) {
    # GSEAResult extends richResult
    # We can either replicate the code, or call the richResult method
    # Example: call the same code
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::head(x@result, n, ...)
  }
)

# ---------------------------
#  TAIL METHODS
# ---------------------------
#' @rdname utilities-methods
#' @aliases tail,Annot-method
setMethod(
  "tail",
  signature = signature(x = "Annot"),
  function(x, n = 6L, ...) {
    cat(
      "=== species is:", x@species,
      "and Annotation is", x@anntype,
      "keytype is", x@keytype, "===\n"
    )
    utils::tail(x@annot, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases tail,richResult-method
setMethod(
  "tail",
  signature = signature(x = "richResult"),
  function(x, n = 6L, ...) {
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::tail(x@result, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases tail,GSEAResult-method
setMethod(
  "tail",
  signature = signature(x = "GSEAResult"),
  function(x, n = 6L, ...) {
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::tail(x@result, n, ...)
  }
)

# ---------------------------
#  DIM METHODS
# ---------------------------
#' @rdname utilities-methods
#' @aliases dim,Annot-method
setMethod(
  "dim",
  signature = signature(x = "Annot"),
  function(x) {
    dim(x@annot)
  }
)

#' @rdname utilities-methods
#' @aliases dim,richResult-method
setMethod(
  "dim",
  signature = signature(x = "richResult"),
  function(x) {
    dim(x@result)
  }
)

#' @rdname utilities-methods
#' @aliases dim,GSEAResult-method
setMethod(
  "dim",
  signature = signature(x = "GSEAResult"),
  function(x) {
    dim(x@result)
  }
)

## summary methods defined in accessors.R (canonical location)
## These are kept here only as fallback definitions without roxygen tags
setMethod(
  f = "summary",
  signature = "richResult",
  definition = function(object, ...) {
    cat("Total input genes is:", length(object@gene),
        "and significant biological term is:", nrow(object@result), "\n")
  }
)

setMethod(
  f = "summary",
  signature = "GSEAResult",
  definition = function(object, ...) {
    cat("Total significant biological term is:",
        sum(object@result$padj < 0.05, na.rm = TRUE), "\n")
  }
)

#' S4 subsetting methods for Annot, richResult, and GSEAResult
#'
#' These methods define `[` subsetting on each respective S4 class.
#'
#' @param x An S4 object (Annot, richResult, or GSEAResult)
#' @param i,j indices for subsetting
#' @param ... further arguments (usually ignored)
#' @param drop logical. Should dimensions be dropped?
#'
#' @name subset-methods
#' @rdname subset-methods
NULL

#' @rdname subset-methods
#' @aliases [,Annot-method
setMethod(
  f = "[",
  signature = signature(x = "Annot", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    x@annot[i, j, drop = drop]
  }
)

#' @rdname subset-methods
#' @aliases [,richResult-method
setMethod(
  f = "[",
  signature = signature(x = "richResult", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    x@result[i, j, drop = drop]
  }
)

#' @rdname subset-methods
#' @aliases [,GSEAResult-method
setMethod(
  f = "[",
  signature = signature(x = "GSEAResult", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    # if GSEAResult extends richResult, you might simply do:
    # callNextMethod() or x@result[i, j, drop=drop]
    x@result[i, j, drop = drop]
  }
)

#' S4 `$` accessors for Annot and richResult
#'
#' Provide `obj$name` access that extracts from the relevant slot columns.
#'
#' @param x An S4 object (Annot or richResult)
#' @param name column name to extract
#'
#' @name dollar-methods
#' @rdname dollar-methods
NULL

#' @rdname dollar-methods
#' @aliases $-Annot
setMethod(
  f = "$",
  signature = signature(x = "Annot"),
  definition = function(x, name) {
    # x@annot is presumably a data.frame or matrix, so
    # use standard R subsetting by column name
    x@annot[[name]]
  }
)

#' @rdname dollar-methods
#' @aliases $-richResult
setMethod(
  f = "$",
  signature = signature(x = "richResult"),
  definition = function(x, name) {
    x@result[[name]]
  }
)

#' @rdname dollar-methods
#' @aliases $-GSEAResult
setMethod(
  f = "$",
  signature = signature(x = "GSEAResult"),
  definition = function(x, name) {
    # If GSEAResult inherits from richResult
    x@result[[name]]
  }
)

#' S4 methods for result() and detail() on richResult and GSEAResult
#'
#' These methods return \code{as.data.frame(x@result)} or \code{x@detail},
#' depending on the class, effectively exposing object contents.
#'
#' @param x A \code{richResult} or \code{GSEAResult} object
#' @return A data frame
#'
#' @name result_detail_methods
#' @rdname result_detail_methods
NULL

#' @rdname result_detail_methods
#' @aliases result,richResult-method
#' @export
setMethod(
  f = "result",
  signature = "richResult",
  definition = function(x) {
    as.data.frame(x@result)
  }
)

#' @rdname result_detail_methods
#' @aliases result,GSEAResult-method
#' @export
setMethod(
  f = "result",
  signature = "GSEAResult",
  definition = function(x) {
    as.data.frame(x@result)
  }
)

#' @rdname result_detail_methods
#' @aliases detail,richResult-method
#' @export
setMethod(
  f = "detail",
  signature = "richResult",
  definition = function(x) {
    as.data.frame(x@detail)
  }
)

#' @rdname result_detail_methods
#' @aliases detail,GSEAResult-method
#' @export
setMethod(
  f = "detail",
  signature = "GSEAResult",
  definition = function(x) {
    as.data.frame(x@detail)
  }
)

#' reverse List
#' @param lhs list with names
#' @export
#' @author Kai Guo
# replace this with c++
reverseList_bk<-function(lhs){
  lhs_n<-rep(names(lhs),times=lapply(lhs,function(x)length(x)))
  res<-sf(as.data.frame(cbind(lhs_n,unlist(lhs))))
  return(res)
}

##'
setAs(from = "data.frame", to = "Annot", def = function(from){
  keytype <- character()
  species <- character()
  anntype <- character()
  GeneID <- as.vector(from[,1])
  Term<-as.vector(from[,2])
  Annot=from$Annot
  annot <- data.frame(GeneID, Term, Annot)
  new("Annot",
      species = species,
      anntype = anntype,
      keytype = keytype,
      annot = annot
  )
})

##'
##'
setAs(from = "data.frame", to = "richResult", def = function(from){
  keytype <- character()
  organism <- character()
  ontology <- character()
  pvalueCutoff <- numeric()
  pAdjustMethod <-character()
  padjCutoff <- numeric()
  Annot <- from$Annot
  Term <- from$Term
  Annotated <- from$Annotated
  Significant <- from$Significant
  Pvalue <- from$Pvalue
  Padj <- from$Padj
  GeneID <- as.vector(from$GeneID)
  gene<-unique(unlist(strsplit(GeneID,",")))
  genenumber <- length(gene)
  resultFis <- data.frame(Annot, Term, Annotated, Significant, Pvalue, Padj, GeneID)
  rownames(resultFis) <- Annot
  new("richResult",
      result=resultFis,
      detail=data.frame(),
      pvalueCutoff   = pvalueCutoff,
      pAdjustMethod  = pAdjustMethod,
      padjCutoff   = padjCutoff,
      genenumber    = genenumber,
      organism       = organism,
      ontology       = ontology,
      gene           = gene,
      keytype        = keytype,
      sep            = ","
  )
})

#' rbind generic function for richResult object
#' @param ... richResult objects to combine
#'@export
#'@author Kai Guo
rbind.richResult<-function(...){
    objects <- list(...)
    objects <- lapply(objects,as.data.frame)
    do.call(rbind, objects)
}

#' rbind generic function for GSEAResult object
#' @param ... GSEAResult objects to combine
#'@export
#'@author Kai Guo
rbind.GSEAResult<-function(...){
  objects <- list(...)
  objects <- lapply(objects,as.data.frame)
  do.call(rbind, objects)
}

#' Remove newlines
#'
#' Replaces any newline characters with a space.
#'
#' @param x Character vector.
#' @return Character vector without newlines.
#' @keywords internal
.clean.char <- function(x) {
  gsub("\\n", " ", x)
}

##'
setAs(from = "richResult", to = "data.frame", def = function(from){
  result <- as.data.frame(from@result)
  result
})
##'
setAs(from = "GSEAResult", to = "data.frame", def = function(from){
  result <- as.data.frame(from@result)
  result
})

##'
setAs(from = "Annot", to = "data.frame", def = function(from){
  result <- as.data.frame(from@annot)
  result
})

##' kappa function
##' @param x comma-separated gene string for first term
##' @param y comma-separated gene string for second term
##' @param geneall character vector of all genes (background)
.kappa<-function(x,y,geneall){
  x<-unlist(strsplit(x,","))
  y<-unlist(strsplit(y,","))
  if(length(intersect(x,y))==0){
    kab=0
  }else{
    tmp<-matrix(0,2,2)
    tmp[1,1]<-length(intersect(x,y))
    tmp[2,1]<-length(setdiff(x,y))
    tmp[1,2]<-length(setdiff(y,x))
    tmp[2,2]<-length(setdiff(geneall,union(x,y)))
    oab<-(tmp[1,1]+tmp[2,2])/sum(tmp)
    aab<-((tmp[1,1]+tmp[2,1])*(tmp[1,1]+tmp[1,2])+(tmp[1,2]+tmp[2,2])*(tmp[2,1]+tmp[2,2]))/(sum(tmp)*sum(tmp))
    if(aab==1){
      kab=0
    }else{
      kab<-(oab-aab)/(1-aab)
    }
  }
  return(kab)
}
##' calculate enrichment score
##' @param x character vector of term identifiers
##' @param df data.frame with enrichment results containing Pvalue column
.calculate_Enrichment_Score<-function(x,df){
  pvalue <- df[x,"Pvalue"]
  esp = ifelse(pvalue==0,16,-log10(pvalue))
  es = sum(esp);
}

##' merge term
##' @param x named list of term clusters
##' @param overlap numeric overlap threshold for merging
.merge_term<-function(x,overlap){
  ml <- x
  res<-list();
  for(i in names(ml)){
    lhs <- setdiff(names(ml),i)
    for(j in lhs){
      ov<-intersect(ml[[i]],ml[[j]])
      un<-union(ml[[i]],ml[[j]])
      ovl<-length(ov)/length(un)
      if(ovl > overlap){
        res[[i]]<-c(i,un)
        ml <- ml[setdiff(names(ml),j)]
      }else{
        res[[i]]<-c(i,ml[[i]])
      }
    }
  }
  return(res)
}

###
.calGSEA<-function (obj,sigpathway, fc, gseaParam = 1, ticksSize = 0.2){
  obj<-split(obj[,1],obj[,ncol(obj)])
  pathway <- obj[[sigpathway]]
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

