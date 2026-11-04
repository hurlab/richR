
# S4 row- and column-related methods (docs in accessors.R)

# --- row.names() methods ---

# --- names() methods ---

# --- colnames() methods ---

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

