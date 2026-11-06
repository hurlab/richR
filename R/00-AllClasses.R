##' Class "richResult"
##' This class represents the result of enrichment analysis.
##'
##'
##' @name richResult-class
##' @aliases richResult-class plot,richResult-method
##'
##' @docType class
##' @slot result enrichment analysis results
##' @slot detail genes included in significant terms and original information
##' @slot tested tested table: one row per tested term whose overlap with the
##'   query is greater than zero, with the columns of \code{result} and a
##'   logical \code{reported} column (\code{TRUE} for the terms \code{result}
##'   holds); a zero-column data.frame when no tested table was recorded. Read
##'   it with \code{\link{testedTerms}}.
##' @slot pvalueCutoff cutoff pvalue
##' @slot pAdjustMethod pvalue adjust method
##' @slot padjCutoff pvalue adjust cutoff value
##' @slot genenumber number of input genes
##' @slot organism organism used
##' @slot ontology biological ontology
##' @slot gene Gene IDs
##' @slot keytype Gene ID type
##' @slot sep character string used to separate the genes when concatenating
##' @exportClass richResult
##' @author Kai Guo
##' @keywords classes
setClass("richResult",
         representation=representation(
           result         = "data.frame",
           detail         = "data.frame",
           tested         = "data.frame",
           pvalueCutoff   = "numeric",
           pAdjustMethod  = "character",
           padjCutoff   = "numeric",
           genenumber    = "numeric",
           organism       = "character",
           ontology       = "character",
           gene           = "character",
           keytype        = "character",
           sep = "character"
         ),
         prototype = prototype(tested = data.frame())
)
##' Class "GSEAResult"
##' This class represents the result from GSEA analysis
##' @name GSEAResult-class
##' @aliases GSEAResult-class plot,GSEAResult-method
##'
##' @docType class
##' @slot result enrichment analysis results
##' @slot tested tested table: one row per pathway fgsea returned, with the
##'   columns of \code{result} and a logical \code{reported} column
##'   (\code{TRUE} for the pathways \code{result} holds); a zero-column
##'   data.frame when no tested table was recorded. Read it with
##'   \code{\link{testedTerms}}.
##' @slot input named numeric vector of input gene scores
##' @slot pvalueCutoff cutoff pvalue
##' @slot pAdjustMethod pvalue adjust method
##' @slot padjCutoff pvalue adjust cutoff value
##' @slot genenumber number of input genes
##' @slot organism organism used
##' @slot ontology biological ontology
##' @slot gene Gene IDs
##' @slot keytype Gene ID type
##' @slot sep character string used to separate the genes when concatenating
##' @exportClass GSEAResult
##' @author Kai Guo
##' @keywords classes
setClass("GSEAResult",
         representation   = representation(
           result          = "data.frame",
           tested          = "data.frame",
           pvalueCutoff   = "numeric",
           pAdjustMethod  = "character",
           padjCutoff   = "numeric",
           genenumber    = "numeric",
           organism       = "character",
           gene           = "character",
           input = 'numeric',
           ontology = "character",
           keytype        = "character",
           sep = "character"
         ),
         prototype = prototype(tested = data.frame())
)
##' Class "Annot"
##' This class represents the Annotation information
##' @name Annot-class
##' @aliases Annot-class
##' @docType class
##' @slot species the species of the annotation file
##' @slot anntype the type of the annotation file
##' @slot keytype Gene ID type
##' @slot annot Annotation information data.frame
##' @exportClass Annot
##' @author Kai Guo
##' @keywords classes
setClass("Annot",
         representation = representation(
           species="character",
           anntype="character",
           keytype="character",
           annot="data.frame"
         ))

