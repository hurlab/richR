#' Plot compare heatmap of Enrichment result among DEG groups
#' @importFrom dplyr full_join
#' @importFrom dplyr arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 coord_equal
#' @importFrom ggplot2 coord_flip
#' @param richRes list of enrichment object
#' @param top the number of Terms you want to display
#' @param colnames the compare DEG group names
#' @param xsize size of group name
#' @param ysize size of Terms name
#' @param usePadj use adjusted p-value instead of p-value
#' @param horizontal flip the heatmap horizontally
#' @param returnData return the data instead of the plot
#' @param ... additional arguments
#' @section Row selection:
#' Draws the union of the first \code{top} identifiers of each result (rows
#' are matched across results by \code{Annot}, else \code{Term}); no
#' significance predicate is applied. When this drops identifiers it signals
#' one \code{\link{richR_plot_selection}} message counted in identifiers over
#' all results (\code{suppressMessages()} silences it). To draw every
#' identifier: \code{richHeatmap(x, top = Inf)}.
#' @section Clustering:
#' Correlation clustering retains the original term order when fewer than
#' two terms or groups are available or any distance is undefined. Raw data
#' requested with \code{returnData = TRUE} does not require clustering.
#' @section Absent terms:
#' An identifier absent from one result has no p-value there: its cell is
#' \code{NA} in the drawn data and in the \code{returnData = TRUE} table and is
#' drawn in the missing-value colour. For the row ordering only it counts as
#' \code{p = 1}, so the order is as before. The fill legend names the column
#' drawn: \code{-log10(Pvalue)} with \code{usePadj = FALSE},
#' \code{-log10(Padj)} with \code{usePadj = TRUE}. Before richR 0.2.0 such a
#' cell held \code{1} and the legend named the other column.
#' @examples
#' \dontrun{
#'   hsago <- buildAnnot(species="human", keytype="SYMBOL", anntype="GO")
#'   hsako <- buildAnnot(species="human", keytype="SYMBOL", anntype="KEGG")
#'   gene <- sample(unique(hsago$GeneID), 1000)
#'   resgo <- richGO(gene, godata = hsago, ontology = "BP")
#'   resko <- richKEGG(gene, kodata = hsako)
#'   richHeatmap(list(GO = resgo, KEGG = resko), top = 50)
#' }
#' @export
#' @author Kai Guo
richHeatmap<-function(richRes, top = 50, colnames = NULL, xsize = 6, ysize = 6,usePadj=FALSE,
                     horizontal=FALSE,returnData=FALSE,...)
{
  # Rows are matched across results by identifier (`Annot`); a legacy result
  # without `Annot` is matched by `Term`.
  key <- if (all(vapply(richRes, function(x) "Annot" %in% colnames(x@result),
                        logical(1)))) "Annot" else "Term"
  sel <- unique(unlist(lapply(richRes, function(x) {
    r <- x@result
    head(as.character(r[[key]]), min(top, nrow(r)))
  })))
  n_ids <- length(unique(unlist(lapply(richRes, function(x) as.character(x@result[[key]])))))
  .plot_selection("richHeatmap", length(sel), n_ids, unit = "identifiers",
                  rule = paste("first", format(top), "identifiers of each result; union"))
  value_col <- if (isTRUE(usePadj)) "Padj" else "Pvalue"
  res <- Reduce(function(x, y) full_join(x, y, by = key),
                lapply(richRes, function(x) x@result[, c(key, value_col)]))
  if (!is.null(colnames)) {
    colnames(res)[2:ncol(res)] <- colnames
  }
  else {
    colnames(res)[2:ncol(res)] <- paste("Group", 1:(ncol(res)-1), sep = "_")
  }
  # An identifier absent from a result stays NA in the drawn data and in the
  # returned table: it has no p-value there (REQ-005). For the row ordering
  # only (ord_values below) it counts as p = 1, as before, so the order is
  # unchanged.
  res <- res[as.character(res[[key]]) %in% sel, , drop = FALSE]
  res<-as.data.frame(res)
  if (key == "Annot") {
    # Display labels over the union of results: one row per identifier.
    pool <- do.call(rbind, lapply(richRes, function(x) x@result[, c("Annot", "Term")]))
    res[[1]] <- .gsea_display_labels(res$Annot, pool$Annot, pool$Term)
    colnames(res)[1] <- "Term"
  }
  rownames(res)<-res$Term
  if (isTRUE(returnData)) return(res)
  ord_values <- res[, 2:ncol(res), drop = FALSE]
  ord_values[is.na(ord_values)] <- 1
  # Correlation clustering needs two terms and varying profiles over groups.
  # Undefined distances retain the input order; cell evidence stays untouched.
  row_order <- seq_len(nrow(res))
  if (nrow(ord_values) >= 2L && ncol(ord_values) >= 2L &&
      all(apply(ord_values, 1, function(x) is.finite(sd(x)) && sd(x) > 0))) {
    cor_mat <- cor(t(ord_values))
    dd <- as.dist((1-cor_mat)/2)
    if (all(is.finite(dd))) row_order <- hclust(dd)$order
  }
  res_ord <- res[row_order, , drop = FALSE]
  value_cols <- setdiff(colnames(res_ord), "Term")
  melted <- data.frame(
    Term = rep(res_ord$Term, times = length(value_cols)),
    variable = rep(value_cols, each = nrow(res_ord)),
    value = unlist(res_ord[, value_cols], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  melted$Term <- factor(melted$Term, levels = res_ord$Term)
  maxp = max(-log10(melted[, 3]), na.rm = TRUE) + 0.5
  # The fill is named after the column it draws (REQ-005).
  if(!isTRUE(usePadj)){
    colnames(melted)[3] <- "Pvalue"
    p<-ggplot(melted, aes(x = variable, y = Term, fill = -log10(Pvalue))) +coord_equal(ratio = 0.8)+
      geom_tile(color = "white") + scale_fill_gradient2(low = "white",high = "red", midpoint = 0, limit = c(0, maxp)) +
      theme_minimal()+theme(axis.text.y = element_text(size = ysize), axis.text.x = element_text(angle = 70,
      vjust = 1, size = xsize, hjust = 1,face = "bold")) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }else{
    colnames(melted)[3] <- "Padj"
    p<-ggplot(melted, aes(x = variable, y = Term, fill = -log10(Padj))) +coord_equal(ratio = 0.8)+
      geom_tile(color = "white") + scale_fill_gradient2(low = "white",high = "red", midpoint = 0, limit = c(0, maxp)) +
      theme_minimal()+theme(axis.text.y = element_text(size = ysize), axis.text.x = element_text(angle = 70,
      vjust = 1, size = xsize, hjust = 1,face = "bold")) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }
  if(isTRUE(horizontal)){
    p<-p+coord_flip()
  }
  if(isTRUE(returnData)){
    return(res)
  }else{
    p
  }
}
