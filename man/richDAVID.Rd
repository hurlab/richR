% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.R
\name{richDAVID}
\alias{richDAVID}
\title{Functional enrichment analysis with DAVID}
\usage{
richDAVID(
  gene,
  keytype = "ENTREZ_GENE_ID",
  species = "human",
  anntype = "GOTERM_BP_FAT",
  universe,
  pvalue = 0.05,
  padj = NULL,
  minSize = 2,
  maxSize = 500,
  keepRich = TRUE,
  filename = NULL,
  padj.method = "BH",
  sep = ",",
  david.user = "richR@und.edu"
)
}
\arguments{
\item{gene}{vector contains gene names}

\item{keytype}{key type export}

\item{anntype}{GOTERM_BP_FAT, KEGG_PATHWAY,GOTERM_CC_FAT}

\item{pvalue}{cutoff pvalue}

\item{padj}{cutoff p adjust value}

\item{david.user}{richR@und.edu}
}
\value{
Annot object
}
\description{
Functional enrichment analysis with DAVID
}
\examples{
\dontrun{
hsako <- buildAnnot(species="human",keytype="ENTREZID",anntype = "KEGG")
hsako <- as.data.frame(hsako)
gene <- sample(unique(hsako$GeneID),1000)
res <- richDAVID(gene,keytype="ENTREZID",species="human")
}
}
\author{
Kai Guo
}
