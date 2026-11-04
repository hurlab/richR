#' GO term to offspring mapping
#'
#' A named list that maps each Gene Ontology term to itself and all of its
#' offspring terms, covering the biological process, molecular function and
#' cellular component ontologies.
#'
#' @format A named list of 42,984 character vectors. Each name is a GO
#'   identifier (for example \code{"GO:0000002"}) and each element holds that
#'   identifier followed by the identifiers of all of its offspring terms.
#' @details richR does not load this object: GO offspring terms are read from
#'   GO.db at run time. It is kept in the package so that user code that loads
#'   it keeps working. The date of the bundled snapshot is unknown, and no
#'   script rebuilds it.
#' @source Gene Ontology offspring relations. The exact source and GO release
#'   are unknown; the mapping has the shape of the \code{GOBPOFFSPRING},
#'   \code{GOCCOFFSPRING} and \code{GOMFOFFSPRING} tables of the Bioconductor
#'   GO.db package.
#' @usage data(godata)
"godata"

#' KEGG pathway database
#'
#' KEGG pathway identifiers and names, used to label pathway terms in
#' \code{\link{richKEGG}}.
#'
#' @format A data.frame with 567 rows and one column. The row names are KEGG
#'   pathway identifiers without the organism prefix (for example
#'   \code{"01100"}); the column \code{annotation} holds the pathway name.
#' @details The date the bundled snapshot was retrieved from KEGG is unknown.
#'   The script \file{data-raw/kegg.R} in the source repository rebuilds the
#'   object and stamps the retrieval date on the line below.
#'
#'   Retrieval date (kegg.db): unknown (bundled snapshot predates data-raw/).
#' @source KEGG database via KEGGREST (\code{keggList("pathway")}).
#' @usage data(kegg)
"kegg.db"

#' KEGG module data
#'
#' KEGG module identifiers and names, used by the module annotation builders.
#'
#' @format A data.frame with 473 rows and one column. The row names are KEGG
#'   module identifiers (for example \code{"M00001"}); the column
#'   \code{annotation} holds the module name.
#' @details The date the bundled snapshot was retrieved from KEGG is unknown.
#'   The script \file{data-raw/module.R} in the source repository rebuilds the
#'   object and stamps the retrieval date on the line below.
#'
#'   Retrieval date (module): unknown (bundled snapshot predates data-raw/).
#' @source KEGG database via KEGGREST (\code{keggList("module")}).
#' @usage data(module)
"module"

#' KEGG pathway hierarchy
#'
#' The KEGG pathway hierarchy: each pathway with its three classification
#' levels, used by \code{\link{richKEGG}}, \code{\link{richLevel}} and
#' \code{\link{richGSEA}} to add and aggregate pathway levels.
#'
#' @format A data.frame with 580 rows and four character columns: \code{ko}
#'   (KEGG pathway identifier, also the row names), \code{Level3} (pathway
#'   name), \code{Level2} (sub-category) and \code{Level1} (category).
#' @details The date the bundled snapshot was retrieved from KEGG is unknown.
#'   The script \file{data-raw/path.R} in the source repository rebuilds the
#'   object and stamps the retrieval date on the line below.
#'
#'   Retrieval date (path): unknown (bundled snapshot predates data-raw/).
#' @source KEGG pathway classification (BRITE hierarchy \code{br08901}). How the
#'   bundled snapshot was retrieved is unknown; \file{data-raw/path.R} rebuilds
#'   it from the KEGG REST API.
#' @usage data(path)
"path"

#' KEGG pathway list
#'
#' KEGG pathway identifiers and names. richR does not load this object. It has
#' 562 rows, five fewer than \code{kegg.db}, and its origin is not documented.
#'
#' @format A data.frame with 562 rows and one column. The row names are KEGG
#'   pathway identifiers (for example \code{"01100"}); the column
#'   \code{annotation} holds the pathway name.
#' @details The date the bundled snapshot was retrieved from KEGG is unknown.
#'   The script \file{data-raw/pathway.R} in the source repository rebuilds the
#'   object and stamps the retrieval date on the line below.
#'
#'   Retrieval date (pathway): unknown (bundled snapshot predates data-raw/).
#' @source KEGG database.
#' @usage data(pathway)
"pathway"
