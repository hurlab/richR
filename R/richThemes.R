#' Group reported terms into themes by the overlap of their evidence genes
#'
#' Groups the terms that one or several enrichment results report into
#' themes of related terms, by the kappa similarity of the genes that support
#' each term, and keeps everything needed to go back from a theme to its
#' terms, their rows, their genes and the input statistics. A theme is a
#' summary of the listed terms, not a statistical test.
#'
#' @section Interpretation:
#' \itemize{
#'   \item Theme grouping summarizes the listed terms and is \strong{not a
#'     statistical test}. No theme-level p-value exists: no element of the
#'     object carries a p-value, adjusted p-value, NES, enrichment score or
#'     any value computed from them, and the representative's p-value is not
#'     a p-value of the theme.
#'   \item Related or overlapping gene sets (for example the same pathway in
#'     two databases, or a parent and a child GO term) are \strong{not
#'     independent replications}: they share genes by construction, so the
#'     number of terms in a theme is not a measure of evidence.
#'   \item Evidence genes are the ORA overlap genes or the GSEA leading-edge
#'     genes of a term. They are the genes behind the term's result,
#'     \strong{not drivers} shown to cause anything.
#' }
#'
#' @section Grouping rule:
#' Every input follows the group rules of \code{\link{bindRich}}: a single
#' result is treated as \code{list(x)} and named \code{Group 1}; an unnamed
#' list is named \code{Group 1}, \code{Group 2}, ...; \code{NULL} entries are
#' failed groups; ORA and GSEA results cannot be mixed. Identifiers are
#' \code{Annot} (ORA) or \code{pathway} (GSEA); \code{Term} is display text
#' only.
#'
#' \describe{
#'   \item{Grouping set}{the identifiers whose \code{bindRich()} state is
#'     \code{"reported"} in at least one group. Terms that are only tested,
#'     not tested, \code{"untested_or_zero_overlap"} or \code{"unavailable"}
#'     are never grouped.}
#'   \item{Evidence genes}{of a term in a group, read from that group's own
#'     record: its \code{result(x)} row when the term is reported there, its
#'     tested-table row when it is tested but not reported there, nothing
#'     otherwise. ORA: the \code{GeneID} cell split on the object's
#'     \code{sep} (\code{","} when the slot is empty). GSEA: the
#'     \code{leadingEdge} cell (a list element, or text split on \code{sep}).
#'     Missing and empty strings are not genes; duplicates count once.}
#'   \item{Similarity genes}{of a grouped term: the union of its evidence
#'     genes over the groups that report it. They pool the supplied groups,
#'     so \strong{themes can change when a group is added or removed}.}
#'   \item{Background}{the union of the similarity genes of every grouped
#'     term; its size is \code{settings$background_size}.}
#'   \item{Similarity}{Cohen's kappa of two terms' similarity genes against
#'     the background (the package's kappa statistic, as used by
#'     \code{\link{richCluster}}). No p-value, adjusted p-value, NES or
#'     enrichment score is used anywhere in grouping.}
#'   \item{Clustering}{average-linkage hierarchical clustering of
#'     \code{1 - kappa}, cut at height \code{1 - cutoff}. The cut is
#'     inclusive: two clusters join when their average kappa is \strong{at
#'     least} \code{cutoff}, whereas \code{richCluster()} keeps pairs with
#'     kappa strictly greater than its cutoff, so the same number does not
#'     mean the same threshold in the two functions. Every grouped term
#'     belongs to exactly one theme; a term joined to no other is a singleton
#'     theme; no theme is removed by size or by any p-value.}
#'   \item{Degenerate kappa}{kappa is \code{0} when two terms share no gene,
#'     and also when the chance agreement is 1 — for example when every
#'     grouped term has the same genes and the background holds no other gene
#'     (kappa cannot be computed there). Such terms stay singletons; the
#'     small \code{background_size} shows why. A term without evidence genes
#'     is a singleton with no \code{genes} row.}
#'   \item{Representative}{the member with the highest mean kappa to the
#'     other members; means within \code{1e-12} of the highest are ties.
#'     Ties and singletons resolve to the first identifier in radix (C-locale)
#'     order.}
#'   \item{Label}{the representative's \code{Term} in the radix-first group
#'     name among the groups that report it with a non-missing, non-blank
#'     \code{Term} (\code{representative_term}); when there is none, the
#'     identifier. Labels are only display text: edit one with
#'     \code{th$themes$label[th$themes$theme_id == "theme_001"] <- "my label"};
#'     every link between tables uses \code{theme_id}.}
#' }
#'
#' Themes, theme ids, representatives, labels and the \code{similarity}
#' table do not depend on the order of the groups or of the rows within a
#' group. The \code{group} factor levels, the order of rows by group, the
#' order of the \code{members} columns after \code{state} and
#' \code{execution} and the levels of factor columns inherited from
#' \code{bindRich()}, and \code{settings$groups} follow the input order.
#'
#' @section Value:
#' A list of class \code{"richThemes"} with six elements:
#' \describe{
#'   \item{\code{themes}}{one row per theme: \code{theme_id} (\code{"theme_"}
#'     and the theme's position, zero-padded to at least three digits),
#'     \code{label}, \code{representative} (identifier),
#'     \code{representative_term} (may be \code{NA}), \code{n_terms},
#'     \code{n_groups_reported} (groups where at least one member is
#'     reported), \code{n_similarity_genes} (size of the union of the
#'     members' similarity genes, so genes of reported records only;
#'     \code{\link{themeGenes}} may list more, because it also shows
#'     tested-not-reported records) and \code{mean_within_kappa} (mean
#'     pairwise kappa among members, \code{NA} for a singleton). Themes are
#'     ordered by \code{n_terms} descending, then by representative.}
#'   \item{\code{members}}{one row per member and group: \code{theme_id},
#'     \code{is_representative}, then every \code{bindRich(x)} column in its
#'     order, with its values and classes unchanged. Rows are ordered by
#'     theme, identifier and group, and are exactly the \code{bindRich(x)}
#'     rows of the members.}
#'   \item{\code{genes}}{one row per member, group and evidence gene, for the
#'     groups where the member is reported or tested but not reported:
#'     \code{theme_id}, \code{group}, the identifier, \code{gene},
#'     \code{evidence_type} (\code{"ora_overlap"} or
#'     \code{"gsea_leading_edge"}) and \code{state}.}
#'   \item{\code{similarity}}{one row per unordered pair of members of one
#'     theme: \code{theme_id}, \code{identifier_a}, \code{identifier_b}
#'     (radix order), \code{kappa} and \code{n_shared} (similarity genes the
#'     two share).}
#'   \item{\code{ranking}}{GSEA: one row per group and ranked gene, the
#'     group's ranked input in stored order: \code{group}, \code{gene},
#'     \code{rank_score}. Zero rows for ORA.}
#'   \item{\code{settings}}{\code{method} (\code{"kappa"}), \code{distance}
#'     (\code{"1 - kappa"}), \code{linkage} (\code{"average"}),
#'     \code{cutoff}, \code{grouped_on}, \code{similarity_genes},
#'     \code{background} (the rules above, as text), \code{background_size},
#'     \code{representative}, \code{label}, \code{groups} (input order),
#'     \code{sources} (each group's \code{\link{enrichmentInfo}} record,
#'     \code{"failed"} for a \code{NULL} group) and \code{richR_version}.}
#' }
#' When nothing is reported in any group, the five tables have zero rows with
#' their columns and classes, and \code{settings} is filled. Every column
#' \code{richThemes()} creates is character, integer, double, logical or
#' factor; \code{members} keeps the classes of \code{bindRich()}.
#'
#' @param x a \code{richResult}, a \code{GSEAResult}, or a list of them (with
#'   \code{NULL} for failed groups), as accepted by \code{\link{bindRich}}.
#' @param cutoff the kappa threshold, one finite number in (0, 1]: clusters
#'   join when their average kappa is at least \code{cutoff}. Default 0.5.
#' @return a \code{richThemes} object, described in the Value section.
#' @seealso \code{\link{themeGenes}}, \code{\link{bindRich}},
#'   \code{\link{richCluster}}
#' @examples
#' \dontrun{
#' res <- batchEnrich(list(S1 = gene1, S2 = gene2), annot = hsako)
#' th <- richThemes(res$results)
#' th$themes
#' th$themes$label[1] <- "my label"
#' themeGenes(th, "theme_001")
#' }
#' @export
richThemes <- function(x, cutoff = 0.5) {
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) ||
      cutoff <= 0 || cutoff > 1) {
    stop("richThemes: cutoff must be one finite number in (0, 1].",
         call. = FALSE)
  }
  if (methods::is(x, "richResult") || methods::is(x, "GSEAResult")) {
    x <- list(x)
  }
  bound <- bindRich(x)
  if (is.null(names(x))) names(x) <- paste("Group", seq_along(x))
  groups <- names(x)
  gsea <- any(vapply(x, methods::is, logical(1), "GSEAResult"))
  id <- if (gsea) "pathway" else "Annot"
  added <- c("theme_id", "is_representative")
  clash <- intersect(names(bound), added)
  if (length(clash)) {
    stop("richThemes: result column(s) ", .bind_quote(clash),
         " clash with the columns richThemes() adds.", call. = FALSE)
  }

  ids <- as.character(bound[[id]])
  state <- bound$state
  ev <- .theme_evidence(bound, x, groups, gsea)

  # Grouping set, similarity genes and background (radix order throughout,
  # so nothing depends on the input order or the session locale).
  grouped <- sort(unique(ids[state %in% "reported"]), method = "radix")
  rep_rows <- which(state %in% "reported")
  sim_genes <- lapply(split(rep_rows, factor(ids[rep_rows], levels = grouped)),
                      function(r) {
                        g <- as.character(unlist(ev[r], use.names = FALSE))
                        sort(unique(g), method = "radix")
                      })
  names(sim_genes) <- grouped
  background <- sort(unique(as.character(unlist(sim_genes, use.names = FALSE))),
                     method = "radix")
  kap <- .theme_kappa(sim_genes, background)

  n <- length(grouped)
  cluster <- if (n >= 2L) {
    stats::cutree(stats::hclust(stats::as.dist(1 - kap), method = "average"),
                  h = 1 - cutoff)
  } else {
    rep(1L, n)
  }
  mem <- unname(split(grouped, cluster))
  rep_of <- vapply(mem, .theme_representative, character(1), kap = kap)
  ord <- order(-lengths(mem), rep_of, method = "radix")
  mem <- mem[ord]
  rep_of <- rep_of[ord]
  n_themes <- length(mem)
  width <- max(3L, nchar(as.character(n_themes)))
  theme_ids <- sprintf(paste0("theme_%0", width, "d"), seq_len(n_themes))

  rep_term <- vapply(rep_of, .theme_rep_term, character(1), ids = ids,
                     state = state, bound = bound)
  rep_term <- unname(rep_term)
  label <- ifelse(is.na(rep_term), rep_of, rep_term)
  themes <- .theme_frame(list(
    theme_id = theme_ids,
    label = as.character(label),
    representative = as.character(rep_of),
    representative_term = as.character(rep_term),
    n_terms = as.integer(lengths(mem)),
    n_groups_reported = vapply(mem, function(m) {
      length(unique(as.character(bound$group[ids %in% m & state %in% "reported"])))
    }, integer(1)),
    n_similarity_genes = vapply(mem, function(m) {
      length(unique(unlist(sim_genes[m], use.names = FALSE)))
    }, integer(1)),
    mean_within_kappa = vapply(mem, function(m) {
      if (length(m) < 2L) return(NA_real_)
      k <- kap[m, m, drop = FALSE]
      mean(k[upper.tri(k)])
    }, numeric(1))
  ))

  # members: the bindRich() rows of every member in every group.
  member_ids <- unlist(mem, use.names = FALSE)
  theme_index <- rep(seq_len(n_themes), lengths(mem))
  theme_of <- function(v) theme_index[match(v, member_ids)]
  sel <- which(ids %in% grouped)
  sel <- sel[order(theme_of(ids[sel]), match(ids[sel], grouped),
                   as.integer(bound$group[sel]))]
  members <- .theme_frame(c(
    list(theme_id = theme_ids[theme_of(ids[sel])],
         is_representative = ids[sel] %in% rep_of),
    as.list(bound[sel, , drop = FALSE])
  ))

  # genes: each member's evidence genes in each group where it is reported
  # or tested but not reported, in the members' row order.
  gsel <- sel[state[sel] %in% c("reported", "tested_not_reported")]
  len <- lengths(ev[gsel])
  genes <- .theme_frame(stats::setNames(list(
    theme_ids[rep(theme_of(ids[gsel]), len)],
    rep(bound$group[gsel], len),
    rep(ids[gsel], len),
    as.character(unlist(ev[gsel], use.names = FALSE)),
    rep(if (gsea) "gsea_leading_edge" else "ora_overlap", sum(len)),
    rep(state[gsel], len)
  ), c("theme_id", "group", id, "gene", "evidence_type", "state")))

  similarity <- .theme_similarity(mem, theme_ids, kap, sim_genes)
  ranking <- .theme_ranking(x, groups, gsea)

  settings <- list(
    method = "kappa",
    distance = "1 - kappa",
    linkage = "average",
    cutoff = cutoff,
    grouped_on = paste("identifiers whose bindRich() state is 'reported'",
                       "in at least one group"),
    similarity_genes = paste(
      "a term's evidence genes unioned over the groups that report it;",
      if (gsea) "evidence genes are the leadingEdge genes (GSEA)" else
        "evidence genes are the GeneID overlap genes (ORA)",
      "of the group's own record"),
    background = "the union of the similarity genes of every grouped term",
    background_size = length(background),
    representative = paste(
      "the member with the highest mean kappa to the other members; means",
      "within 1e-12 of the highest are ties; ties and singletons resolve to",
      "the first identifier in radix order"),
    label = paste(
      "the representative's Term in the radix-first group that reports it",
      "with a non-missing, non-blank Term; otherwise the identifier"),
    groups = groups,
    sources = lapply(x, function(o) if (is.null(o)) "failed" else enrichmentInfo(o)),
    richR_version = as.character(utils::packageVersion("richR"))
  )

  structure(list(themes = themes, members = members, genes = genes,
                 similarity = similarity, ranking = ranking,
                 settings = settings),
            class = "richThemes")
}

#' A data.frame from a named list of equal-length columns, with compact row
#' names and list columns kept as they are
#' @keywords internal
#' @noRd
.theme_frame <- function(cols) {
  n <- if (length(cols)) length(cols[[1L]]) else 0L
  structure(cols, names = names(cols), row.names = .set_row_names(n),
            class = "data.frame")
}

#' The separator a result object records, or "," when it records none
#' @keywords internal
#' @noRd
.theme_sep <- function(obj) {
  s <- obj@sep
  if (length(s) == 1L && !is.na(s) && nzchar(s)) s else ","
}

#' Evidence genes of every bindRich() row that is reported or tested but not
#' reported, read from the row's own record (bindRich() takes each value from
#' the group's result or tested table); character(0) for every other row
#' @keywords internal
#' @noRd
.theme_evidence <- function(bound, x, groups, gsea) {
  out <- rep(list(character(0)), nrow(bound))
  col <- if (gsea) "leadingEdge" else "GeneID"
  if (!col %in% names(bound)) return(out)
  v <- bound[[col]]
  seps <- vapply(x, function(o) if (is.null(o)) "," else .theme_sep(o), character(1))
  gi <- as.integer(bound$group)
  for (r in which(bound$state %in% c("reported", "tested_not_reported"))) {
    if (is.list(v)) {
      g <- as.character(v[[r]])
    } else {
      cell <- as.character(v[[r]])
      g <- if (is.na(cell)) character(0) else
        strsplit(cell, seps[[gi[[r]]]], fixed = TRUE)[[1L]]
    }
    g <- g[!is.na(g) & nzchar(g)]
    out[[r]] <- unique(g)
  }
  out
}

#' Cohen's kappa between every pair of gene sets against a background,
#' computed from the 2 x 2 table exactly as .kappa() does, on gene vectors
#' @keywords internal
#' @noRd
.theme_kappa <- function(sets, background) {
  n <- length(sets)
  N <- length(background)
  inc <- matrix(0, N, n)
  for (i in seq_len(n)) inc[match(sets[[i]], background), i] <- 1
  n11 <- crossprod(inc)
  n10 <- matrix(colSums(inc), n, n) - n11
  n01 <- t(n10)
  n00 <- N - n11 - n10 - n01
  oab <- (n11 + n00) / N
  aab <- ((n11 + n10) * (n11 + n01) + (n01 + n00) * (n10 + n00)) / (N * N)
  kap <- (oab - aab) / (1 - aab)
  kap[n11 == 0 | aab == 1] <- 0
  diag(kap) <- 1
  dimnames(kap) <- list(names(sets), names(sets))
  kap
}

#' The member with the highest mean kappa to the others (ties within 1e-12
#' resolve to the first identifier; members arrive in radix order)
#' @keywords internal
#' @noRd
.theme_representative <- function(m, kap) {
  if (length(m) == 1L) return(m)
  k <- kap[m, m, drop = FALSE]
  diag(k) <- NA
  mk <- rowMeans(k, na.rm = TRUE)
  m[mk >= max(mk) - 1e-12][[1L]]
}

#' The representative's Term in the radix-first group that reports it with a
#' non-missing, non-blank Term, or NA
#' @keywords internal
#' @noRd
.theme_rep_term <- function(r, ids, state, bound) {
  rows <- which(ids == r & state %in% "reported")
  term <- as.character(bound$Term[rows])
  grp <- as.character(bound$group[rows])
  ok <- !is.na(term) & nzchar(trimws(term))
  if (!any(ok)) return(NA_character_)
  term[ok][order(grp[ok], method = "radix")][[1L]]
}

#' One row per unordered pair of members of one theme
#' @keywords internal
#' @noRd
.theme_similarity <- function(mem, theme_ids, kap, sim_genes) {
  parts <- lapply(seq_along(mem), function(t) {
    m <- mem[[t]]
    if (length(m) < 2L) return(NULL)
    p <- which(upper.tri(diag(length(m))), arr.ind = TRUE)
    p <- p[order(p[, 1L], p[, 2L]), , drop = FALSE]
    a <- m[p[, 1L]]
    b <- m[p[, 2L]]
    list(theme_id = rep(theme_ids[[t]], length(a)), identifier_a = a,
         identifier_b = b, kappa = kap[cbind(a, b)],
         n_shared = vapply(seq_along(a), function(j) {
           length(intersect(sim_genes[[a[[j]]]], sim_genes[[b[[j]]]]))
         }, integer(1)))
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  .theme_frame(list(
    theme_id = as.character(unlist(lapply(parts, `[[`, "theme_id"))),
    identifier_a = as.character(unlist(lapply(parts, `[[`, "identifier_a"))),
    identifier_b = as.character(unlist(lapply(parts, `[[`, "identifier_b"))),
    kappa = as.numeric(unlist(lapply(parts, `[[`, "kappa"))),
    n_shared = as.integer(unlist(lapply(parts, `[[`, "n_shared")))
  ))
}

#' Every GSEA group's ranked input in stored order; zero rows for ORA
#' @keywords internal
#' @noRd
.theme_ranking <- function(x, groups, gsea) {
  keep <- if (gsea) which(!vapply(x, is.null, logical(1))) else integer(0)
  inputs <- lapply(x[keep], function(o) o@input)
  gene <- lapply(inputs, function(v) {
    if (is.null(names(v))) rep(NA_character_, length(v)) else names(v)
  })
  .theme_frame(list(
    group = factor(rep(groups[keep], lengths(inputs)), levels = groups),
    gene = as.character(unlist(gene, use.names = FALSE)),
    rank_score = as.numeric(unlist(lapply(inputs, unname), use.names = FALSE))
  ))
}
