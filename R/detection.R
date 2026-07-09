#' Cross-species homolog correlations over the metacell grid
#'
#' For each focal gene in species A, correlates its expression across matched
#' metacell pairs against its ortholog and against each paralog of that ortholog
#' in species B. Correlations run in real gene space over the fine metacell grid
#' (\code{\link{match_metacells}}), giving many paired observations per gene.
#'
#' @param metacells_a,metacells_b \code{metacells} objects for species A and B.
#' @param metacell_correspondence A \code{metacell_correspondence}
#'   (\code{\link{match_metacells}}) — the grid of matched pairs.
#' @param homology_graph A \code{homology_graph} with orthologs AND paralogs.
#' @param species_a,species_b Species identifiers matching the graph.
#' @param focal_genes Optional character vector of species-A genes to restrict
#'   to. NULL = all genes with an ortholog and at least one paralog.
#' @param cor_method \code{"spearman"} (default) or \code{"pearson"}.
#' @param min_pairs Minimum matched pairs required to compute a correlation.
#' @param verbose Print progress. Default TRUE.
#'
#' @export
compute_homolog_correlations <- function(metacells_a, metacells_b,
                                         metacell_correspondence,
                                         homology_graph,
                                         species_a = NULL, species_b = NULL,
                                         focal_genes = NULL,
                                         cor_method = c("spearman", "pearson"),
                                         min_pairs = 10,
                                         verbose = TRUE) {
  cor_method <- match.arg(cor_method)
  mmc <- metacell_correspondence
  if (nrow(mmc) < min_pairs)
    stop("Too few matched metacell pairs (", nrow(mmc), ") for correlation.",
         call. = FALSE)

  Da <- metacells_a$data   # genes x metacells (A)
  Db <- metacells_b$data   # genes x metacells (B)

  # expression vectors aligned to the grid: for each matched pair, the A-gene
  # value from metacell_a and the B-gene value from metacell_b.
  mca <- mmc$metacell_a; mcb <- mmc$metacell_b

  hg <- homology_graph
  # orient so gene_a = species A
  flip_needed <- !is.null(species_a) && !is.null(species_b) &&
    any(hg$relationship=="ortholog" & hg$species_a==species_b)
  orient <- function(e) if (flip_needed) .flip_edges(e) else e

  orth <- orient(hg[hg$relationship == "ortholog" &
                    !is.na(hg$ortholog_type) & hg$ortholog_type=="one2one", ])
  # paralogs are within-species; we need species-B paralogs (partners of the
  # ortholog on the B side). Keep B-side paralog edges.
  para_b <- hg[hg$relationship == "paralog" & hg$species_a == species_b, ]

  # focal set: A genes that have a one2one ortholog in B, present in both grids
  cand <- orth[orth$gene_a %in% rownames(Da) & orth$gene_b %in% rownames(Db), ]
  if (!is.null(focal_genes)) cand <- cand[cand$gene_a %in% focal_genes, ]
  if (nrow(cand) == 0) stop("No focal genes with usable orthologs.", call.=FALSE)

  cor_vec <- function(ga, gb) {
    if (!(ga %in% rownames(Da)) || !(gb %in% rownames(Db))) return(c(NA, 0))
    xa <- Da[ga, mca]; xb <- Db[gb, mcb]
    ok <- is.finite(xa) & is.finite(xb)
    if (sum(ok) < min_pairs) return(c(NA, sum(ok)))
    if (stats::sd(xa[ok])==0 || stats::sd(xb[ok])==0) return(c(NA, sum(ok)))
    c(stats::cor(xa[ok], xb[ok], method = cor_method), sum(ok))
  }

  rows <- list()
  for (i in seq_len(nrow(cand))) {
    ga <- cand$gene_a[i]; gb_ortho <- cand$gene_b[i]

    # ortholog correlation: A gene vs its B ortholog
    ro <- cor_vec(ga, gb_ortho)
    rows[[length(rows)+1L]] <- data.frame(
      gene_a=ga, gene_b=gb_ortho, relationship="ortholog",
      r=ro[1], n_pairs=ro[2], stringsAsFactors=FALSE)

    # paralog correlations: A gene vs each B-side paralog of the ortholog
    b_paralogs <- unique(c(
      para_b$gene_b[para_b$gene_a == gb_ortho],
      para_b$gene_a[para_b$gene_b == gb_ortho]))
    b_paralogs <- setdiff(b_paralogs, gb_ortho)
    for (gp in b_paralogs) {
      rp <- cor_vec(ga, gp)
      rows[[length(rows)+1L]] <- data.frame(
        gene_a=ga, gene_b=gp, relationship="paralog",
        r=rp[1], n_pairs=rp[2], stringsAsFactors=FALSE)
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (isTRUE(verbose)) message(nrow(out), " correlations for ",
                               length(unique(out$gene_a)), " focal genes over ",
                               nrow(mmc), " metacell pairs.")
    attr(out, "params") <- list(min_pairs = min_pairs, cor_method = cor_method)  
    class(out) <- c("homolog_correlations", "data.frame")
  out
}
#' Detect paralog substitution candidates
#'
#' For each focal gene, computes delta_r = max(paralog r) - ortholog r and
#' flags candidates where a paralog tracks the cross-species pattern better than
#' the annotated ortholog. The default threshold (0.3) matches Tarashansky et
#' al. 2021; plot_substitution_elbow exposes the distribution behind it.
#'
#' @param homolog_correlations A \code{homolog_correlations}
#'   (\code{\link{compute_homolog_correlations}}).
#' @param delta_threshold Minimum delta_r to flag. Default 0.3.
#' @param require_ortholog If TRUE (default), only score genes that have an
#'   ortholog correlation (drops paralog-only or NA-ortholog genes).
#'
#' @return A \code{substitutions} data.frame: gene, ortholog, n_orthologs,
#'   best_paralog, r_ortholog, r_paralog, delta_r, flagged. Sorted by
#'   descending delta_r. \code{n_orthologs} is the number of ortholog edges the
#'   focal gene has in the homology graph; \code{ortholog} and
#'   \code{r_ortholog} report the first of them, so \code{n_orthologs > 1}
#'   marks a one-to-many or many-to-many focal gene whose delta_r depends on
#'   which ortholog was taken. Ambiguity is reported rather than resolved.
#'   Parameters are recorded in \code{attr(x, "params")}. 
#'   Focal genes that cannot be scored are not silently discarded: they are
#'   recorded in \code{attr(x, "dropped")} with a reason. \code{"ortholog_invariant"}
#'   means the ortholog was observed across the full grid but shows no variance —
#'   typically unexpressed in species B. Because \code{delta_r} subtracts a
#'   measured ortholog correlation, a substitution in which the ortholog has been
#'   wholly lost or silenced is undefined rather than maximal, and appears here
#'   rather than among the flagged genes.
#'
#' @export
detect_substitutions <- function(homolog_correlations,
                                 delta_threshold = 0.3,
                                 require_ortholog = TRUE) {
hc <- homolog_correlations
  hcp <- attr(hc, "params")
  min_pairs <- if (!is.null(hcp) && !is.null(hcp$min_pairs)) hcp$min_pairs else NA_integer_
  genes <- unique(hc$gene_a)
  rows <- list(); dropped <- list()

  .drop <- function(g, reason, orth, para) {
    dropped[[length(dropped) + 1L]] <<- data.frame(
      gene         = g,
      reason       = reason,
      ortholog     = if (nrow(orth)) orth$gene_b[1]  else NA_character_,
      n_pairs      = if (nrow(orth)) orth$n_pairs[1] else NA_integer_,
      best_paralog = if (nrow(para)) para$gene_b[which.max(para$r)] else NA_character_,
      r_paralog    = if (nrow(para)) max(para$r)     else NA_real_,
      stringsAsFactors = FALSE)
  }

  for (g in genes) {
    sub <- hc[hc$gene_a == g, ]
    orth <- sub[sub$relationship == "ortholog", ]
    para <- sub[sub$relationship == "paralog" & is.finite(sub$r), ]
    r_orth <- if (nrow(orth) && is.finite(orth$r[1])) orth$r[1] else NA

    if (require_ortholog && is.na(r_orth)) {
      reason <- if (nrow(orth) == 0) {
        "no_ortholog_edge"
      } else if (!is.na(min_pairs) && orth$n_pairs[1] < min_pairs) {
        "ortholog_too_sparse"
      } else {
        # full pair coverage but no correlation: one side is invariant across
        # the grid. Typically the ortholog is unexpressed in species B — the
        # limiting case of substitution, which delta_r cannot score.
        "ortholog_invariant"
      }
      .drop(g, reason, orth, para); next
    }
    if (nrow(para) == 0) { .drop(g, "no_paralog_correlation", orth, para); next }
    best <- para[which.max(para$r), ]
    delta <- best$r - (if (is.na(r_orth)) 0 else r_orth)
    rows[[length(rows)+1L]] <- data.frame(
      gene = g,
      ortholog = if (nrow(orth)) orth$gene_b[1] else NA,
      n_orthologs = nrow(orth),
      best_paralog = best$gene_b,
      r_ortholog = r_orth,
      r_paralog = best$r,
      delta_r = delta,
      flagged = delta >= delta_threshold,
      stringsAsFactors = FALSE)
  }
out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) {
    out <- data.frame(
      gene = character(0), ortholog = character(0),
      n_orthologs = integer(0),
      best_paralog = character(0), r_ortholog = numeric(0),
      r_paralog = numeric(0), delta_r = numeric(0),
      flagged = logical(0), stringsAsFactors = FALSE)
  } else {
    rownames(out) <- NULL
    out <- out[order(-out$delta_r), ]
  }
attr(out, "dropped") <- if (length(dropped))
    do.call(rbind, dropped) else
    data.frame(gene = character(0), reason = character(0), ortholog = character(0),
               n_pairs = integer(0), best_paralog = character(0),
               r_paralog = numeric(0), stringsAsFactors = FALSE) 
attr(out, "params") <- list(delta_threshold = delta_threshold,
                              require_ortholog = require_ortholog)
  class(out) <- c("substitutions", "data.frame")
  out
}

#' @export
summary.substitutions <- function(object, ...) {
  cat("substitutions:", nrow(object), "focal genes scored\n")
  cat("  flagged (delta_r >= threshold):", sum(object$flagged), "\n")
  cat("  delta_r range:", sprintf("%.2f", min(object$delta_r)), "-",
      sprintf("%.2f", max(object$delta_r)), "\n")
d <- attr(object, "dropped")
  if (!is.null(d) && nrow(d)) {
    cat("\n  ", nrow(d), " focal gene(s) not scored:\n", sep = "")
    for (r in unique(d$reason))
      cat("    ", r, ": ", paste(d$gene[d$reason == r], collapse = ", "), "\n", sep = "")
  }  
  invisible(object)
}
