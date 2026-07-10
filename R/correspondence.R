#' Cross-species cluster correspondence from expression
#'
#' Determines which cluster in one species corresponds to which in another,
#' from expression alone, without cell-type labels. Cluster mean-expression
#' profiles are computed over the one-to-one ortholog set, gene axes aligned
#' across species through the homology graph, clusters correlated, and
#' reciprocal best matches retained.
#'
#' @param clusters_a,clusters_b Clustered Seurat objects (single objects, or
#'   single-element results of \code{\link{cluster_species}}) for species A and B.
#' @param homology_graph A \code{homology_graph} (\code{\link{build_homology_graph}}).
#' @param species_a,species_b Species identifiers matching the graph's
#'   \code{species_a}/\code{species_b} columns. If NULL, inferred from the graph
#'   when it contains exactly one cross-species pair.
#' @param cor_method \code{"pearson"} (default) or \code{"spearman"}. Spearman
#'   was formerly the default, on the grounds that absolute expression scales
#'   differ across species and platforms. \code{center_genes} and
#'   \code{scale_genes} now remove that difference directly, so Pearson on
#'   standardized profiles is both the closer analog of SAMap's preprocessing
#'   and the cheaper computation. \code{\link{compute_homolog_correlations}}
#'   keeps Spearman by default: it correlates one gene across metacell pairs,
#'   an axis standardization does not touch, where a single high-expressing
#'   metacell can dominate a Pearson fit.
#' @param min_correspondence Minimum correlation for a match to be reported.
#'   Default 0.3; see \code{plot_correspondence_elbow}. Note that this gates the
#'   \code{reciprocal} flag only; every cluster pair is scored and returned.
#' @param center_genes Subtract each gene's across-cluster mean within each
#'   species before correlating (default TRUE). Uncentered profiles are dominated
#'   by the shared abundance baseline and yield near-uniform, uninformative
#'   scores; set FALSE only to reproduce that behaviour.
#' @param scale_genes Divide each gene by its across-cluster standard deviation
#'   within each species (default TRUE), after centering. Centering and scaling
#'   together reproduce the standardization SAMap applies before manifold
#'   projection. Genes with zero variance in either species are dropped.
#' @param assay,slot Expression to average: the log-normalized \code{data} slot
#'   of the \code{RNA} assay by default.
#' @param verbose Print progress. Default TRUE.
#'
#' @return A \code{cluster_correspondence}: data.frame of cluster_a, cluster_b,
#'   score, and three logicals. \code{mutual_best} is TRUE when each cluster is
#'   the other's highest-scoring partner, a structural property of the
#'   correlation matrix. \code{passes} is TRUE when the score reaches
#'   \code{min_correspondence}. \code{reciprocal}, which downstream stages
#'   consume, is their conjunction. One row per A-cluster: every A-cluster's best
#'   B-cluster is reported, so absence from the table means the cluster does not
#'   exist, not that it scored poorly. The full correlation matrix is attached as
#'   \code{attr(x, "cor_matrix")}.
#' @export
match_clusters <- function(clusters_a, clusters_b, homology_graph,
                           species_a = NULL, species_b = NULL,
                           cor_method = c("pearson", "spearman"),
                           min_correspondence = 0.3,
                           center_genes = TRUE,
                           scale_genes = TRUE,
                           assay = "RNA", slot = "data",
                           verbose = TRUE) {
  cor_method <- match.arg(cor_method)

  obj_a <- .unwrap_clustered(clusters_a)
  obj_b <- .unwrap_clustered(clusters_b)

  # ---- pick the one-to-one ortholog bridge ----------------------------------
  hg <- homology_graph
  o2o <- hg[hg$relationship == "ortholog" &
            !is.na(hg$ortholog_type) & hg$ortholog_type == "one2one", ,
            drop = FALSE]
  if (nrow(o2o) == 0) stop("No one-to-one orthologs in the homology graph.",
                           call. = FALSE)

  # orient graph so gene_a is species A of THIS call. The graph is directed
  # (focal species_a); if the caller's species_a matches graph species_b, flip.
  if (!is.null(species_a) && !is.null(species_b)) {
    if (all(o2o$species_a == species_b) && all(o2o$species_b == species_a)) {
      o2o <- .flip_edges(o2o)
    }
  }

  # keep orthologs whose genes are present in BOTH objects
  ga <- rownames(obj_a); gb <- rownames(obj_b)
  keep <- o2o$gene_a %in% ga & o2o$gene_b %in% gb
  o2o <- o2o[keep, , drop = FALSE]
  # de-duplicate: one row per gene_a and per gene_b (defensive; o2o should be 1:1)
  o2o <- o2o[!duplicated(o2o$gene_a) & !duplicated(o2o$gene_b), , drop = FALSE]
if (nrow(o2o) == 0) {
    stop("No one-to-one orthologs bridge the two objects. Check that the ",
         "homology graph's gene IDs match the objects' rownames and that ",
         "species_a/species_b are correct.", call. = FALSE)
  }
if (isTRUE(verbose)) message("Bridging on ", nrow(o2o),
                               " one-to-one orthologs present in both objects.")
  if (nrow(o2o) < 100) warning("Fewer than 100 shared orthologs; ",
                               "correspondence may be unstable.", call. = FALSE)

  # ---- cluster mean-expression profiles over the bridge ---------------------
  prof_a <- .cluster_profiles(obj_a, o2o$gene_a, assay, slot)  # genes x clusters
  prof_b <- .cluster_profiles(obj_b, o2o$gene_b, assay, slot)
  if (isTRUE(center_genes)) {
    # Subtract each gene's mean across that species' clusters. Without this,
    # correlations are dominated by the shared abundance baseline (housekeeping
    # genes high everywhere), every cluster pair scores 0.5-0.8, and best-match
    # assignment is decided by third-decimal differences. Centering leaves each
    # cluster's deviation from its own species' average, which is the quantity
    # that is comparable across species.
    prof_a <- prof_a - rowMeans(prof_a)
    prof_b <- prof_b - rowMeans(prof_b)
  }
  if (isTRUE(scale_genes)) {
    # Divide each gene by its across-cluster SD, so that a gene's vote is its
    # pattern across clusters rather than its dynamic range. Together with
    # centering this is the zero-mean/unit-variance standardization SAMap
    # applies before projecting into the joint space. Genes invariant across
    # clusters have SD 0 and carry no signal; drop rather than divide by zero.
    sd_a <- apply(prof_a, 1, stats::sd)
    sd_b <- apply(prof_b, 1, stats::sd)
    keep <- sd_a > 0 & sd_b > 0
    if (sum(keep) < 50)
      stop("Fewer than 50 bridge genes vary across clusters in both species.",
           call. = FALSE)
    if (isTRUE(verbose) && any(!keep))
      message("Dropping ", sum(!keep), " bridge genes invariant across clusters.")
    prof_a <- prof_a[keep, , drop = FALSE] / sd_a[keep]
    prof_b <- prof_b[keep, , drop = FALSE] / sd_b[keep]
  }
  # rows are aligned: row i of prof_a is gene_a[i], row i of prof_b is its
  # ortholog gene_b[i]. So the matrices share a (species-bridged) gene axis.
  stopifnot(nrow(prof_a) == nrow(prof_b))

  # ---- correlate every A-cluster against every B-cluster --------------------
  cormat <- stats::cor(prof_a, prof_b, method = cor_method)
  rownames(cormat) <- colnames(prof_a)
  colnames(cormat) <- colnames(prof_b)

  # ---- reciprocal best matches ----------------------------------------------
  best_b_for_a <- colnames(cormat)[max.col(cormat, ties.method = "first")]
  best_a_for_b <- rownames(cormat)[max.col(t(cormat), ties.method = "first")]
  names(best_b_for_a) <- rownames(cormat)
  names(best_a_for_b) <- colnames(cormat)

  pairs <- data.frame(
    cluster_a = rownames(cormat),
    cluster_b = best_b_for_a,
    score     = cormat[cbind(seq_len(nrow(cormat)),
                             match(best_b_for_a, colnames(cormat)))],
    stringsAsFactors = FALSE
  )
  # Two independent properties, kept separate. mutual_best is structural: A's
  # best B also names A as its best. passes is a magnitude judgement against
  # min_correspondence. Bundling them makes it impossible to ask whether the
  # matching found the right pair separately from whether the score was large,
  # which is the question when the correlation scale changes (e.g. under
  # standardization). `reciprocal` retains the old meaning: both must hold.
  pairs$mutual_best <- mapply(function(a, b) identical(best_a_for_b[[b]], a),
                              pairs$cluster_a, pairs$cluster_b)
  pairs$passes     <- pairs$score >= min_correspondence
  pairs$reciprocal <- pairs$mutual_best & pairs$passes

  pairs <- pairs[order(-pairs$score), ]
  rownames(pairs) <- NULL
  attr(pairs, "cor_matrix") <- cormat
  attr(pairs, "n_bridge_genes") <- nrow(o2o)
  class(pairs) <- c("cluster_correspondence", "data.frame")
  pairs
}

# ---- helpers ---------------------------------------------------------------
.unwrap_clustered <- function(x) {
  if (inherits(x, "Seurat")) return(x)
  if (inherits(x, "species_clusters") || is.list(x)) {
    if (length(x) != 1) stop("Pass a single clustered object per species.",
                             call. = FALSE)
    return(x[[1]])
  }
  stop("Expected a Seurat object or single-element cluster result.", call. = FALSE)
}

.flip_edges <- function(e) {
  e[c("gene_a","gene_b")]       <- e[c("gene_b","gene_a")]
  e[c("species_a","species_b")] <- e[c("species_b","species_a")]
  e
}

.cluster_profiles <- function(obj, genes, assay, slot) {
  m <- Seurat::GetAssayData(obj, assay = assay, layer = slot)[genes, , drop = FALSE]
  cl <- droplevels(Seurat::Idents(obj))  # drop empty levels: rowMeans over
  # zero cells yields NaN, which poisons an entire row/column of cormat
  clusters <- levels(cl)
  prof <- vapply(clusters, function(k) {
    cells <- which(cl == k)
    Matrix::rowMeans(m[, cells, drop = FALSE])
  }, numeric(nrow(m)))
  colnames(prof) <- clusters
  prof
}

#' @export
summary.cluster_correspondence <- function(object, ...) {
  n_recip <- sum(object$reciprocal)
  cat("cluster_correspondence:\n")
  cat("  ortholog bridge:", .bridge_n(object), "one2one genes\n")
  cat("  A-clusters:", nrow(object), "| reciprocal matches:", n_recip, "\n")
  if (!is.null(object$mutual_best))
    cat("  mutually best:", sum(object$mutual_best),
        "| of those, passing threshold:", sum(object$mutual_best & object$passes), "\n")
  cat("  score range:", sprintf("%.2f", min(object$score)), "-",
      sprintf("%.2f", max(object$score)), "\n")
  invisible(object)
}

# Bridge size, tolerant of objects serialized before this attribute was renamed
# from "n_orthologs" to "n_bridge_genes". `n_orthologs` is now a per-gene column
# on `substitutions` and means something different.
.bridge_n <- function(x) {
  n <- attr(x, "n_bridge_genes")
  if (is.null(n)) attr(x, "n_orthologs") else n
}
