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
#' @param cor_method \code{"spearman"} (default, robust to cross-species scale
#'   differences) or \code{"pearson"}.
#' @param min_correspondence Minimum correlation for a match to be reported.
#'   Default 0.3; see \code{plot_correspondence_elbow}.
#' @param assay,slot Expression to average: the log-normalized \code{data} slot
#'   of the \code{RNA} assay by default.
#' @param verbose Print progress. Default TRUE.
#'
#' @return A \code{cluster_correspondence}: data.frame of cluster_a, cluster_b,
#'   score, reciprocal (logical). Reciprocal best matches passing
#'   \code{min_correspondence} are marked \code{reciprocal = TRUE}. The full
#'   correlation matrix is attached as \code{attr(x, "cor_matrix")}.
#'
#' @export
match_clusters <- function(clusters_a, clusters_b, homology_graph,
                           species_a = NULL, species_b = NULL,
                           cor_method = c("spearman", "pearson"),
                           min_correspondence = 0.3,
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
  pairs$reciprocal <- mapply(function(a, b) identical(best_a_for_b[[b]], a),
                             pairs$cluster_a, pairs$cluster_b) &
                      pairs$score >= min_correspondence

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
  cl <- Seurat::Idents(obj)
  # mean expression per cluster
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
