#' Build metacells per species by coarse-graining within clusters
#'
#' Coarse-grains each species into metacells with SuperCell (Bilous et al.
#' 2022), preserving within-cluster manifold structure rather than collapsing
#' each cluster to a single pseudobulk profile. Graining is performed
#' \emph{within each cluster} so every metacell is pure to one cluster by
#' construction — the substrate \code{\link{match_metacells}} requires.
#'
#' @param clusters A clustered Seurat object, or a single-element
#'   \code{\link{cluster_species}} result.
#' @param gamma Graining level: approximately \code{gamma} cells per metacell.
#'   Default 20 (SuperCell convention). Select deliberately with
#'   \code{plot_graining_curve}, not by substitution count.
#' @param assay,slot Expression to grain: log-normalized \code{data} of the
#'   \code{RNA} assay by default.
#' @param n_var_genes Variable genes used for the kNN graph SuperCell builds.
#' @param k_knn Nearest neighbors for SuperCell's graph. Default 5.
#' @param min_cells Clusters smaller than this are kept whole (one metacell)
#'   rather than grained, avoiding degenerate tiny-cluster graining.
#' @param cluster_col Metadata column holding cluster identity. Defaults to the
#'   active idents.
#' @param verbose Print progress. Default TRUE.
#'
#' @return A \code{metacells} object: a list with
#'   \item{counts}{genes x metacells aggregated (summed) raw counts matrix}
#'   \item{data}{genes x metacells averaged log-normalized expression}
#'   \item{meta}{data.frame per metacell: metacell_id, cluster, n_cells}
#'   \item{membership}{named vector mapping each original cell to its metacell}
#'   Parameters recorded in \code{attr(x, "params")}.
#'
#' @export
build_metacells <- function(clusters,
                            gamma = 20,
                            assay = "RNA", slot = "data",
                            n_var_genes = 2000, k_knn = 5,
                            min_cells = 10,
                            cluster_col = NULL,
                            verbose = TRUE) {
  if (!requireNamespace("SuperCell", quietly = TRUE)) {
    stop("build_metacells requires the 'SuperCell' package.", call. = FALSE)
  }
  obj <- .unwrap_clustered(clusters)  # reuse helper from correspondence.R

  cl <- if (is.null(cluster_col)) Seurat::Idents(obj) else
    factor(obj[[cluster_col, drop = TRUE]])
  clusters_lvl <- levels(cl)

  data_mat   <- Seurat::GetAssayData(obj, assay = assay, layer = slot)
  counts_mat <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")

  all_membership <- character(0)
  mc_counts_list <- list()
  mc_data_list   <- list()
  mc_meta_list   <- list()
  mc_counter <- 0L

  for (k in clusters_lvl) {
    cells_k <- names(cl)[cl == k]
    n_k <- length(cells_k)
    if (n_k == 0) next

    if (n_k < min_cells) {
      # too small to grain: keep as a single metacell
      mc_counter <- mc_counter + 1L
      mid <- paste0("mc", mc_counter)
      all_membership[cells_k] <- mid
      mc_counts_list[[mid]] <- Matrix::rowSums(counts_mat[, cells_k, drop = FALSE])
      mc_data_list[[mid]]   <- Matrix::rowMeans(data_mat[, cells_k, drop = FALSE])
      mc_meta_list[[mid]]   <- data.frame(metacell_id = mid, cluster = k,
                                          n_cells = n_k, stringsAsFactors = FALSE)
      next
    }

    # cap gamma so we don't request more coarsening than cells allow
    gamma_k <- min(gamma, floor(n_k / 2))
    if (isTRUE(verbose)) message("  cluster ", k, ": ", n_k,
                                 " cells -> gamma ", gamma_k)

    ge_k <- data_mat[, cells_k, drop = FALSE]
    sc <- SuperCell::SCimplify(
      ge_k, gamma = gamma_k, k.knn = k_knn, n.var.genes = n_var_genes,
      genes.use = NULL
    )
    memb <- sc$membership  # integer per cell, 1..n_metacells_in_cluster

    # aggregate within this cluster
    for (m in sort(unique(memb))) {
      mc_counter <- mc_counter + 1L
      mid <- paste0("mc", mc_counter)
      cells_m <- cells_k[memb == m]
      all_membership[cells_m] <- mid
      mc_counts_list[[mid]] <- Matrix::rowSums(counts_mat[, cells_m, drop = FALSE])
      mc_data_list[[mid]]   <- Matrix::rowMeans(data_mat[, cells_m, drop = FALSE])
      mc_meta_list[[mid]]   <- data.frame(metacell_id = mid, cluster = k,
                                          n_cells = length(cells_m),
                                          stringsAsFactors = FALSE)
    }
  }

  mc_counts <- do.call(cbind, mc_counts_list)
  mc_data   <- do.call(cbind, mc_data_list)
  colnames(mc_counts) <- names(mc_counts_list)
  colnames(mc_data)   <- names(mc_data_list)
  mc_meta <- do.call(rbind, mc_meta_list); rownames(mc_meta) <- NULL

  if (isTRUE(verbose)) message(nrow(mc_meta), " metacells from ",
                               length(all_membership), " cells across ",
                               length(clusters_lvl), " clusters.")

  out <- list(counts = mc_counts, data = mc_data, meta = mc_meta,
              membership = all_membership)
  attr(out, "params") <- list(gamma = gamma, k_knn = k_knn,
                              n_var_genes = n_var_genes, min_cells = min_cells,
                              per_cluster = TRUE)
  class(out) <- c("metacells", "list")
  out
}

#' @export
summary.metacells <- function(object, ...) {
  p <- attr(object, "params")
  m <- object$meta
  cat("metacells:", nrow(m), "from", sum(m$n_cells), "cells\n")
  cat("  gamma:", p$gamma, "| per-cluster graining:", p$per_cluster, "\n")
  cat("  clusters:", length(unique(m$cluster)),
      "| metacells/cluster range:",
      min(table(m$cluster)), "-", max(table(m$cluster)), "\n")
  cat("  cells/metacell: median", median(m$n_cells),
      "(", min(m$n_cells), "-", max(m$n_cells), ")\n")
  invisible(object)
}
#' Fine cross-species correspondence by metacell matching
#'
#' Within each matched cluster pair (from \code{\link{match_clusters}}), pairs
#' metacells across species by reciprocal best correlation in gene space over
#' the one-to-one ortholog bridge. Produces a fine grid of matched cross-species
#' metacell pairs — the analog of SAMap's cell-level linking (Tarashansky et al.
#' 2021), but constrained by the interpretable cluster correspondence,
#' coarse-grained for robustness, and computed in gene space rather than a joint
#' embedding.
#'
#' @param metacells_a,metacells_b \code{metacells} objects (\code{\link{build_metacells}})
#'   for species A and B.
#' @param correspondence A \code{cluster_correspondence} (\code{\link{match_clusters}}).
#'   Only reciprocal matches are used.
#' @param homology_graph A \code{homology_graph}; its one-to-one orthologs bridge
#'   the gene axes.
#' @param species_a,species_b Species identifiers matching the graph.
#' @param cor_method \code{"spearman"} (default) or \code{"pearson"}.
#' @param verbose Print progress. Default TRUE.
#' @param mutual_k Neighborhood size for mutual k-nearest-neighbor metacell
#'   matching. Larger values yield more pairs (a denser grid). Default 5.
#'
#' @return A \code{metacell_correspondence}: data.frame of metacell_a,
#'   metacell_b, cluster_a, cluster_b, score. Only reciprocal-best metacell
#'   pairs within matched clusters are retained. Unmatched metacells are dropped
#'   and their count reported in \code{attr(x, "unmatched")}.
#'
#' @export
match_metacells <- function(metacells_a, metacells_b, correspondence,
                            homology_graph,
                            species_a = NULL, species_b = NULL,
                            cor_method = c("spearman", "pearson"),
		            mutual_k = 5,                            
			    verbose = TRUE) {
  cor_method <- match.arg(cor_method)
  stopifnot(inherits(metacells_a, "metacells"),
            inherits(metacells_b, "metacells"))

  # ---- ortholog bridge, oriented so gene_a is species A ---------------------
  hg <- homology_graph
  o2o <- hg[hg$relationship == "ortholog" &
            !is.na(hg$ortholog_type) & hg$ortholog_type == "one2one", ,
            drop = FALSE]
  if (!is.null(species_a) && !is.null(species_b) &&
      all(o2o$species_a == species_b) && all(o2o$species_b == species_a)) {
    o2o <- .flip_edges(o2o)
  }
  ga <- rownames(metacells_a$data); gb <- rownames(metacells_b$data)
  o2o <- o2o[o2o$gene_a %in% ga & o2o$gene_b %in% gb, , drop = FALSE]
  o2o <- o2o[!duplicated(o2o$gene_a) & !duplicated(o2o$gene_b), , drop = FALSE]
  if (nrow(o2o) == 0)
    stop("No one-to-one orthologs bridge the two metacell sets.", call. = FALSE)
  if (isTRUE(verbose)) message("Bridging metacells on ", nrow(o2o), " orthologs.")

  # bridged expression: rows aligned by ortholog
  Xa <- as.matrix(metacells_a$data[o2o$gene_a, , drop = FALSE])  # genes x mcA
  Xb <- as.matrix(metacells_b$data[o2o$gene_b, , drop = FALSE])  # genes x mcB

  meta_a <- metacells_a$meta; meta_b <- metacells_b$meta
  recip <- correspondence[correspondence$reciprocal, , drop = FALSE]
  if (nrow(recip) == 0)
    stop("No reciprocal cluster matches in the correspondence.", call. = FALSE)

  pairs_list <- list()
  n_unmatched_a <- 0L; n_unmatched_b <- 0L

  for (i in seq_len(nrow(recip))) {
    ca <- as.character(recip$cluster_a[i])
    cb <- as.character(recip$cluster_b[i])

    mc_a <- meta_a$metacell_id[as.character(meta_a$cluster) == ca]
    mc_b <- meta_b$metacell_id[as.character(meta_b$cluster) == cb]
    if (length(mc_a) == 0 || length(mc_b) == 0) next

    # correlate metacells of cluster ca (A) against metacells of cluster cb (B)
    cm <- stats::cor(Xa[, mc_a, drop = FALSE],
                     Xb[, mc_b, drop = FALSE], method = cor_method)
    rownames(cm) <- mc_a; colnames(cm) <- mc_b

# mutual k-nearest-neighbors within this cluster pair.
    # For each A-metacell, its top-k B partners; for each B, its top-k A partners.
    # Keep pairs that are mutually in each other's top-k.
    k <- min(mutual_k, ncol(cm), nrow(cm))

    topk_b <- apply(cm, 1, function(r) colnames(cm)[order(r, decreasing = TRUE)[seq_len(k)]])
    topk_a <- apply(cm, 2, function(cl) rownames(cm)[order(cl, decreasing = TRUE)[seq_len(k)]])
    # topk_b: k x nA matrix (B partners per A); topk_a: k x nB (A partners per B)
    topk_b <- if (is.matrix(topk_b)) topk_b else matrix(topk_b, nrow = 1)
    topk_a <- if (is.matrix(topk_a)) topk_a else matrix(topk_a, nrow = 1)
    colnames(topk_b) <- rownames(cm)
    colnames(topk_a) <- colnames(cm)

    matched_a_here <- character(0)
    for (ma in rownames(cm)) {
      b_partners <- topk_b[, ma]
      for (mb in b_partners) {
        if (ma %in% topk_a[, mb]) {                # mutual: each in other's top-k
          pairs_list[[length(pairs_list) + 1L]] <- data.frame(
            metacell_a = ma, metacell_b = mb,
            cluster_a = ca, cluster_b = cb,
            score = cm[ma, mb], stringsAsFactors = FALSE)
          matched_a_here <- c(matched_a_here, ma)
        }
      }
    }
    n_unmatched_a <- n_unmatched_a + length(setdiff(rownames(cm), matched_a_here))
  }

  if (length(pairs_list) == 0)
    stop("No reciprocal metacell pairs found within matched clusters.",
         call. = FALSE)
  out <- do.call(rbind, pairs_list); rownames(out) <- NULL
  out <- out[order(out$cluster_a, -out$score), ]

  if (isTRUE(verbose)) message(nrow(out), " matched metacell pairs across ",
                               nrow(recip), " cluster pairs (",
                               n_unmatched_a, " A-metacells unmatched).")

  attr(out, "unmatched") <- c(a = n_unmatched_a)
  attr(out, "n_orthologs") <- nrow(o2o)
  class(out) <- c("metacell_correspondence", "data.frame")
  out
}

#' @export
summary.metacell_correspondence <- function(object, ...) {
  cat("metacell_correspondence:", nrow(object), "matched pairs\n")
  cat("  ortholog bridge:", attr(object, "n_orthologs"), "genes\n")
  cat("  cluster pairs:", length(unique(paste(object$cluster_a, object$cluster_b))), "\n")
  cat("  score range:", sprintf("%.2f", min(object$score)), "-",
      sprintf("%.2f", max(object$score)), "\n")
  by_cl <- table(paste0(object$cluster_a, "->", object$cluster_b))
  cat("  pairs/cluster range:", min(by_cl), "-", max(by_cl), "\n")
  invisible(object)
}
