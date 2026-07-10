#' Assess robustness of substitution calls to clustering and graining
#'
#' Re-runs the detection pipeline across a grid of metacell graining factors
#' (\code{gammas}) and, optionally, clustering resolutions, returning one row per
#' (resolution, gamma, focal gene) so that \code{delta_r} can be read against the
#' grid directly.
#'
#' Two modes. With \code{resolutions = NULL} (default) the supplied clustered
#' objects are taken as given and only Stages 3-6 are re-run: correspondence is
#' computed once, since it depends on cluster profiles and not on graining, then
#' metacells are rebuilt at each gamma. This is the cheap, testable path.
#' Supplying \code{resolutions} additionally requires \code{seurat_a} and
#' \code{seurat_b} and re-runs \code{\link{cluster_species}} at each resolution,
#' which is expensive: the cost is \code{length(resolutions) * length(gammas)}
#' pipelines.
#'
#' A call that survives the grid is a claim about the data. A call that appears
#' only at one gamma is a claim about the graining.
#'
#' @param clusters_a,clusters_b Clustered Seurat objects. Ignored, and may be
#'   NULL, when \code{resolutions} is supplied.
#' @param homology_graph A \code{homology_graph}.
#' @param gammas Numeric vector of graining factors. Default \code{c(10, 20, 30)}.
#' @param resolutions Numeric vector of Louvain resolutions, or NULL (default)
#'   to use the clustering already present in \code{clusters_a}/\code{clusters_b}.
#' @param seurat_a,seurat_b Raw Seurat objects; required iff \code{resolutions}
#'   is supplied.
#' @param focal_genes Passed to \code{\link{compute_homolog_correlations}}.
#' @param species_a,species_b Passed through to the staged functions.
#' @param delta_threshold Passed to \code{\link{detect_substitutions}}.
#' @param min_correspondence Passed to \code{\link{match_clusters}}.
#' @param verbose Print progress. Default TRUE.
#' @param ... Further arguments to \code{\link{build_metacells}}.
#'
#' @return A \code{stability_grid}: a data.frame carrying every column returned
#'   by \code{\link{detect_substitutions}}, prefixed by \code{resolution} and
#'   \code{gamma}, with \code{n_metacells_a}, \code{n_metacells_b} and
#'   \code{n_metacell_pairs} recording the size of the run that produced it.
#'   Grid cells that error contribute zero rows and emit a warning; the failure
#'   is recorded in \code{attr(x, "failures")}.
#'
#' @seealso \code{\link{plot_graining_curve}}
#' @export
assess_stability <- function(clusters_a = NULL, clusters_b = NULL,
                             homology_graph,
                             gammas = c(10, 20, 30),
                             resolutions = NULL,
                             seurat_a = NULL, seurat_b = NULL,
                             focal_genes = NULL,
                             species_a = NULL, species_b = NULL,
                             delta_threshold = 0.3,
                             min_correspondence = 0.3,
                             verbose = TRUE,
                             ...) {
  recluster <- !is.null(resolutions)
  if (recluster) {
    if (is.null(seurat_a) || is.null(seurat_b))
      stop("resolutions supplied, so seurat_a and seurat_b are required.", call. = FALSE)
  } else {
    if (is.null(clusters_a) || is.null(clusters_b))
      stop("Supply clusters_a and clusters_b, or resolutions with seurat_a/seurat_b.",
           call. = FALSE)
    resolutions <- NA_real_   # single pass, resolution unknown/inherited
  }
  if (length(gammas) == 0L) stop("gammas must be non-empty.", call. = FALSE)

  out <- list(); failures <- list()

  for (res in resolutions) {
    if (recluster) {
      if (isTRUE(verbose)) message("resolution ", res, ": clustering both species")
      ca <- cluster_species(seurat_a, resolution = res, verbose = FALSE)
      cb <- cluster_species(seurat_b, resolution = res, verbose = FALSE)
    } else {
      ca <- clusters_a; cb <- clusters_b
    }

    # Stage 3 is gamma-invariant: cluster profiles precede metacells.
    corr <- match_clusters(ca, cb, homology_graph,
                           species_a = species_a, species_b = species_b,
                           min_correspondence = min_correspondence,
                           verbose = FALSE)

    for (g in gammas) {
      tag <- paste0("resolution=", res, ", gamma=", g)
      if (isTRUE(verbose)) message("  ", tag)
      cell <- tryCatch({
        mc_a <- build_metacells(ca, gamma = g, verbose = FALSE, ...)
        mc_b <- build_metacells(cb, gamma = g, verbose = FALSE, ...)
        mm   <- match_metacells(mc_a, mc_b, corr, homology_graph,
                                species_a = species_a, species_b = species_b,
                                verbose = FALSE)
        hc   <- compute_homolog_correlations(mc_a, mc_b, mm, homology_graph,
                                             species_a = species_a,
                                             species_b = species_b,
                                             focal_genes = focal_genes,
                                             verbose = FALSE)
        sub  <- as.data.frame(detect_substitutions(hc, delta_threshold = delta_threshold))
        if (nrow(sub) == 0L) NULL else cbind(
          resolution = res, gamma = g,
          n_metacells_a = .n_metacells(mc_a), n_metacells_b = .n_metacells(mc_b),
          n_metacell_pairs = nrow(as.data.frame(mm)),
          sub, stringsAsFactors = FALSE)
      }, error = function(e) {
        warning("Grid cell failed (", tag, "): ", conditionMessage(e), call. = FALSE)
        failures[[tag]] <<- conditionMessage(e)
        NULL
      })
      if (!is.null(cell)) out[[length(out) + 1L]] <- cell
    }
  }

  if (length(out) == 0L)
    stop("Every grid cell failed or returned no rows. See warnings.", call. = FALSE)

  res_df <- do.call(rbind, out)
  rownames(res_df) <- NULL
  attr(res_df, "gammas") <- gammas
  attr(res_df, "resolutions") <- resolutions
  attr(res_df, "failures") <- failures
  class(res_df) <- c("stability_grid", "data.frame")
  res_df
}

# Metacell count, without assuming the container's class
.n_metacells <- function(mc) {
  if (!is.null(mc$membership)) return(length(unique(mc$membership)))
  if (!is.null(mc$SC$membership)) return(length(unique(mc$SC$membership)))
  if (inherits(mc, "Seurat")) return(ncol(mc))
  NA_integer_
}

#' @export
summary.stability_grid <- function(object, ...) {
  d <- as.data.frame(object)
  n_cells <- length(unique(paste(d$resolution, d$gamma)))
  cat("stability_grid:\n")
  cat("  grid cells with calls:", n_cells, "\n")
  cat("  gammas:", paste(attr(object, "gammas"), collapse = ", "), "\n")
  if (!all(is.na(attr(object, "resolutions"))))
    cat("  resolutions:", paste(attr(object, "resolutions"), collapse = ", "), "\n")
  f <- attr(object, "failures")
  if (length(f)) cat("  failed cells:", length(f), "\n")
  gcol <- intersect(c("focal_gene", "gene"), names(d))[1]
  if (!is.na(gcol)) {
    tab <- sort(table(d[[gcol]]), decreasing = TRUE)
    cat("  genes called in all", n_cells, "cells:",
        paste(names(tab)[tab == n_cells], collapse = ", "), "\n")
  }
  invisible(object)
}
