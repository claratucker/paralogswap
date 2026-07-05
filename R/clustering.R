#' Independent per-species clustering
#'
#' Clusters each species independently using the standard Seurat workflow
#' (Stuart et al. 2019; Hao et al. 2021): \code{NormalizeData} ->
#' \code{FindVariableFeatures} -> \code{ScaleData} -> \code{RunPCA} ->
#' \code{FindNeighbors} -> \code{FindClusters}. No cross-species information is
#' used at this stage; cross-species correspondence is derived later from
#' expression (\code{\link{match_clusters}}).
#'
#' @param object Either a single Seurat object containing all species (split by
#'   \code{species_col}), or a named list of per-species Seurat objects keyed by
#'   species identifier.
#' @param species_col Column in the object's metadata giving species identity.
#'   Used only when \code{object} is a single merged object; ignored for a list.
#' @param assay Input counts assay. With \code{normalization = "sctransform"},
#'   clustering is performed on the resulting \code{SCT} assay.
#' @param normalization \code{"lognorm"} (default) or \code{"sctransform"}.
#' @param resolution Clustering resolution. A scalar (default 0.8, the Seurat
#'   convention), or a vector to produce multiple clusterings for the Stage 7
#'   stability check.
#' @param n_pcs Number of principal components carried into neighbor-graph
#'   construction and clustering. Inspect with \code{plot_pca_elbow}.
#' @param n_var_genes Number of variable features per species.
#' @param verbose Print per-species progress. Default \code{TRUE}.
#'
#' @return A named list of the per-species Seurat objects, each carrying its
#'   PCA reduction, neighbor graph, and cluster assignments. For a scalar
#'   \code{resolution}, assignments are in \code{seurat_clusters}; for a vector,
#'   in per-resolution columns \code{clusters_res<r>}. Clustering parameters are
#'   recorded in \code{attr(x, "params")}.
#'
#' @export
cluster_species <- function(object,
                            species_col   = "species",
                            assay         = "RNA",
                            normalization = c("lognorm", "sctransform"),
                            resolution    = 0.8,
                            n_pcs         = 30,
                            n_var_genes   = 2000,
                            verbose       = TRUE) {
  normalization <- match.arg(normalization)

  # ---- resolve input to a named list of per-species objects ------------------
  obj_list <- .as_species_list(object, species_col)

  # ---- cluster each species independently ------------------------------------
  clustered <- lapply(names(obj_list), function(sp) {
    if (isTRUE(verbose)) message("Clustering species: ", sp,
                                 " (", ncol(obj_list[[sp]]), " cells)")
    .cluster_one(
      obj_list[[sp]],
      assay = assay, normalization = normalization,
      resolution = resolution, n_pcs = n_pcs,
      n_var_genes = n_var_genes, verbose = verbose
    )
  })
  names(clustered) <- names(obj_list)

  attr(clustered, "params") <- list(
    normalization = normalization, resolution = resolution,
    n_pcs = n_pcs, n_var_genes = n_var_genes, assay = assay
  )
  class(clustered) <- c("species_clusters", "list")
  clustered
}

# ---- input normalization: single merged object OR named list -> named list ---
.as_species_list <- function(object, species_col) {
  if (is.list(object) && !inherits(object, "Seurat")) {
    if (is.null(names(object)) || any(names(object) == "")) {
      stop("When `object` is a list, it must be NAMED (one name per species).",
           call. = FALSE)
    }
    if (!all(vapply(object, inherits, logical(1), "Seurat"))) {
      stop("All elements of the `object` list must be Seurat objects.",
           call. = FALSE)
    }
    return(object)
  }
  if (inherits(object, "Seurat")) {
    if (!species_col %in% colnames(object[[]])) {
      stop("`species_col` = '", species_col,
           "' not found in object metadata. Provide a named list instead, ",
           "or set species_col to an existing column.", call. = FALSE)
    }
    sp <- as.character(object[[species_col, drop = TRUE]])
    parts <- split(colnames(object), sp)
    return(lapply(parts, function(cells) subset(object, cells = cells)))
  }
  stop("`object` must be a Seurat object or a named list of Seurat objects.",
       call. = FALSE)
}

# ---- cluster a single Seurat object ----------------------------------------
.cluster_one <- function(obj, assay, normalization, resolution, n_pcs,
                         n_var_genes, verbose) {
  # cap n_pcs and n_var_genes for small objects (avoids irlba/PCA failures)
  n_cells <- ncol(obj)
  n_feat  <- nrow(obj)
  n_pcs_use <- min(n_pcs, n_cells - 1L, n_feat - 1L)
  if (n_pcs_use < n_pcs && isTRUE(verbose)) {
    message("  n_pcs reduced to ", n_pcs_use, " (small object)")
  }

  if (normalization == "sctransform") {
    obj <- Seurat::SCTransform(obj, assay = assay, verbose = FALSE,
                               variable.features.n = n_var_genes)
    # SCTransform sets the default assay to "SCT" and scales internally
  } else {
    Seurat::DefaultAssay(obj) <- assay
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_var_genes,
                                        verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
  }

  obj <- Seurat::RunPCA(obj, npcs = n_pcs_use, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = 1:n_pcs_use, verbose = FALSE)

  # scalar or vector resolution
  if (length(resolution) == 1L) {
    obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)
    # seurat_clusters is set automatically
  } else {
    for (r in resolution) {
      obj <- Seurat::FindClusters(obj, resolution = r, verbose = FALSE)
      obj[[paste0("clusters_res", r)]] <- Seurat::Idents(obj)
    }
    # leave Idents at the last resolution; per-res columns hold all of them
  }
  obj
}

#' @export
summary.species_clusters <- function(object, ...) {
  p <- attr(object, "params")
  cat("species_clusters:", length(object), "species\n")
  cat("  normalization:", p$normalization, "| n_pcs:", p$n_pcs,
      "| n_var_genes:", p$n_var_genes, "\n")
  cat("  resolution(s):", paste(p$resolution, collapse = ", "), "\n\n")
  for (sp in names(object)) {
    o <- object[[sp]]
    if (length(p$resolution) == 1L) {
      nclust <- length(unique(Seurat::Idents(o)))
      cat(sprintf("  %-12s %6d cells, %3d clusters\n", sp, ncol(o), nclust))
    } else {
      counts <- vapply(p$resolution, function(r)
        length(unique(o[[paste0("clusters_res", r)]][, 1])), integer(1))
      cat(sprintf("  %-12s %6d cells, clusters: %s\n", sp, ncol(o),
                  paste(sprintf("res%.1f=%d", p$resolution, counts),
                        collapse = ", ")))
    }
  }
  invisible(object)
}
