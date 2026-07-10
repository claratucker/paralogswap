#' @importFrom rlang .data
NULL

# =============================================================================
# Stage 9 - diagnostic plots
#
# ElbowPlot idiom (Hao et al. 2021): a cheap plot that surfaces the decision
# behind a threshold and hands it back to the user, rather than applying a
# hidden cutoff. All functions return a ggplot object so they compose with
# + theme()/+ labs() like Seurat::ElbowPlot.
# =============================================================================

#' Plot the delta_r distribution as a ranked elbow curve
#'
#' Ranks every scored gene by delta_r = r_paralog - r_ortholog (descending) and
#' plots the sorted curve, with the flagging threshold drawn as a reference
#' line. If delta_r has a natural break between substitution candidates and the
#' smooth bulk, that bend is the candidate set, read off visually; if it is
#' smooth, the plot shows that honestly and marks where the default lands.
#'
#' @param substitutions A \code{substitutions} object
#'   (\code{\link{detect_substitutions}}).
#' @param delta_threshold Threshold line to draw. NULL (default) reads it from
#'   \code{attr(substitutions, "params")$delta_threshold} if present, else 0.3.
#' @param n_label Label the top \code{n_label} genes by delta_r. Default 5;
#'   0 suppresses labels.
#' @return A ggplot object: gene rank (x) vs delta_r (y), threshold line marked.
#'   \code{delta_z} is the Fisher-transformed difference and \code{delta_r} the
#'   raw one; both are always reported, and \code{flagged} follows whichever
#'   \code{delta_scale} names. Rows are sorted by the chosen statistic.
#'   \code{ortholog_undefined} marks rows reached only under
#'   \code{require_ortholog = FALSE}, where the ortholog correlation could not be
#'   computed and zero was substituted for it. For those rows \code{delta_r}
#'   equals \code{r_paralog} and is not a difference between two measured
#'   correlations; read them as the \code{ortholog_invariant} class, where the
#'   substitution is most complete and \code{delta_r} is least meaningful.
#' @export
plot_substitution_elbow <- function(substitutions,
                                    delta_threshold = NULL,
                                    n_label = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_substitution_elbow requires 'ggplot2'.", call. = FALSE)

  df <- as.data.frame(substitutions)
  if (nrow(df) == 0)
    return(ggplot2::ggplot() + ggplot2::theme_bw() +
             ggplot2::labs(title = "No substitutions to plot",
                           x = "gene (ranked by delta_r)",
                           y = expression(Delta * r)))

  if (is.null(delta_threshold)) {
    p <- attr(substitutions, "params")
    delta_threshold <- if (!is.null(p) && !is.null(p$delta_threshold))
      p$delta_threshold else 0.3
  }

  # Plot the statistic that was actually thresholded, and trust the object's
  # `flagged`. Recomputing it here from delta_r silently disagreed with the table
  # whenever detect_substitutions() thresholded on the z scale.
  pp <- attr(substitutions, "params")
  scale_used <- if (!is.null(pp) && !is.null(pp$delta_scale)) pp$delta_scale else "r"
  stat_col <- if (scale_used == "z" && !is.null(df$delta_z)) "delta_z" else "delta_r"
  df$stat <- df[[stat_col]]
  df <- df[order(-df$stat), ]
  df$rank <- seq_len(nrow(df))
  if (is.null(df$flagged)) df$flagged <- df$stat >= delta_threshold
  lab <- df[seq_len(min(n_label, nrow(df))), , drop = FALSE]

  ylab <- if (stat_col == "delta_z")
    expression(Delta * z == atanh(r[paralog]) - atanh(r[ortholog])) else
    expression(Delta * r == r[paralog] - r[ortholog])

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$stat)) +
    ggplot2::geom_hline(yintercept = delta_threshold,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_line(colour = "grey65") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$flagged), size = 1.7) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#c0392b", `FALSE` = "grey55"),
      labels = c(`TRUE` = "flagged", `FALSE` = "not flagged"),
      name = NULL) +
    ggplot2::annotate("text", x = nrow(df), y = delta_threshold,
                      label = sprintf("threshold = %.2f (%s scale)",
                                      delta_threshold, scale_used),
                      hjust = 1, vjust = -0.6, size = 3, colour = "grey40") +
    ggplot2::labs(x = paste0("gene (ranked by ", stat_col, ")"),
                  y = ylab, title = "Substitution elbow") +
    ggplot2::theme_bw()

  if (n_label > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      g <- g + ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(label = .data$gene), size = 3,
        min.segment.length = 0, max.overlaps = Inf)
    } else {
      g <- g + ggplot2::geom_text(
        data = lab, ggplot2::aes(label = .data$gene),
        hjust = -0.15, size = 3)
    }
  }
  g
}

#' Plot ortholog vs. best-paralog correlation, with the diagonal
#'
#' For each focal gene, plots its ortholog correlation (x) against its best
#' paralog correlation (y). The y = x diagonal is the neutral line; a gene's
#' vertical distance above it is delta_r, so a substitution candidate is
#' literally a point sitting above the diagonal. A second, dashed line at
#' y = x + delta_threshold marks the flagging boundary.
#'
#' @param homolog_correlations A \code{homolog_correlations}
#'   (\code{\link{compute_homolog_correlations}}). A \code{substitutions}
#'   object also works (its precomputed r_ortholog / r_paralog are used).
#' @param delta_threshold Boundary line offset above the diagonal. Default 0.3.
#' @param delta_scale Scale on which \code{delta_threshold} is applied. Default
#'   \code{"z"}: Fisher's transform, \code{atanh(r_paralog) - atanh(r_ortholog)}.
#'   Correlation is bounded and its sampling variance depends on its own value,
#'   so a raw difference of 0.3 near r = 0 is not the same quantity as one near
#'   r = 0.9; the transform makes the two comparable. \code{"r"} thresholds the
#'   untransformed difference and reproduces the behaviour of versions before
#'   this argument existed.
#'
#' @param show_invariant Draw genes whose ortholog correlation is undefined in a
#'   strip left of the panel. Default TRUE. These are the \code{ortholog_invariant}
#'   class -- the completest substitutions, which no delta can score. Omitting
#'   them silently drops the genes that matter most.
#' @param show_r_boundary When \code{delta_scale = "z"}, also draw the raw-scale
#'   boundary \code{y = x + delta_threshold} as a dotted line. It exits the unit
#'   square at \code{x = 1 - delta_threshold}. Default TRUE.
#' @param label_genes Genes to label (e.g. \code{"RAMP1"}). NULL = auto-label
#'   the strongest candidates; FALSE = no labels.
#' @param max_auto_labels Cap on auto-labelled candidates. Default 10.
#'
#' @return A ggplot object: r_ortholog (x) vs max r_paralog (y).
#' @export
plot_ortholog_vs_paralog <- function(homolog_correlations,
                                     delta_threshold = NULL,
                                     delta_scale = NULL,
                                     show_invariant = TRUE,
                                     show_r_boundary = TRUE,
                                     label_genes = NULL,
                                     max_auto_labels = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_ortholog_vs_paralog requires 'ggplot2'.", call. = FALSE)

  hc <- as.data.frame(homolog_correlations)

  # Read the threshold and scale the object was actually built with. A bare
  # default here would draw a boundary the calls were never made against.
  pp <- attr(homolog_correlations, "params")
  if (is.null(delta_threshold))
    delta_threshold <- if (!is.null(pp$delta_threshold)) pp$delta_threshold else 0.3
  if (is.null(delta_scale))
    delta_scale <- if (!is.null(pp$delta_scale)) pp$delta_scale else "r"

  # Accept a precomputed substitutions object directly.
  if (all(c("r_ortholog", "r_paralog", "gene") %in% names(hc))) {
    df <- data.frame(gene = hc$gene, r_ortholog = hc$r_ortholog,
                     r_paralog = hc$r_paralog, stringsAsFactors = FALSE)
    df <- df[is.finite(df$r_paralog), , drop = FALSE]
  } else {
    # Collapse homolog_correlations to one row per focal gene.
    genes <- unique(hc$gene_a)
    recs <- lapply(genes, function(gene) {
      sub  <- hc[hc$gene_a == gene, ]
      orth <- sub[sub$relationship == "ortholog" & is.finite(sub$r), ]
      para <- sub[sub$relationship == "paralog"  & is.finite(sub$r), ]
      # A gene whose ortholog correlation is undefined but whose paralogs are
      # finite is the ortholog_invariant class: the completest substitution, and
      # the one delta cannot score. Retain it with r_ortholog = NA; it is drawn
      # in the strip, not the scatter.
      if (nrow(para) == 0) return(NULL)
      best <- para[which.max(para$r), ]
      data.frame(gene = gene,
                 r_ortholog = if (nrow(orth)) orth$r[1] else NA_real_,
                 r_paralog = best$r, stringsAsFactors = FALSE)
    })
    df <- do.call(rbind, recs)
  }

  if (is.null(df) || nrow(df) == 0)
    return(ggplot2::ggplot() + ggplot2::theme_bw() +
             ggplot2::labs(title = "No gene has both an ortholog and a paralog correlation"))

  inv <- df[!is.finite(df$r_ortholog), , drop = FALSE]   # ortholog_invariant
  df  <- df[is.finite(df$r_ortholog), , drop = FALSE]
  if (nrow(df) == 0)
    return(ggplot2::ggplot() + ggplot2::theme_bw() +
             ggplot2::labs(title = "No gene has a measurable ortholog correlation"))

  df$delta_r <- df$r_paralog - df$r_ortholog
  df$delta_z <- .fisher_z(df$r_paralog) - .fisher_z(df$r_ortholog)
  df$candidate <- if (delta_scale == "z") df$delta_z >= delta_threshold else
                                          df$delta_r >= delta_threshold
  if (!is.null(hc$flagged) && !is.null(hc$gene))
    df$candidate <- hc$flagged[match(df$gene, hc$gene)]  # genes in `inv` are absent from df; match() drops them

  stat_col <- if (delta_scale == "z") "delta_z" else "delta_r"
  rng <- range(c(df$r_ortholog, df$r_paralog), na.rm = TRUE)
  .pad <- 0.12 * diff(rng)
  .strip_x <- rng[1] - .pad / 2
  .xlim <- c(rng[1] - .pad,
             rng[2])   # rng spans BOTH axes: clipping to max(r_paralog) hides
                       # genes whose ortholog correlation exceeds every paralog

  if (is.null(label_genes)) {
    lab <- df[df$candidate, , drop = FALSE]
    lab <- lab[order(-lab[[stat_col]]), , drop = FALSE]
    if (nrow(lab) > max_auto_labels) lab <- lab[seq_len(max_auto_labels), , drop = FALSE]
  } else if (isFALSE(label_genes)) {
    lab <- df[0, , drop = FALSE]
  } else {
    lab <- df[df$gene %in% label_genes, , drop = FALSE]
  }

  g <- ggplot2::ggplot(df, ggplot2::aes(.data$r_ortholog, .data$r_paralog)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40") +
    # Flagging boundary. On the z scale it is a curve, y = tanh(atanh(x) + t),
    # which stays inside the unit square. The straight r-scale line y = x + t is
    # drawn faintly for comparison: it leaves the square at x = 1 - t, above
    # which no gene can be flagged however well its paralog tracks.
    ggplot2::geom_line(
      data = local({
        xx <- seq(max(-0.999, rng[1]), min(0.999, rng[2]), length.out = 400)
        yy <- if (delta_scale == "z") tanh(atanh(xx) + delta_threshold) else
                                      xx + delta_threshold
        data.frame(r_ortholog = xx, r_paralog = yy)[yy <= 1, ]
      }),
      linetype = "dashed", colour = "grey55") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$candidate),
                        size = 1.7, alpha = 0.85) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#c0392b", `FALSE` = "grey55"),
      labels = c(`TRUE` = "candidate", `FALSE` = "conserved"), name = NULL) +
    # coord_fixed, not coord_cartesian + aspect.ratio: the latter squares the
    # panel, not the data units, so the y = x diagonal is drawn off 45 degrees
    # whenever the reserved strip makes the x-range wider than the y-range.
    ggplot2::coord_fixed(ratio = 1, xlim = .xlim, ylim = rng) +

    ggplot2::labs(
      x = expression(r[ortholog]),
      y = expression(max ~ r[paralog]),
      title = "Ortholog vs. best paralog correlation",
      subtitle = sprintf("%d candidate%s; %d gene%s with undefined ortholog correlation%s",
                         sum(df$candidate, na.rm = TRUE),
                         if (sum(df$candidate, na.rm = TRUE) == 1) "" else "s",
                         nrow(inv), if (nrow(inv) == 1) "" else "s",
                         if (!isTRUE(show_invariant) && nrow(inv) > 0) " (hidden)" else "")) +
    ggplot2::theme_bw()

  # ---- ortholog_invariant genes: a strip, not a coordinate -----------------
  # These have no r_ortholog. Placing them at x = 0 would assert exactly what
  # detect_substitutions() marks as an assertion -- that a silent ortholog is
  # uncorrelated with the focal gene. They belong outside the panel.
  if (isTRUE(show_invariant) && nrow(inv) > 0) {
    inv$r_ortholog <- .strip_x
    inv$candidate <- TRUE
    g <- g +
      ggplot2::geom_vline(xintercept = rng[1] - .pad * 0.10,
                          colour = "grey75", linewidth = 0.3) +
      ggplot2::geom_point(data = inv, shape = 21, fill = NA,
                          colour = "#c0392b", size = 1.9, stroke = 0.6) +
      ggplot2::annotate("text", x = .strip_x, y = rng[2],
                        label = "ortholog\nundefined", size = 2.6,
                        colour = "grey40", vjust = 1, lineheight = 0.9)
  }

  # The raw-scale boundary, for comparison. It leaves the unit square at
  # x = 1 - delta_threshold: above that, no gene can be flagged on the r scale
  # however well its paralog tracks. Genes between the two lines are exactly
  # those the Fisher transform recovers.
  if (isTRUE(show_r_boundary) && delta_scale == "z") {
    g <- g + ggplot2::geom_abline(slope = 1, intercept = delta_threshold,
                                  linetype = "dotted", colour = "grey80")
  }

  if (nrow(lab) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      g <- g + ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(label = .data$gene), size = 3,
        min.segment.length = 0, max.overlaps = Inf)
    } else {
      g <- g + ggplot2::geom_text(
        data = lab, ggplot2::aes(label = .data$gene),
        hjust = -0.15, size = 3)
    }
  }
  g
}

#' Mean within-metacell expression coherence
#'
#' For each metacell, the mean pairwise correlation (in single-cell expression
#' space) among the cells assigned to it. A continuous coherence measure: high
#' when a metacell groups genuinely similar cells, lower when graining lumps in
#' heterogeneity. Cells-per-metacell are subsampled to \code{max_cells} to keep
#' the pairwise correlation cheap at coarse graining.
#'
#' @param data_mat genes x cells expression matrix (single-cell, not aggregated).
#' @param membership Named vector mapping each cell to its metacell id.
#' @param cor_method \code{"pearson"} (default) or \code{"spearman"}.
#' @param max_cells Cap on cells sampled per metacell for the calc. Default 30.
#' @param seed RNG seed for the subsample. Default 1.
#' @return Named numeric vector: metacell id -> mean within-metacell correlation
#'   (NA for single-cell metacells).
#' @keywords internal
.within_metacell_cor <- function(data_mat, membership,
                                 cor_method = "pearson",
                                 max_cells = 30, seed = 1) {
  mcs <- unique(membership)
  set.seed(seed)
  out <- vapply(mcs, function(m) {
    cells <- names(membership)[membership == m]
    if (length(cells) < 2) return(NA_real_)
    if (length(cells) > max_cells) cells <- sample(cells, max_cells)
    X <- as.matrix(data_mat[, cells, drop = FALSE])
    keep <- apply(X, 1, stats::sd) > 0          # drop zero-variance genes
    if (sum(keep) < 2) return(NA_real_)
    cm <- stats::cor(X[keep, , drop = FALSE], method = cor_method)
    mean(cm[upper.tri(cm)])
  }, numeric(1))
  names(out) <- mcs
  out
}

#' Plot metacell fidelity across SuperCell graining levels
#'
#' Sweeps gamma and plots two single-species fidelity quantities that stay flat
#' while coarse-graining is safe and degrade once metacells begin blurring real
#' structure: within-metacell \emph{coherence} (mean pairwise correlation among
#' a metacell's cells) and correlation \emph{preservation} (how faithfully the
#' gene-gene correlation structure at a given gamma reproduces the finest
#' graining). Choose the largest gamma still on the plateau. Read gamma off
#' fidelity, never off substitution count (see the anti-circularity note in the
#' validation vignette).
#'
#' @param clusters A clustered Seurat object or \code{\link{cluster_species}}
#'   result (the same input \code{\link{build_metacells}} takes).
#' @param gammas Graining levels to sweep. Default c(5, 10, 20, 50, 100). The
#'   smallest is the reference for preservation.
#' @param metric \code{"both"} (default), \code{"coherence"}, or
#'   \code{"preservation"}.
#' @param assay,slot Expression to grain. Default RNA/data (matches build_metacells).
#' @param cor_method \code{"pearson"} (default) or \code{"spearman"}.
#' @param n_panel Number of gene pairs (drawn from top-variance genes) used for
#'   the preservation panel. Default 500.
#' @param seed RNG seed. Default 1.
#' @param verbose Print progress. Default TRUE.
#'
#' @return A ggplot of the fidelity curve. The underlying per-gamma values
#'   (\code{gamma}, \code{coherence}, \code{preservation}) are attached as
#'   \code{attr(x, "curve")} so the chosen gamma can be justified from the
#'   record rather than from the rendered figure. The function selects no
#'   gamma; read the largest value still on the preservation plateau and pass
#'   it to \code{\link{build_metacells}}.
#' @export
plot_graining_curve <- function(clusters,
                                gammas = c(5, 10, 20, 50, 100),
                                metric = c("both", "coherence", "preservation"),
                                assay = "RNA", slot = "data",
                                cor_method = "pearson",
                                n_panel = 500, seed = 1,
                                verbose = TRUE) {
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_graining_curve requires 'ggplot2'.", call. = FALSE)

  obj <- .unwrap_clustered(clusters)   # shared helper (correspondence.R)
  data_mat <- Seurat::GetAssayData(obj, assay = assay, layer = slot)
  gammas <- sort(unique(gammas))

  rows <- list(); grids <- list()
  for (g in gammas) {
    if (isTRUE(verbose)) message("  graining at gamma ", g, " ...")
    mc <- build_metacells(clusters, gamma = g, assay = assay, verbose = FALSE)
    coh <- .within_metacell_cor(data_mat, mc$membership,
                               cor_method = cor_method, seed = seed)
    rows[[as.character(g)]] <- data.frame(
      gamma = g, coherence = mean(coh, na.rm = TRUE),
      n_metacells = ncol(mc$data), stringsAsFactors = FALSE)
    grids[[as.character(g)]] <- mc$data
  }
  curve <- do.call(rbind, rows)

  if (metric %in% c("both", "preservation")) {
    finest <- as.matrix(grids[[as.character(min(gammas))]])
    v <- apply(finest, 1, stats::var)
    panel <- names(sort(v, decreasing = TRUE))[seq_len(min(n_panel * 2, length(v)))]
    set.seed(seed); panel <- sample(panel)
    npair <- floor(length(panel) / 2)
    pairs <- cbind(panel[seq_len(npair)], panel[npair + seq_len(npair)])

    panel_corr <- function(grid) {
      G <- as.matrix(grid)
      vapply(seq_len(nrow(pairs)), function(i) {
        a <- pairs[i, 1]; b <- pairs[i, 2]
        if (!(a %in% rownames(G)) || !(b %in% rownames(G))) return(NA_real_)
        xa <- G[a, ]; xb <- G[b, ]
        if (stats::sd(xa) == 0 || stats::sd(xb) == 0) return(NA_real_)
        stats::cor(xa, xb, method = cor_method)
      }, numeric(1))
    }
    ref <- panel_corr(finest)
    curve$preservation <- vapply(gammas, function(g) {
      cg <- panel_corr(grids[[as.character(g)]])
      ok <- is.finite(ref) & is.finite(cg)
      if (sum(ok) < 3) return(NA_real_)
      stats::cor(ref[ok], cg[ok], method = "pearson")
    }, numeric(1))
  }

  keep <- switch(metric,
                 both         = c("coherence", "preservation"),
                 coherence    = "coherence",
                 preservation = "preservation")
  long <- do.call(rbind, lapply(keep, function(q) data.frame(
    gamma = curve$gamma, quantity = q, value = curve[[q]],
    stringsAsFactors = FALSE)))

p <- ggplot2::ggplot(long, ggplot2::aes(.data$gamma, .data$value,
                                     colour = .data$quantity)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::scale_x_log10(breaks = gammas) +
    ggplot2::scale_colour_manual(
      values = c(coherence = "#2c7fb8", preservation = "#c0392b"),
      labels = c(coherence = "within-metacell coherence",
                 preservation = "correlation preservation"),
      name = NULL) + 
    ggplot2::labs(x = "gamma (cells per metacell, log scale)", y = "fidelity",
                  title = "Graining fidelity curve",
                  subtitle = "choose the largest gamma still on the plateau") +
    ggplot2::theme_bw()

  attr(p, "curve") <- curve
  p
}

#' Heatmap of cross-species cluster correspondence
#'
#' Renders the full cluster-by-cluster correlation matrix computed by
#' \code{\link{match_clusters}} and carried on the result as
#' \code{attr(x, "cor_matrix")}, outlining reciprocal best matches. This is the
#' visual form of the claim that cross-species correspondence can be read from
#' expression without cell-type labels.
#'
#' Every cell carries a score. \code{match_clusters} correlates every A-cluster
#' against every B-cluster and uses \code{min_correspondence} only to gate the
#' \code{reciprocal} flag, so no pair is discarded and no cell is empty. Read a
#' pale cell as a low correlation, not as missing data.
#'
#' With \code{cluster_order = TRUE} rows are permuted so that each A-cluster sits
#' near the B-cluster it matches best. Block structure that survives this
#' reordering is real; vertical stripes mean many A-clusters are attracted to one
#' B-cluster, which reordering cannot hide.
#'
#' @param correspondence A \code{cluster_correspondence} from \code{match_clusters}.
#' @param mark_reciprocal Outline reciprocal best matches. Default TRUE.
#' @param low,high Fill scale endpoints.
#' @param show_values Print the score in each tile. Default FALSE.
#' @param cluster_order Reorder rows by best match. Default TRUE.
#'
#' @return A ggplot. Composable with \code{+ theme()} as usual.
#' @export
plot_correspondence_heatmap <- function(correspondence,
                                        mark_reciprocal = TRUE,
                                        low = "#f7fbff", high = "#08306b",
                                        show_values = FALSE,
                                        cluster_order = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_correspondence_heatmap requires the 'ggplot2' package.", call. = FALSE)

  cormat <- attr(correspondence, "cor_matrix")
  if (is.null(cormat))
    stop("No 'cor_matrix' attribute; this object predates it. Re-run match_clusters().",
         call. = FALSE)

  recip <- as.data.frame(correspondence)
  recip <- recip[which(recip$reciprocal), , drop = FALSE]

  # numeric cluster labels sort numerically, not lexically
  .ord <- function(x) {
    n <- suppressWarnings(as.numeric(x))
    if (anyNA(n)) order(x) else order(n)
  }
  col_lv <- colnames(cormat)[.ord(colnames(cormat))]

  if (isTRUE(cluster_order)) {
    best <- match(colnames(cormat)[max.col(cormat, ties.method = "first")], col_lv)
    row_lv <- rownames(cormat)[order(best, -apply(cormat, 1, max))]
  } else {
    row_lv <- rownames(cormat)[.ord(rownames(cormat))]
  }

  d <- expand.grid(cluster_a = rownames(cormat), cluster_b = colnames(cormat),
                   stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  d$score <- as.vector(cormat)   # both column-major; rows vary fastest
  d$cluster_a <- factor(d$cluster_a, levels = rev(row_lv))
  d$cluster_b <- factor(d$cluster_b, levels = col_lv)

  bridge <- .bridge_n(correspondence)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$cluster_b, y = .data$cluster_a)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$score), colour = "grey92") +
    ggplot2::scale_fill_gradient(low = low, high = high, name = "score") +
    ggplot2::labs(
      x = "species B cluster", y = "species A cluster",
      title = "Cross-species cluster correspondence",
      subtitle = paste0(nrow(cormat), " x ", ncol(cormat), " clusters, ",
                        nrow(recip), " reciprocal best matches",
                        if (!is.null(bridge)) paste0("; bridge of ", bridge,
                                                     " one2one orthologs") else "")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  if (isTRUE(mark_reciprocal) && nrow(recip) > 0) {
    recip$cluster_a <- factor(recip$cluster_a, levels = rev(row_lv))
    recip$cluster_b <- factor(recip$cluster_b, levels = col_lv)
    p <- p + ggplot2::geom_tile(data = recip, fill = NA,
                                colour = "firebrick", linewidth = 0.7)
  }
  if (isTRUE(show_values)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$score)),
                                size = 2.2, colour = "grey20")
  }
  p
}
