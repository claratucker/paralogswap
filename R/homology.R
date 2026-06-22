#' Homology graph objects
#'
#' Construct and inspect a \code{homology_graph}: a tidy edge table describing
#' ortholog and paralog relationships between genes of two species. Most
#' paralogswap functions take a \code{homology_graph} as input.
#'
#' A \code{homology_graph} is a \code{data.frame} (with class
#' \code{"homology_graph"}) of one row per homologous gene pair, with columns:
#' \describe{
#'   \item{gene_a}{character. Gene identifier in species A (the focal species;
#'     edges are directed A to B).}
#'   \item{gene_b}{character. Gene identifier in species B.}
#'   \item{relationship}{factor with levels \code{"ortholog"},
#'     \code{"paralog"}, and \code{"unclassified"}. Unclassified pairs are
#'     sequence-similar homologs that are neither assigned ortholog nor
#'     paralog; they are carried through the pipeline but not scored for
#'     substitution (Stage 5).}
#'   \item{duplication_level}{character or \code{NA}. Taxonomic level of the most
#'     recent common ancestor at which the pair shares an orthology group;
#'     used to stratify paralogs by duplication age.}
#'   \item{seq_similarity}{numeric in \code{[0, 1]} or \code{NA}. Optional
#'     sequence-similarity weight for the edge.}
#' }
#'
#' Edges are directed from species A to species B: every row represents a focal
#' gene in A and a homologous partner in B. To assess substitution in both
#' lineages, score in both directions via the \code{symmetric} argument of
#' \code{\link{detect_substitutions}} (or run the pipeline twice with
#' \code{species_a} and \code{species_b} swapped).
#'
#' @param edge_df A \code{data.frame} containing at least \code{gene_a},
#'   \code{gene_b}, and \code{relationship}. \code{duplication_level} and
#'   \code{seq_similarity} are added as \code{NA} if absent.
#' @param verbose Print a brief construction message. Default \code{TRUE}.
#'
#' @return An object of class \code{homology_graph}.
#'
#' @examples
#' edges <- data.frame(
#'   gene_a = c("Actb", "Gapdh", "Myl2"),
#'   gene_b = c("ACTB", "GAPDH", "MYL7"),
#'   relationship = c("ortholog", "ortholog", "paralog")
#' )
#' hg <- as_homology_graph(edges)
#' summary(hg)
#'
#' @export
as_homology_graph <- function(edge_df, verbose = TRUE) {
  if (inherits(edge_df, "homology_graph")) {
    return(edge_df)
  }
  if (!is.data.frame(edge_df)) {
    stop("`edge_df` must be a data.frame.", call. = FALSE)
  }

  required <- c("gene_a", "gene_b", "relationship")
  missing_cols <- setdiff(required, colnames(edge_df))
  if (length(missing_cols) > 0) {
    stop(
      "`edge_df` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Coerce gene IDs to character to guard against factor inputs.
  edge_df$gene_a <- as.character(edge_df$gene_a)
  edge_df$gene_b <- as.character(edge_df$gene_b)

  # relationship must be one of the three recognized levels.
  rel <- as.character(edge_df$relationship)
  valid_rel <- c("ortholog", "paralog", "unclassified")
  bad_rel <- setdiff(unique(rel), valid_rel)
  if (length(bad_rel) > 0) {
    stop(
      "`relationship` must contain only 'ortholog', 'paralog', or ",
      "'unclassified'; found: ",
      paste(bad_rel, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  edge_df$relationship <- factor(rel, levels = valid_rel)

  # Fill optional columns with NA if absent.
  if (!"duplication_level" %in% colnames(edge_df)) {
    edge_df$duplication_level <- rep(NA_character_, nrow(edge_df))
  } else {
    edge_df$duplication_level <- as.character(edge_df$duplication_level)
  }
  if (!"seq_similarity" %in% colnames(edge_df)) {
    edge_df$seq_similarity <- rep(NA_real_, nrow(edge_df))
  } else {
    edge_df$seq_similarity <- as.numeric(edge_df$seq_similarity)
    if (any(edge_df$seq_similarity < 0 | edge_df$seq_similarity > 1,
            na.rm = TRUE)) {
      stop("`seq_similarity` must lie in [0, 1].", call. = FALSE)
    }
  }

  # Drop rows with missing gene IDs; these cannot participate in correlation.
  n_before <- nrow(edge_df)
  edge_df <- edge_df[!is.na(edge_df$gene_a) & !is.na(edge_df$gene_b), ,
                     drop = FALSE]
  n_dropped <- n_before - nrow(edge_df)

  # Reorder to the canonical column order, keeping any extra user columns.
  canonical <- c("gene_a", "gene_b", "relationship",
                 "duplication_level", "seq_similarity")
  extra <- setdiff(colnames(edge_df), canonical)
  edge_df <- edge_df[, c(canonical, extra), drop = FALSE]
  rownames(edge_df) <- NULL

  class(edge_df) <- c("homology_graph", "data.frame")

  if (isTRUE(verbose)) {
    message(
      "Built homology_graph: ", nrow(edge_df), " edges (",
      sum(edge_df$relationship == "ortholog"), " ortholog, ",
      sum(edge_df$relationship == "paralog"), " paralog, ",
      sum(edge_df$relationship == "unclassified"), " unclassified)",
      if (n_dropped > 0) {
        paste0("; dropped ", n_dropped, " edge(s) with missing IDs")
      } else {
        ""
      },
      "."
    )
  }

  edge_df
}

#' @rdname as_homology_graph
#' @param object A \code{homology_graph}.
#' @param ... Ignored, for S3 compatibility.
#' @export
summary.homology_graph <- function(object, ...) {
  n_edges <- nrow(object)
  n_ortho <- sum(object$relationship == "ortholog")
  n_para <- sum(object$relationship == "paralog")
  n_uncl <- sum(object$relationship == "unclassified")

  # Composition percentages, computable from edges alone.
  pct_ortholog <- if (n_edges > 0) 100 * n_ortho / n_edges else NA_real_
  pct_paralog <- if (n_edges > 0) 100 * n_para / n_edges else NA_real_
  pct_unclassified <- if (n_edges > 0) 100 * n_uncl / n_edges else NA_real_

  # Many-to-many: genes appearing in more than one edge, per side.
  a_multi <- sum(table(object$gene_a) > 1)
  b_multi <- sum(table(object$gene_b) > 1)

  n_unique_a <- length(unique(object$gene_a))
  n_unique_b <- length(unique(object$gene_b))

  # Unmapped-ID rate is only knowable at build time; the source builder
  # records it as an attribute. Report it when present.
  unmapped_rate <- attr(object, "unmapped_rate")

  cat("homology_graph summary\n")
  cat("  edges:                ", n_edges, "\n", sep = "")
  cat("  ortholog edges:       ", n_ortho,
      sprintf(" (%.1f%%)", pct_ortholog), "\n", sep = "")
  cat("  paralog edges:        ", n_para,
      sprintf(" (%.1f%%)", pct_paralog), "\n", sep = "")
  cat("  unclassified edges:   ", n_uncl,
      sprintf(" (%.1f%%)", pct_unclassified), "\n", sep = "")
  cat("  unique genes (A):     ", n_unique_a, "\n", sep = "")
  cat("  unique genes (B):     ", n_unique_b, "\n", sep = "")
  cat("  multi-mapping A genes:", a_multi, "\n", sep = "")
  cat("  multi-mapping B genes:", b_multi, "\n", sep = "")
  if (!is.null(unmapped_rate)) {
    cat("  unmapped-ID rate:     ",
        sprintf("%.1f%%", 100 * unmapped_rate), "\n", sep = "")
  } else {
    cat("  unmapped-ID rate:      not available (built from edge table)\n")
  }

  invisible(list(
    n_edges = n_edges,
    n_ortholog = n_ortho,
    n_paralog = n_para,
    n_unclassified = n_uncl,
    pct_ortholog = pct_ortholog,
    pct_paralog = pct_paralog,
    pct_unclassified = pct_unclassified,
    n_unique_a = n_unique_a,
    n_unique_b = n_unique_b,
    multi_a = a_multi,
    multi_b = b_multi,
    unmapped_rate = unmapped_rate
  ))
}
