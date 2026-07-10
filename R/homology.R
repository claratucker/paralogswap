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
#'   \item{relationship}{factor with levels \code{"ortholog"} and
#'     \code{"paralog"}, read directly from Ensembl's orthology type.}
#'   \item{ortholog_type}{factor with levels \code{"one2one"},
#'     \code{"one2many"}, \code{"many2many"}, and \code{NA} for within-species
#'     paralog edges. Taken from Ensembl Compara; drives \code{ambiguous}.}
#'   \item{duplication_level}{character or \code{NA}. Taxonomic level of the
#'     duplication node (from the gene tree / Compara); used to stratify
#'     paralogs by duplication age.}
#'   \item{perc_id}{numeric in \code{[0, 100]} or \code{NA}. Ensembl homolog
#'     percent identity.}
#'   \item{confidence}{integer/character or \code{NA}. Ensembl orthology
#'     confidence flag (1 = high, 0 = low).}
#'   \item{ambiguous}{logical. \code{TRUE} where \code{ortholog_type ==
#'     "many2many"}; derived, not supplied.}
#'   \item{source}{character. Provenance of the edge: one of
#'     \code{"ensembl"}, \code{"eggnog"}, \code{"genetree"}, \code{"manual"}.}
#' }
#'
#' Edges are directed from species A to species B: every row represents a focal
#' gene in A and a homologous partner in B. To assess substitution in both
#' lineages, run the pipeline in both directions (swap \code{species_a} and
#' \code{species_b}).
#'
#' @param edge_df A \code{data.frame} with at least \code{gene_a},
#'   \code{gene_b}, \code{relationship}, and \code{ortholog_type}. Optional
#'   columns (\code{duplication_level}, \code{perc_id}, \code{confidence},
#'   \code{source}) are added as \code{NA}/\code{"manual"} if absent.
#'   \code{ambiguous} is always (re)derived from \code{ortholog_type}.
#' @param verbose Print a brief construction message. Default \code{TRUE}.
#'
#' @return An object of class \code{homology_graph}.
#'
#' @examples
#' edges <- data.frame(
#'   gene_a = c("Actb", "Gapdh", "Myl2"),
#'   gene_b = c("ACTB", "GAPDH", "MYL7"),
#'   relationship = c("ortholog", "ortholog", "paralog"),
#'   ortholog_type = c("one2one", "one2one", NA)
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

  required <- c("gene_a", "gene_b", "relationship", "ortholog_type")
  missing_cols <- setdiff(required, colnames(edge_df))
  if (length(missing_cols) > 0) {
    stop(
      "`edge_df` is missing required column(s): ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }

  edge_df$gene_a <- as.character(edge_df$gene_a)
  edge_df$gene_b <- as.character(edge_df$gene_b)

  # relationship: two levels, read from the source's orthology typing.
  rel <- as.character(edge_df$relationship)
  valid_rel <- c("ortholog", "paralog")
  bad_rel <- setdiff(unique(rel), valid_rel)
  if (length(bad_rel) > 0) {
    stop(
      "`relationship` must be 'ortholog' or 'paralog'; found: ",
      paste(bad_rel, collapse = ", "), ".",
      call. = FALSE
    )
  }
  edge_df$relationship <- factor(rel, levels = valid_rel)

  # `action` is meaningful only in a corrections table consumed by
  # build_homology_graph(). Silently treating a removal as an addition would
  # inject the very edge the user asked to retract.
  if ("action" %in% colnames(edge_df)) {
    act <- tolower(trimws(as.character(edge_df$action)))
    if (any(act != "add")) {
      stop("`edge_df` contains removals (action = 'remove'). `as_homology_graph()` ",
           "only builds a graph from edges; pass removals to ",
           "`build_homology_graph(corrections = )`.", call. = FALSE)
    }
  }

  # ortholog_type: three levels for orthologs; NA allowed (e.g. paralog edges).
  ot <- as.character(edge_df$ortholog_type)
  valid_ot <- c("one2one", "one2many", "many2many")
  bad_ot <- setdiff(unique(ot[!is.na(ot)]), valid_ot)
  if (length(bad_ot) > 0) {
    stop(
      "`ortholog_type` must be one of 'one2one', 'one2many', 'many2many' ",
      "(or NA); found: ", paste(bad_ot, collapse = ", "), ".",
      call. = FALSE
    )
  }
  edge_df$ortholog_type <- factor(ot, levels = valid_ot)

  # Optional columns -> NA if absent.
  if (!"duplication_level" %in% colnames(edge_df)) {
    edge_df$duplication_level <- NA_character_
  } else {
    edge_df$duplication_level <- as.character(edge_df$duplication_level)
  }
  if (!"perc_id" %in% colnames(edge_df)) {
    edge_df$perc_id <- NA_real_
  } else {
    edge_df$perc_id <- as.numeric(edge_df$perc_id)
    if (any(edge_df$perc_id < 0 | edge_df$perc_id > 100, na.rm = TRUE)) {
      stop("`perc_id` must lie in [0, 100].", call. = FALSE)
    }
  }
  if (!"confidence" %in% colnames(edge_df)) {
    edge_df$confidence <- NA
  }
  # source: default to "manual" for user-supplied edges lacking provenance.
  if (!"source" %in% colnames(edge_df)) {
    edge_df$source <- "manual"
  } else {
    edge_df$source <- as.character(edge_df$source)
  }

  # ambiguous is ALWAYS derived from ortholog_type; never trusted from input.
  edge_df$ambiguous <- !is.na(edge_df$ortholog_type) &
    edge_df$ortholog_type == "many2many"

# Invariant: paralog edges are within-species; ortholog edges cross-species.
  para <- edge_df$relationship == "paralog"
  if (any(para & edge_df$species_a != edge_df$species_b)) {
    stop("paralog edges must have species_a == species_b (within-species).",
         call. = FALSE)
  }
  if (any(!para & edge_df$species_a == edge_df$species_b)) {
    stop("ortholog edges must have species_a != species_b (cross-species).",
         call. = FALSE)
  }

  # Drop rows with missing gene IDs; they can't participate in correlation.
  n_before <- nrow(edge_df)
  edge_df <- edge_df[!is.na(edge_df$gene_a) & !is.na(edge_df$gene_b), ,
                     drop = FALSE]
  n_dropped <- n_before - nrow(edge_df)

  canonical <- c("gene_a", "gene_b", "relationship", "ortholog_type",
                 "duplication_level", "perc_id", "confidence",
                 "ambiguous", "source")
  extra <- setdiff(colnames(edge_df), canonical)
  edge_df <- edge_df[, c(canonical, extra), drop = FALSE]
  rownames(edge_df) <- NULL

  class(edge_df) <- c("homology_graph", "data.frame")

  if (isTRUE(verbose)) {
    message(
      "Built homology_graph: ", nrow(edge_df), " edges (",
      sum(edge_df$relationship == "ortholog"), " ortholog, ",
      sum(edge_df$relationship == "paralog"), " paralog); ",
      sum(edge_df$ambiguous), " ambiguous (many2many)",
      if (n_dropped > 0) paste0("; dropped ", n_dropped, " with missing IDs") else "",
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
  n_para  <- sum(object$relationship == "paralog")
  ex <- attr(object, "excluded")
  if (!is.null(ex) && nrow(ex))
    cat("  excluded by correction:", nrow(ex), "edge(s)\n")

  pct <- function(n) if (n_edges > 0) 100 * n / n_edges else NA_real_

  # ortholog_type breakdown (orthologs only carry a type).
  ot <- object$ortholog_type
  n_o2o <- sum(ot == "one2one",  na.rm = TRUE)
  n_o2m <- sum(ot == "one2many", na.rm = TRUE)
  n_m2m <- sum(ot == "many2many", na.rm = TRUE)
  n_ambig <- sum(object$ambiguous)

  # confidence breakdown (Ensembl: 1 = high, 0 = low).
  conf <- suppressWarnings(as.numeric(as.character(object$confidence)))
  n_high <- sum(conf == 1, na.rm = TRUE)
  n_low  <- sum(conf == 0, na.rm = TRUE)

  # source breakdown.
  src_tab <- table(object$source)

  unmapped_rate <- attr(object, "unmapped_rate")

  cat("homology_graph summary\n")
  cat("  edges:              ", n_edges, "\n", sep = "")
  cat("  ortholog / paralog: ", n_ortho, " / ", n_para,
      sprintf(" (%.1f%% / %.1f%%)", pct(n_ortho), pct(n_para)), "\n", sep = "")
  cat("  ortholog_type:      ", n_o2o, " one2one, ", n_o2m, " one2many, ",
      n_m2m, " many2many\n", sep = "")
  cat("  ambiguous (m2m):    ", n_ambig,
      sprintf(" (%.1f%%)", pct(n_ambig)), "\n", sep = "")
  cat("  confidence:         ", n_high, " high, ", n_low, " low\n", sep = "")
  cat("  source:             ",
      paste(sprintf("%s=%d", names(src_tab), as.integer(src_tab)),
            collapse = ", "), "\n", sep = "")
  if (!is.null(unmapped_rate)) {
    cat("  unmapped-ID rate:   ", sprintf("%.1f%%", 100 * unmapped_rate),
        "\n", sep = "")
  } else {
    cat("  unmapped-ID rate:   not available (built from edge table)\n")
  }

  invisible(list(
    n_edges = n_edges, n_ortholog = n_ortho, n_paralog = n_para,
    n_one2one = n_o2o, n_one2many = n_o2m, n_many2many = n_m2m,
    n_ambiguous = n_ambig, n_high_conf = n_high, n_low_conf = n_low,
    source_counts = src_tab, unmapped_rate = unmapped_rate
  ))
}

#' Assemble homology edges from raw ortholog/paralog rows (internal)
#'
#' Pure logic: takes the two directional ortholog pulls and the per-species
#' paralog pulls (each already normalized to the edge-table columns), unions and
#' dedupes them, normalizes ortholog direction to the focal species, applies
#' filters, and merges corrections. No network. Separated from
#' \code{build_homology_graph} so the assembly logic is unit-testable.
#'
#' @param ortholog_edges data.frame of cross-species ortholog edges (both
#'   directions stacked), with the standard edge columns.
#' @param paralog_edges data.frame of within-species paralog edges (both
#'   species stacked), or NULL.
#' @param species_a,species_b Focal and target species prefixes.
#' @param corrections Optional user edge table.
#' @param mirror Optional Ensembl mirror: one of "useast", "uswest", "asia".
#'   Can be faster/more reliable than the default www.ensembl.org for large pulls.
#' @param min_perc_id,min_confidence Filters (see build_homology_graph).
#' @param verbose Passed to as_homology_graph.
#' @return A homology_graph.
#' @keywords internal
assemble_homology_edges <- function(ortholog_edges, paralog_edges,
                                    species_a, species_b,
                                    corrections = NULL,
                                    min_perc_id = NULL,
                                    min_confidence = "low",
                                    verbose = TRUE) {

  edges <- do.call(rbind, Filter(Negate(is.null),
                                 list(ortholog_edges, paralog_edges)))
  if (is.null(edges) || nrow(edges) == 0) {
    stop("No homology edges to assemble.", call. = FALSE)
  }

  # ---- union: dedupe pairs seen in both directions ---------------------------
  make_key <- function(df) {
    with(df, ifelse(
      species_a == species_b,
      paste(relationship, species_a, pmin(gene_a, gene_b), pmax(gene_a, gene_b)),
      paste(relationship,
            pmin(paste(species_a, gene_a), paste(species_b, gene_b)),
            pmax(paste(species_a, gene_a), paste(species_b, gene_b)))
    ))
  }
  key <- make_key(edges)
  ord <- order(-(!is.na(edges$confidence) & edges$confidence == 1),
               -ifelse(is.na(edges$perc_id), -1, edges$perc_id))
  edges <- edges[ord, , drop = FALSE]
  edges <- edges[!duplicated(key[ord]), , drop = FALSE]

  # ---- normalize ortholog direction to focal species (species_a) -------------
  flip <- edges$relationship == "ortholog" & edges$species_a == species_b
  if (any(flip)) {
    edges[flip, c("gene_a", "gene_b")]       <- edges[flip, c("gene_b", "gene_a")]
    edges[flip, c("species_a", "species_b")] <- edges[flip, c("species_b", "species_a")]
  }

  # ---- optional filters ------------------------------------------------------
  if (!is.null(min_perc_id)) {
    edges <- edges[is.na(edges$perc_id) | edges$perc_id >= min_perc_id, ,
                   drop = FALSE]
  }
  if (min_confidence == "high") {
    conf <- suppressWarnings(as.numeric(as.character(edges$confidence)))
    edges <- edges[is.na(conf) | conf == 1, , drop = FALSE]
  }

# ---- apply corrections: removals first, then additions ---------------------
  excluded <- NULL
  if (!is.null(corrections)) {
    corrections <- as.data.frame(corrections, stringsAsFactors = FALSE)
    if (!"action" %in% names(corrections)) corrections$action <- "add"
    corrections$action <- tolower(trimws(corrections$action))
    bad <- setdiff(unique(corrections$action), c("add", "remove"))
    if (length(bad))
      stop("`corrections$action` must be \"add\" or \"remove\"; got: ",
           paste(bad, collapse = ", "), call. = FALSE)

    rm_df  <- corrections[corrections$action == "remove", , drop = FALSE]
    add_df <- corrections[corrections$action == "add",    , drop = FALSE]

    if (isTRUE(verbose))
      message("Corrections: removing ", nrow(rm_df), ", adding ", nrow(add_df), " edge(s) ...")

if (nrow(rm_df)) {
      # Validate schema and compute keys, but strip `action` first: the guard in
      # as_homology_graph() exists to stop users passing removals to it, and must
      # not fire on this internal key-computation call.
      rm_key <- make_key(as_homology_graph(rm_df[, setdiff(names(rm_df), "action"),
                                                 drop = FALSE], verbose = FALSE))
      ek     <- make_key(edges)
      hit    <- ek %in% rm_key
      excluded <- edges[hit, , drop = FALSE]
      unmatched <- setdiff(rm_key, ek)
if (length(unmatched) && isTRUE(verbose))
        message("  ", length(unmatched), " removal edge(s) matched nothing ",
                "(expected when a bounded pull never fetched them).")

      edges <- edges[!hit, , drop = FALSE]
    }

    if (nrow(add_df)) {
      corr     <- as_homology_graph(add_df, verbose = FALSE)
      edges    <- edges[!(make_key(edges) %in% make_key(corr)), , drop = FALSE]
      common   <- intersect(colnames(edges), colnames(corr))
      edges    <- rbind(edges[, common, drop = FALSE], corr[, common, drop = FALSE])
    }
  }

  out <- as_homology_graph(edges, verbose = verbose)
  attr(out, "excluded") <- excluded    # NULL when nothing removed
  out
}

#' Build a homology graph from Ensembl Compara
#'
#' Fetches cross-species ortholog edges and within-species paralog edges for a
#' pair of species from Ensembl Compara via biomaRt (Durinck et al. 2009), and
#' assembles them into a symmetric \code{homology_graph}. Cross-species
#' orthologs are pulled in \emph{both} directions and unioned, because Ensembl's
#' directional homolog tables can disagree (e.g. in release 116, human RAMP2
#' returns no mouse ortholog while mouse Ramp2 is mis-assigned to the VPS25
#' orthogroup). Within-species paralogs are pulled for \emph{both} species, so
#' the graph supports detection in either direction without rebuilding.
#'
#' The package fetches curated orthology but does not attempt to detect or
#' repair its errors automatically. Where curated calls are known to be
#' incomplete or wrong, supply a \code{corrections} edge table (e.g. derived
#' from eggNOG orthogroups or a reconciled gene tree); it is merged with the
#' Ensembl edges and tagged in the \code{source} column. See the validation
#' vignette for the RAMP2 correction.
#'
#' @param species_a,species_b Ensembl species prefixes, e.g. \code{"hsapiens"},
#'   \code{"mmusculus"}, \code{"mmurinus"}.
#' @param genes Optional character vector of species-A gene identifiers to
#'   restrict the query to (recommended for targeted analyses; a whole-genome
#'   pull is large). \code{NULL} pulls all genes with homologs.
#' @param ensembl_release Optional integer Ensembl release to pin for
#'   reproducibility. \code{NULL} uses the current release.
#' @param id_type Either \code{"ensembl"} (gene stable IDs) or \code{"symbol"}
#'   (external gene names); controls both the filter and the returned IDs.
#' @param corrections Optional edge table of user-supplied homology, merged with
#'   the pulled Ensembl edges and tagged in the \code{source} column (e.g.
#'   \code{"eggnog"}, \code{"genetree"}). Requires the columns \code{gene_a},
#'   \code{gene_b}, \code{species_a}, \code{species_b}, \code{relationship} and
#'   \code{ortholog_type}, and must be written in the same identifier space as
#'   \code{id_type}: edges join by string match and no conversion is performed,
#'   so an Ensembl-ID correction will not attach to a symbol-space graph.
#'
#'   An optional \code{action} column takes \code{"add"} (the default when the
#'   column is absent) or \code{"remove"}. Removals are applied before
#'   additions, so a pair may be retracted and re-supplied in one table; they
#'   match irrespective of edge orientation, and the retracted rows are kept on
#'   the result in \code{attr(x, "excluded")} with their original
#'   \code{source}, so an override is recorded rather than hidden. A removal
#'   that matches no pulled edge is reported and otherwise ignored, which is
#'   expected whenever a bounded pull never fetched the offending edge.
#'
#'   Removal exists because curation cannot be inferred from expression. Unlike
#'   SAMap, which prunes its homology graph to positively correlated gene pairs,
#'   paralogswap must not select edges using the correlations detection
#'   measures -- that would tune the apparatus on the result. A wrong edge is
#'   therefore retracted by declaration, with provenance, not by data. This
#'   matters because \code{detect_substitutions} takes \code{max(r_paralog)}:
#'   supplying correct edges cannot out-vote an incorrect one, since a spurious
#'   paralog may win the maximum by chance.
#' @param min_perc_id Optional numeric floor on Ensembl homolog percent
#'   identity. \code{NULL} keeps all edges and relies on \code{min_confidence}.
#' @param min_confidence Character; keep orthologs at or above this Ensembl
#'   confidence. \code{"low"} keeps all; \code{"high"} keeps only confidence-1
#'   calls. Default \code{"low"} (keep everything, let the user filter).
#' @param verbose Print progress and the construction message. Default TRUE.
#' @param mirror Optional biomaRt mirror (e.g. "useast"); NULL uses the default host.
#' @param include_paralogs Logical; if FALSE, skip the paralog pull (orthologs only).
#' @return A \code{homology_graph}.
#'
#' @examples
#' \dontrun{
#' ramp <- c("ENSG00000132329", "ENSG00000131477", "ENSG00000122679")
#' hg <- build_homology_graph("hsapiens", "mmusculus", genes = ramp)
#' summary(hg)
#' }
#'
#' @export
build_homology_graph <- function(species_a, species_b,
                                 genes = NULL,
                                 ensembl_release = NULL,
                                 mirror = NULL,
                                 id_type = c("ensembl", "symbol"),
                                 corrections = NULL,
                                 include_paralogs = TRUE,
                                 min_perc_id = NULL,
                                 min_confidence = c("low", "high"),
                                 verbose = TRUE) {
  id_type <- match.arg(id_type)
  min_confidence <- match.arg(min_confidence)

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("build_homology_graph requires the 'biomaRt' package.", call. = FALSE)
  }

  filt <- if (id_type == "ensembl") "ensembl_gene_id" else "external_gene_name"
  id_attr <- filt
  name_attr <- "external_gene_name"

connect <- function(sp) {
    args <- list(biomart = "genes", dataset = paste0(sp, "_gene_ensembl"))
    if (!is.null(ensembl_release)) args$version <- ensembl_release
    if (!is.null(mirror))          args$mirror  <- mirror
    do.call(biomaRt::useEnsembl, args)
  }

  if (isTRUE(verbose)) message("Connecting to Ensembl (", species_a,
                               ", ", species_b, ") ...")
  mart_a <- connect(species_a)
  mart_b <- connect(species_b)

  # ---- helper: pull cross-species ortholog edges in ONE direction ------------
  # Queries the FROM mart for homologs in the TO species. Returns edges directed
  # from -> to, with ortholog_type / perc_id / confidence.
  pull_orthologs <- function(from_mart, from_sp, to_sp, restrict_ids) {
    pre <- to_sp   # biomaRt homolog attribute prefix is the target species
    attrs <- c(
      id_attr, name_attr,
      paste0(pre, "_homolog_ensembl_gene"),
      paste0(pre, "_homolog_associated_gene_name"),
      paste0(pre, "_homolog_orthology_type"),
      paste0(pre, "_homolog_perc_id"),
      paste0(pre, "_homolog_orthology_confidence")
    )
    args <- list(attributes = attrs, mart = from_mart)
    if (!is.null(restrict_ids)) {
      args$filters <- filt
      args$values  <- restrict_ids
    }
    res <- do.call(biomaRt::getBM, args)
    if (nrow(res) == 0) return(NULL)
    names(res) <- c("from_id", "from_name", "to_id", "to_name",
                    "otype_raw", "perc_id", "confidence")
    # Drop rows with no homolog (blank target id).
    res <- res[!is.na(res$to_id) & res$to_id != "", , drop = FALSE]
    if (nrow(res) == 0) return(NULL)
    # Map Ensembl orthology_type -> our ortholog_type levels.
    ot <- rep(NA_character_, nrow(res))
    ot[res$otype_raw == "ortholog_one2one"]  <- "one2one"
    ot[res$otype_raw == "ortholog_one2many"] <- "one2many"
    ot[res$otype_raw == "ortholog_many2many"] <- "many2many"
    keep <- !is.na(ot)   # keep only ortholog rows (drops any stray paralog rows)
    data.frame(
      gene_a = if (id_type == "ensembl") res$from_id[keep] else res$from_name[keep],
      gene_b = if (id_type == "ensembl") res$to_id[keep]   else res$to_name[keep],
      species_a = from_sp, species_b = to_sp,
      relationship = "ortholog",
      ortholog_type = ot[keep],
      duplication_level = NA_character_,
      perc_id = res$perc_id[keep],
      confidence = res$confidence[keep],
      source = "ensembl",
      stringsAsFactors = FALSE
    )
  }

  # ---- helper: pull WITHIN-species paralog edges -----------------------------
  pull_paralogs <- function(mart, sp, restrict_ids) {
    pre <- sp
    attrs <- c(
      id_attr, name_attr,
      paste0(pre, "_paralog_ensembl_gene"),
      paste0(pre, "_paralog_associated_gene_name"),
      paste0(pre, "_paralog_perc_id"),
      paste0(pre, "_paralog_subtype")   # taxonomic node of the duplication event
    )
    args <- list(attributes = attrs, mart = mart)
    if (!is.null(restrict_ids)) {
      args$filters <- filt
      args$values  <- restrict_ids
    }
    res <- do.call(biomaRt::getBM, args)
    if (nrow(res) == 0) return(NULL)
    names(res) <- c("from_id", "from_name", "to_id", "to_name", "perc_id",
                    "subtype")
    res <- res[!is.na(res$to_id) & res$to_id != "", , drop = FALSE]
    if (nrow(res) == 0) return(NULL)
    res$subtype[is.na(res$subtype) | !nzchar(trimws(res$subtype))] <- NA_character_
    data.frame(
      gene_a = if (id_type == "ensembl") res$from_id else res$from_name,
      gene_b = if (id_type == "ensembl") res$to_id   else res$to_name,
      species_a = sp, species_b = sp,
      relationship = "paralog",
      ortholog_type = NA_character_,   # paralogs carry no cross-species type
      duplication_level = res$subtype,  # e.g. "Vertebrata"; the duplication node
      perc_id = res$perc_id,
      confidence = NA,
      source = "ensembl",
      stringsAsFactors = FALSE
    )
  }

  # ---- fetch: orthologs BOTH directions, paralogs BOTH species ---------------
  if (isTRUE(verbose)) message("Pulling orthologs (both directions) ...")
  o_ab <- pull_orthologs(mart_a, species_a, species_b, genes)
  # Reverse direction restricted to the B genes we actually found (keeps it small).
# Genes with no ortholog return an empty gene_b, and biomaRt treats an empty
  # filter value as "no filter" -- which would make the reverse pull unrestricted.
  b_ids <- if (!is.null(o_ab)) .clean_ids(o_ab$gene_b) else NULL

  o_ba <- pull_orthologs(mart_b, species_b, species_a,
                         if (is.null(genes)) NULL else b_ids)


if (isTRUE(include_paralogs)) {
    if (isTRUE(verbose)) message("Pulling paralogs (both species) ...")
    a_ids <- genes
    p_a <- pull_paralogs(mart_a, species_a, a_ids)
    p_b <- pull_paralogs(mart_b, species_b, b_ids)
  } else {
    p_a <- NULL; p_b <- NULL
  }

  ortholog_edges <- do.call(rbind, Filter(Negate(is.null), list(o_ab, o_ba)))
  paralog_edges  <- do.call(rbind, Filter(Negate(is.null), list(p_a, p_b)))

  assemble_homology_edges(
    ortholog_edges = ortholog_edges,
    paralog_edges  = paralog_edges,
    species_a      = species_a,
    species_b      = species_b,
    corrections    = corrections,
    min_perc_id    = min_perc_id,
    min_confidence = min_confidence,
    verbose        = verbose
  )
}


# biomaRt returns "" (not NA, not absent) for a gene with no homolog, and treats
# an empty filter value as "no filter" -- so passing it back unsanitized makes the
# reverse pull unrestricted, sweeping in unrelated genes. Drop empties and NAs
# before using IDs as filters.
.clean_ids <- function(x) {
  if (is.null(x)) return(NULL)
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0) NULL else x
}
