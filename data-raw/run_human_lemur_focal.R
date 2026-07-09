#!/usr/bin/env Rscript
# =====================================================================
# data-raw/run_human_lemur_focal.R
# ---------------------------------------------------------------------
# Human-anchored focal RAMP run. Recovers the SECOND published triad,
# which the lemur-anchored arm structurally cannot flag.
#
# Why this exists
#   compute_homolog_correlations() scores focal genes of species A against
#   homologs in species B, so a lemur-anchored run can only flag lemur genes.
#   Triad 1 (human RAMP2 - lemur RAMP1* - mouse Ramp2) has its outlier on
#   lemur: recovered, delta_r = 0.602. Triad 2 (human RAMP3* - lemur RAMP2 -
#   mouse Ramp2) has its outlier on HUMAN. Anchoring on human addresses it.
#
# What this run establishes (and what it cannot)
#   Human RAMP3 is expressed in lung (101/192 matched metacells, sd = 0.669).
#   Its annotated ortholog, lemur RAMP3, is silent: nonzero in 8 of 1,704
#   lemur metacells atlas-wide, sd = 0.0019, and invariant across all 192
#   matched pairs. So r_ortholog is UNDEFINED, not low, and
#   delta_r = r_paralog - r_ortholog cannot be computed. detect_substitutions
#   correctly declines to score the gene (require_ortholog = TRUE).
#
#   The finding is nonetheless present in the correlation table: human RAMP3
#   tracks lemur RAMP2 (r = 0.486) above lemur RAMP1 (r = 0.407) and not its
#   own ortholog at all -- consistent with Ezran et al. 2025. This script
#   therefore asserts Triad 2 as an UNSCOREABLE substitution: an ortholog that
#   has been lost or silenced, with the paralog carrying the pattern.
#
#   This is a real limitation of delta_r, not a defect of this run: the method
#   is least sensitive at the completest substitutions, where the ortholog is
#   gone rather than merely diverged. Report it; do not paper over it by
#   setting require_ortholog = FALSE, which would silently substitute an
#   assumed-zero baseline for a measured one.
#
# Cost: no clustering, no metacell rebuild. Reuses the metacells and metacell
#   grid already in S3 from the lemur-human arm.
#
# Validity of reusing the lemur-anchored grid
#   match_metacells() uses MUTUAL k-nearest-neighbor matching: a pair (a,b) is
#   kept iff each is in the other's top-k. That condition is symmetric, so the
#   pair set is identical from either side and may be relabelled by swapping
#   columns. (Reciprocal-best-per-A would NOT be safely flippable.) The script
#   asserts this: ortholog correlations must reproduce the lemur-anchored run.
#
#     conda activate r-paralogswap
#     Rscript data-raw/run_human_lemur_focal.R
# =====================================================================

set.seed(1)
stopifnot(dir.exists("R"))
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# ---- CONFIG ---------------------------------------------------------
HUMAN <- "hsapiens"    # species_a / anchor for THIS run
LEMUR <- "mmurinus"    # species_b

ENSEMBL_RELEASE <- 116
MIRRORS         <- c("useast", "asia", NA)   # NA = default host (www)
ID_TYPE         <- "symbol"
MIN_CONFIDENCE  <- "low"
MIN_PERC_ID     <- NULL
COR_METHOD      <- "spearman"
MIN_PAIRS       <- 10
DELTA_THRESHOLD <- 0.3
RESUME          <- TRUE

FOCAL <- c("RAMP1", "RAMP2", "RAMP3")   # human and lemur both uppercase
CORRECTIONS <- NULL                     # Compara is correct for lemur-human

# Ortholog r values from the lemur-anchored run. Correlation is symmetric, so
# a correct flip must reproduce them exactly. Guards silent grid corruption.
EXPECT_R_ORTHOLOG <- c(RAMP1 = -0.05839964, RAMP2 = 0.62541372)

S3_BUCKET   <- "s3://paralogswap-data"
S3_RESULTS  <- file.path(S3_BUCKET, "results")
S3_HOMOLOGY <- file.path(S3_BUCKET, "homology")

# ---- helpers --------------------------------------------------------
s3_read_rds <- function(key) {
  tmp <- tempfile(fileext = ".rds"); on.exit(unlink(tmp), add = TRUE)
  st <- system2("aws", c("s3", "cp", key, tmp), stdout = TRUE, stderr = TRUE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0)
    stop("S3 download failed for ", key, "\n", paste(st, collapse = "\n"))
  message("  \u2193 read  ", key); readRDS(tmp)
}
s3_write_rds <- function(obj, key) {
  tmp <- tempfile(fileext = ".rds"); on.exit(unlink(tmp), add = TRUE)
  saveRDS(obj, tmp)
  st <- system2("aws", c("s3", "cp", tmp, key), stdout = TRUE, stderr = TRUE)
  status <- attr(st, "status")
  if (!is.null(status) && status != 0)
    stop("S3 upload failed for ", key, "\n", paste(st, collapse = "\n"))
  message("  \u2191 saved ", key); invisible(key)
}
s3_exists <- function(key) {
  st <- suppressWarnings(system2("aws", c("s3", "ls", key), stdout = TRUE, stderr = TRUE))
  isTRUE(length(st) > 0) && !any(grepl("error|not *found", st, ignore.case = TRUE))
}
resume_or <- function(key, build_fn) {
  if (isTRUE(RESUME) && s3_exists(key)) { message("  RESUME: loading ", key); return(s3_read_rds(key)) }
  obj <- build_fn(); s3_write_rds(obj, key); obj
}
banner <- function(...) message("\n", strrep("=", 62), "\n", ..., "\n", strrep("=", 62))

# biomaRt forbids version + mirror together (it drops the mirror silently), so
# mirror attempts drop the release pin; the current release IS 116 anyway.
build_graph_resilient <- function(..., mirrors = MIRRORS, tries_per = 2, wait = 8) {
  args <- list(...); pinned <- args$ensembl_release; last <- NULL; n <- 0L
  for (m in mirrors) for (t in seq_len(tries_per)) {
    n <- n + 1L; a <- args
    if (is.na(m)) { a$mirror <- NULL; a$ensembl_release <- pinned }
    else          { a$mirror <- m;    a$ensembl_release <- NULL }
    out <- tryCatch(do.call(build_homology_graph, a), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    last <- out
    message(sprintf("  attempt %d failed%s: %s", n,
                    if (is.na(m)) " [default]" else paste0(" [mirror=", m, "]"),
                    conditionMessage(out))); Sys.sleep(wait)
  }
  stop(last)
}

# Relabel a lemur-anchored metacell_correspondence as human-anchored.
# Valid only because mutual-kNN pairing is symmetric (see header).
flip_correspondence <- function(mmc) {
  out <- data.frame(
    metacell_a = mmc$metacell_b,   # human metacells become side A
    metacell_b = mmc$metacell_a,   # lemur metacells become side B
    cluster_a  = mmc$cluster_b,
    cluster_b  = mmc$cluster_a,
    score      = mmc$score,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$cluster_a, -out$score), ]; rownames(out) <- NULL
  bn <- attr(mmc, "n_bridge_genes"); if (is.null(bn)) bn <- attr(mmc, "n_orthologs")
  attr(out, "n_bridge_genes") <- bn   # A-side unmatched count no longer applies
  class(out) <- c("metacell_correspondence", "data.frame")
  out
}

# Why a focal gene could not be scored. Reads attr(subs,"dropped") when the
# package provides it; otherwise derives the same verdict from the correlation
# table, so this script is correct before and after R/detection.R is patched.
diagnose <- function(gene, hc, subs, n_grid, min_pairs) {
  d <- attr(subs, "dropped")
  if (!is.null(d) && gene %in% d$gene) {
    r <- d[d$gene == gene, ][1, ]
    return(list(reason = r$reason, ortholog = r$ortholog, n_pairs = r$n_pairs,
                best_paralog = r$best_paralog, r_paralog = r$r_paralog,
                source = "attr(subs, \"dropped\")"))
  }
  sub  <- hc[hc$gene_a == gene, ]
  orth <- sub[sub$relationship == "ortholog", ]
  para <- sub[sub$relationship == "paralog" & is.finite(sub$r), ]
  reason <- if (nrow(orth) == 0) "no_ortholog_edge"
    else if (is.finite(orth$r[1])) "scored"
    else if (orth$n_pairs[1] < min_pairs) "ortholog_too_sparse"
    else "ortholog_invariant"   # full coverage, no variance -> ortholog silent
  list(reason = reason,
       ortholog     = if (nrow(orth)) orth$gene_b[1]  else NA_character_,
       n_pairs      = if (nrow(orth)) orth$n_pairs[1] else NA_integer_,
       best_paralog = if (nrow(para)) para$gene_b[which.max(para$r)] else NA_character_,
       r_paralog    = if (nrow(para)) max(para$r)     else NA_real_,
       source = "derived from homolog_correlations")
}

# =====================================================================
# STAGE 0 — focal homology graph, HUMAN-anchored
#   Built human-first, so gene_a is human and compute_homolog_correlations
#   needs no reorientation. para_b then resolves to the LEMUR paralogs, which
#   is what human-anchored detection must compare against. genes= bounds the
#   paralog pull (unbounded, it exceeds biomaRt's request ceiling).
# =====================================================================
banner("STAGE 0  focal homology graph  (", HUMAN, " -> ", LEMUR, ")")
hg_focal <- resume_or(file.path(S3_HOMOLOGY, "human_lemur_hg_focal.rds"), function()
  build_graph_resilient(
    species_a = HUMAN, species_b = LEMUR,
    genes = FOCAL, ensembl_release = ENSEMBL_RELEASE,
    id_type = ID_TYPE, corrections = CORRECTIONS,
    include_paralogs = TRUE,
    min_perc_id = MIN_PERC_ID, min_confidence = MIN_CONFIDENCE, verbose = TRUE))

.para_lemur <- hg_focal[hg_focal$relationship == "paralog" & hg_focal$species_a == LEMUR, ]
message("  lemur-side RAMP paralog edges (para_b): ", nrow(.para_lemur))
stopifnot("focal graph must carry lemur-side paralog edges" = nrow(.para_lemur) > 0)

# =====================================================================
# STAGE 1 — reuse metacells + grid from the lemur-human arm
# =====================================================================
banner("STAGE 1  reuse metacells and metacell grid from S3")
mc_lemur <- s3_read_rds(file.path(S3_RESULTS, "lemur_lung_metacells.rds"))
mc_human <- s3_read_rds(file.path(S3_RESULTS, "human_lung_metacells.rds"))
mmc_lemur_anchored <- s3_read_rds(file.path(S3_RESULTS, "lemur_human_metacell_correspondence.rds"))

mmc <- flip_correspondence(mmc_lemur_anchored)
stopifnot(
  "flip must preserve the pair count"          = nrow(mmc) == nrow(mmc_lemur_anchored),
  "flipped metacell_a must index human grid"   = all(mmc$metacell_a %in% colnames(mc_human$data)),
  "flipped metacell_b must index lemur grid"   = all(mmc$metacell_b %in% colnames(mc_lemur$data))
)
message("  ", nrow(mmc), " matched metacell pairs, relabelled human-anchored")
s3_write_rds(mmc, file.path(S3_RESULTS, "human_lemur_metacell_correspondence.rds"))

# expression sanity: which RAMPs are measurable on this grid, per species?
banner("STAGE 1b  expression on the matched grid")
for (g in FOCAL) {
  xh <- mc_human$data[g, mmc$metacell_a]; xl <- mc_lemur$data[g, mmc$metacell_b]
  message(sprintf("  %-5s human sd=%.4f nonzero=%3d/%d | lemur sd=%.4f nonzero=%3d/%d",
                  g, stats::sd(xh), sum(xh > 0), length(xh),
                     stats::sd(xl), sum(xl > 0), length(xl)))
}

# =====================================================================
# STAGE 2 — correlations + detection, human-anchored
# =====================================================================
banner("STAGE 2  compute_homolog_correlations + detect_substitutions")
hc <- compute_homolog_correlations(
  metacells_a = mc_human, metacells_b = mc_lemur,   # A = human
  metacell_correspondence = mmc, homology_graph = hg_focal,
  species_a = HUMAN, species_b = LEMUR,
  focal_genes = FOCAL, cor_method = COR_METHOD,
  min_pairs = MIN_PAIRS, verbose = TRUE)
s3_write_rds(hc, file.path(S3_RESULTS, "human_lemur_homolog_correlations.rds"))

subs <- detect_substitutions(hc, delta_threshold = DELTA_THRESHOLD,
                             require_ortholog = TRUE)
s3_write_rds(subs, file.path(S3_RESULTS, "human_lemur_substitutions.rds"))

cat("\n--- homolog_correlations (all edges) ---\n"); print(as.data.frame(hc))
cat("\n--- substitutions (scoreable genes) ---\n");   print(as.data.frame(subs))

# =====================================================================
# GATE 1 — the flip is sound
#   Correlation is symmetric: ortholog r must reproduce the lemur-anchored run.
# =====================================================================
banner("GATE 1  flip fidelity (ortholog r must match the lemur-anchored run)")
hcd <- as.data.frame(hc)
for (g in names(EXPECT_R_ORTHOLOG)) {
  got <- hcd$r[hcd$gene_a == g & hcd$relationship == "ortholog"][1]
  message(sprintf("  %-5s r_ortholog = %+.6f (expected %+.6f)", g, got, EXPECT_R_ORTHOLOG[[g]]))
  stopifnot("ortholog r must match the lemur-anchored run" =
              isTRUE(abs(got - EXPECT_R_ORTHOLOG[[g]]) < 1e-6))
}

# =====================================================================
# GATE 2 — Triad 1 mirror, reported not asserted
#   Human RAMP2's ortholog (lemur RAMP2, r=0.625) still beats its best paralog
#   (lemur RAMP1, r=0.544), so delta_r < 0 and it is NOT flagged. Correct: the
#   substitution is a property of lemur RAMP1 having swapped, not of human
#   RAMP2. The relation is directional, and this demonstrates it.
# =====================================================================
banner("GATE 2  Triad 1 mirror (directionality check)")
d <- as.data.frame(subs)
r2 <- d[d$gene == "RAMP2", , drop = FALSE]
if (nrow(r2) == 1) {
  message(sprintf("  human RAMP2: best_paralog=%s r_ort=%.3f r_par=%.3f dr=%.3f flagged=%s",
                  r2$best_paralog, r2$r_ortholog, r2$r_paralog, r2$delta_r, r2$flagged))
  stopifnot("human RAMP2's best paralog must be lemur RAMP1" =
              identical(as.character(r2$best_paralog[1]), "RAMP1"))
  if (!isTRUE(r2$flagged[1]))
    message("  not flagged, as expected: substitution is directional (lemur RAMP1 -> human RAMP2)")
}

# =====================================================================
# GATE 3 — Triad 2: human RAMP3, an UNSCOREABLE substitution
#   Assert what the data support: the ortholog is invariant across the full
#   grid (silenced in lemur), delta_r is undefined, and the paralog carrying
#   the pattern is lemur RAMP2. This is Ezran et al.'s claim, recovered as far
#   as delta_r permits.
# =====================================================================
banner("GATE 3  Triad 2 (human RAMP3 -> lemur RAMP2, ortholog silenced)")
dx <- diagnose("RAMP3", hcd, subs, n_grid = nrow(mmc), min_pairs = MIN_PAIRS)
message("  verdict source: ", dx$source)
message(sprintf("  reason=%s  ortholog=%s  n_pairs=%s", dx$reason, dx$ortholog, dx$n_pairs))
message(sprintf("  best paralog=%s at r=%.3f", dx$best_paralog, dx$r_paralog))

r3_par <- hcd[hcd$gene_a == "RAMP3" & hcd$relationship == "paralog" & is.finite(hcd$r), ]
r3_par <- r3_par[order(-r3_par$r), ]

stopifnot(
  "human RAMP3 must NOT appear among scoreable substitutions" = !("RAMP3" %in% d$gene),
  "its ortholog must be lemur RAMP3"                           = identical(dx$ortholog, "RAMP3"),
  "the ortholog must be invariant, not sparse"                 = identical(dx$reason, "ortholog_invariant"),
  "ortholog observed across the full grid"                     = dx$n_pairs == nrow(mmc),
  "best paralog must be lemur RAMP2"                           = identical(dx$best_paralog, "RAMP2"),
  "RAMP2 must beat RAMP1 as paralog partner"                   = identical(r3_par$gene_b[1], "RAMP2"),
  "paralog correlation must be substantial"                    = r3_par$r[1] > 0.3
)
message(sprintf("\n  human RAMP3 tracks lemur RAMP2 (r=%.3f) above lemur RAMP1 (r=%.3f),",
                r3_par$r[1], r3_par$r[2]))
message("  and its own ortholog not at all (invariant). Consistent with Ezran et al. 2025.")
message("  delta_r is undefined here: the method cannot score a substitution whose")
message("  ortholog has been wholly silenced. Reported, not flagged.")

banner("GATE PASSED  \u2014 Triad 2 recovered as an unscoreable substitution.")
banner("DONE  \u2014 human-anchored focal results persisted under ", S3_RESULTS)
