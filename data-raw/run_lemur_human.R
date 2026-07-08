#!/usr/bin/env Rscript
# =====================================================================
# validation/01_lemur_human.R
# ---------------------------------------------------------------------
# Regenerate the lemur–human paralogswap pipeline (Stages 0–6), persist
# every intermediate to S3, and assert the published RAMP checkpoint.
#
# Signatures/defaults below were verified against R/ in the repo (not the
# V1 API doc, which abstracts species-pairs into single objects and is
# stale on match_metacells / normalization tokens / id_type).
#
# Run from repo root inside the r-paralogswap conda env:
#     conda activate r-paralogswap      # confirm prompt shows (r-paralogswap)
#     Rscript validation/01_lemur_human.R
#
# Assumes (same as your existing atlas pulls): aws CLI configured for
# s3://paralogswap-data, and data-raw/lemur_lung.rds present locally.
# =====================================================================

set.seed(1)

# ---- package load (no devtools) -------------------------------------
stopifnot(dir.exists("R"))
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# =====================================================================
# CONFIG  (values confirmed against R/ unless marked CONFIRM)
# =====================================================================
# Species tokens are Ensembl biomaRt dataset PREFIXES (connect() does
# paste0(sp,"_gene_ensembl")) — NOT friendly names. lemur = mmurinus.
LEMUR <- "mmurinus"   # species_a / anchor
HUMAN <- "hsapiens"   # species_b

ENSEMBL_RELEASE  <- 116          # matches ramp_check.log ("Ensembl Genes 116")
MIRRORS          <- c("useast", "asia", NA)  # tried in order; NA = default host (www).
                                 # Main site flaps hard; mirrors are usually up but biomaRt
                                 # won't combine a mirror with the release pin, so mirror pulls
                                 # use the CURRENT release (== 116 now). Put NA first instead if
                                 # you want the explicit release pin and www is healthy.
RESUME           <- TRUE         # reuse clustered objects already in S3 instead of re-clustering
ID_TYPE          <- "symbol"     # atlases are SYMBOL-keyed (var/_index = RAMP1 etc.),
                                 # so the graph must be too, or match_clusters finds
                                 # no bridge. "ensembl" would need remapping both atlases.
                                 # Unnamed lemur loci (LOC…) drop out — expected.
MIN_CONFIDENCE   <- "low"        # "low" keeps all edges; "high" = confidence==1 only
MIN_PERC_ID      <- NULL

NORMALIZATION <- "lognorm"       # c("lognorm","sctransform") — NOT "LogNormalize"
RESOLUTION    <- 0.8             # revisit if RAMP triad doesn't reproduce (Section 4 hook)
N_PCS         <- 30
N_VAR_GENES   <- 2000
ASSAY         <- "RNA"
SLOT          <- "data"
COR_METHOD    <- "spearman"      # code default (c("spearman","pearson"))

GAMMA         <- 20
K_KNN         <- 5
MIN_CELLS     <- 10
MUTUAL_K      <- 5               # match_metacells is mutual-kNN; reciprocal-best -> ~9 pairs (wrong)
MIN_PAIRS     <- 10
MIN_CORRESPONDENCE <- 0.3
DELTA_THRESHOLD    <- 0.3

MAKE_PLOTS <- TRUE
FIG_DIR    <- "tests/manual/plots"  # already .gitignored; change if you add validation/figures/

# --- Stage 0 focal genes. In SYMBOL space (ID_TYPE="symbol") the graph, atlas
#     rownames, and focal ids all key on symbols, so RAMP1/2/3 are just the
#     symbols. Ensembl ids kept for reference (verified 3 ways: REST r116,
#     reverse lookup, ramp_check.log; all one2one, RAMP2 clean at 82.86%):
#       lemur  RAMP1 ENSMICG00000047662  RAMP2 ENSMICG00000028500  RAMP3 ENSMICG00000032638
#       human  RAMP1 ENSG00000132329     RAMP2 ENSG00000131477     RAMP3 ENSG00000122679
FOCAL_LEMUR <- c(RAMP1 = "RAMP1", RAMP2 = "RAMP2", RAMP3 = "RAMP3")
FOCAL_HUMAN <- c(RAMP1 = "RAMP1", RAMP2 = "RAMP2", RAMP3 = "RAMP3")

# lemur–human needs no gene-tree correction (verified). The RAMP2/VPS25 fix is
# mouse-only; keep this template for run_lemur_mouse.R. NOTE: in symbol space the
# correction edge is symbol-keyed too (human "RAMP2" -> mouse "Ramp2"), and the
# anchored pair is lemur–mouse, so gene_a/species_a are the lemur side there.
CORRECTIONS <- NULL
# e.g. for the lemur–mouse arm (id_type="symbol"), the human->mouse RAMP2 gap
# is repaired on the mouse side; build that correction against whichever pair
# actually loses the edge and tag source="eggnog".

# ---- S3 layout + helpers --------------------------------------------
S3_BUCKET   <- "s3://paralogswap-data"
S3_RESULTS  <- file.path(S3_BUCKET, "results")
S3_HOMOLOGY <- file.path(S3_BUCKET, "homology")

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
banner <- function(...) message("\n", strrep("=", 62), "\n", ..., "\n", strrep("=", 62))
soft_check <- function(label, value, expected, tol = 0.15) {
  ok <- is.numeric(value) && !is.na(value) &&
        abs(value - expected) <= tol * expected
  message(sprintf("  [%s] %s: got %s (expected ~%s)",
                  if (ok) "ok" else "DRIFT", label, value, expected)); invisible(ok)
}
# biomaRt has no retry, ignores options(timeout), and — critically — FORBIDS
# version + mirror together (it silently drops the mirror and hits the main
# host). Ensembl BioMart also flaps. So: on mirror attempts we drop the release
# pin. That's lossless here because the current Ensembl release IS 116 (verified
# via REST and ramp_check.log), so a mirror pull == release 116 anyway. The
# release pin is only honored on the default-host attempt.
build_graph_resilient <- function(..., mirrors = MIRRORS, tries_per = 2, wait = 8) {
  args <- list(...); pinned <- args$ensembl_release; last <- NULL; n <- 0L
  for (m in mirrors) for (t in seq_len(tries_per)) {
    n <- n + 1L; a <- args
    if (is.na(m)) {                        # default host: keep the exact release pin
      a$mirror <- NULL; a$ensembl_release <- pinned
    } else {                               # mirror: version+mirror is illegal, drop the pin
      a$mirror <- m;    a$ensembl_release <- NULL
    }
    out <- tryCatch(do.call(build_homology_graph, a), error = function(e) e)
    if (!inherits(out, "error")) {
      if (n > 1) message("  ... ok on attempt ", n,
                         if (is.na(m)) " (default host, release pinned)"
                         else paste0(" (mirror=", m, ", current release)"))
      return(out)
    }
    last <- out
    message(sprintf("  attempt %d failed%s: %s", n,
                    if (is.na(m)) " [default]" else paste0(" [mirror=", m, "]"),
                    conditionMessage(out)))
    Sys.sleep(wait)
  }
  stop(last)
}
# is an S3 object present? (for RESUME)
s3_exists <- function(key) {
  st <- suppressWarnings(system2("aws", c("s3", "ls", key), stdout = TRUE, stderr = TRUE))
  isTRUE(length(st) > 0) && !any(grepl("error|not *found", st, ignore.case = TRUE))
}

# Load the lemur atlas now — its gene set bounds the Stage 1a pull (below), and
# it is reused for clustering in Stage 2.
lemur_raw <- readRDS("data-raw/lemur_lung.rds")

# =====================================================================
# STAGE 1 — homology graphs (two of them, by design)
#   hg_broad : one2one orthologs over the atlas gene set, NO paralogs. Feeds
#              cluster + metacell matching. Restricted to rownames(lemur_raw);
#              this atlas is ~31.5k genes (near genome-wide), so the restriction
#              is modest — the real timeout defense is mirror-cycling retry.
#   hg_focal : RAMP orthologs + RAMP paralogs (both species). Feeds detection.
#              genes= restricts the pull, so include_paralogs=TRUE is safe here.
#   Both builds cycle mirrors via build_graph_resilient() (see MIRRORS in CONFIG).
# =====================================================================
lemur_genes <- rownames(lemur_raw)
message("lemur atlas genes available to bridge: ", length(lemur_genes))

banner("STAGE 1a  build_homology_graph  (atlas-restricted orthologs, ", LEMUR, " <-> ", HUMAN, ")")
hg_broad <- build_graph_resilient(
  species_a = LEMUR, species_b = HUMAN,
  genes = lemur_genes,                    # restrict to atlas (bounds pull; note atlas ~= genome here)
  ensembl_release = ENSEMBL_RELEASE,
  id_type = ID_TYPE, corrections = CORRECTIONS,
  include_paralogs = FALSE,               # broad graph: orthologs only
  min_perc_id = MIN_PERC_ID, min_confidence = MIN_CONFIDENCE, verbose = TRUE
)
s3_write_rds(hg_broad, file.path(S3_HOMOLOGY, "lemur_human_hg_broad.rds"))

banner("STAGE 1b  build_homology_graph  (focal RAMP orthologs + paralogs)")
hg_focal <- build_graph_resilient(
  species_a = LEMUR, species_b = HUMAN,
  genes = unname(FOCAL_LEMUR),            # restrict -> paralog pull is bounded
  ensembl_release = ENSEMBL_RELEASE,
  id_type = ID_TYPE, corrections = CORRECTIONS,
  include_paralogs = TRUE,                # focal graph: paralogs needed for r_paralog
  min_perc_id = MIN_PERC_ID, min_confidence = MIN_CONFIDENCE, verbose = TRUE
)
s3_write_rds(hg_focal, file.path(S3_HOMOLOGY, "lemur_human_hg_focal.rds"))
# sanity: RAMP paralog edges must exist on the HUMAN side (para_b in detection)
.para_human <- hg_focal[hg_focal$relationship == "paralog" & hg_focal$species_a == HUMAN, ]
message("  human-side RAMP paralog edges in hg_focal: ", nrow(.para_human),
        if (nrow(.para_human) == 0) "  <-- WARNING: detection will find no paralogs" else "")

# =====================================================================
# STAGE 2 — load atlases + cluster each species independently
#   list keys (lemur/human) are cosmetic; species tokens are passed separately.
#   Atlas rownames are symbols (RAMP1 etc.); the symbol-space graph matches them.
# =====================================================================
banner("STAGE 2  cluster_species")
lc_key <- file.path(S3_RESULTS, "lemur_lung_clustered.rds")
hc_key <- file.path(S3_RESULTS, "human_lung_clustered.rds")
if (isTRUE(RESUME) && s3_exists(lc_key) && s3_exists(hc_key)) {
  message("  RESUME: loading clustered objects from S3 (skipping re-cluster)")
  sc_lemur <- s3_read_rds(lc_key)
  sc_human <- s3_read_rds(hc_key)
} else {
  human_raw <- s3_read_rds(file.path(S3_RESULTS, "human_lung.rds"))  # lemur_raw loaded early
  sc <- cluster_species(
    object = list(lemur = lemur_raw, human = human_raw),
    species_col = "species", assay = ASSAY, normalization = NORMALIZATION,
    resolution = RESOLUTION, n_pcs = N_PCS, n_var_genes = N_VAR_GENES, verbose = TRUE
  )
  stopifnot(all(c("lemur", "human") %in% names(sc)))
  sc_lemur <- sc[["lemur"]]; sc_human <- sc[["human"]]
  s3_write_rds(sc_lemur, lc_key)
  s3_write_rds(sc_human, hc_key)
}

# =====================================================================
# STAGE 3 — cross-species cluster correspondence (uses hg_broad)
# =====================================================================
banner("STAGE 3  match_clusters")
correspondence <- match_clusters(
  clusters_a = sc_lemur, clusters_b = sc_human, homology_graph = hg_broad,
  species_a = LEMUR, species_b = HUMAN,
  cor_method = COR_METHOD, min_correspondence = MIN_CORRESPONDENCE,
  assay = ASSAY, slot = SLOT, verbose = TRUE
)
s3_write_rds(correspondence, file.path(S3_RESULTS, "lemur_human_cluster_correspondence.rds"))
message("  matched cluster pairs: ",
        tryCatch(nrow(correspondence), error = function(e) "inspect structure"))
cat("\n--- summary(cluster_correspondence) ---\n"); print(summary(correspondence))

# =====================================================================
# STAGE 4 — metacells, per species (one clustered object at a time)
# =====================================================================
banner("STAGE 4  build_metacells  (gamma=", GAMMA, ")")
mc_lemur <- build_metacells(sc_lemur, gamma = GAMMA, assay = ASSAY, slot = SLOT,
                            n_var_genes = N_VAR_GENES, k_knn = K_KNN,
                            min_cells = MIN_CELLS, cluster_col = NULL, verbose = TRUE)
mc_human <- build_metacells(sc_human, gamma = GAMMA, assay = ASSAY, slot = SLOT,
                            n_var_genes = N_VAR_GENES, k_knn = K_KNN,
                            min_cells = MIN_CELLS, cluster_col = NULL, verbose = TRUE)
s3_write_rds(mc_lemur, file.path(S3_RESULTS, "lemur_lung_metacells.rds"))
s3_write_rds(mc_human, file.path(S3_RESULTS, "human_lung_metacells.rds"))
soft_check("lemur metacells", ncol(mc_lemur$data), 1704)
soft_check("human metacells", ncol(mc_human$data), 2000)

# =====================================================================
# STAGE 5 — fine metacell correspondence (mutual-kNN; uses hg_broad)
#   NB arg order: correspondence BEFORE homology_graph
# =====================================================================
banner("STAGE 5  match_metacells  (mutual_k=", MUTUAL_K, ")")
metacell_correspondence <- match_metacells(
  metacells_a = mc_lemur, metacells_b = mc_human,
  correspondence = correspondence, homology_graph = hg_broad,
  species_a = LEMUR, species_b = HUMAN,
  cor_method = COR_METHOD, mutual_k = MUTUAL_K, verbose = TRUE
)
s3_write_rds(metacell_correspondence,
             file.path(S3_RESULTS, "lemur_human_metacell_correspondence.rds"))
soft_check("metacell pairs", nrow(metacell_correspondence), 192, tol = 0.20)

# =====================================================================
# STAGE 6 — homolog correlations (uses hg_focal) + substitution detection
# =====================================================================
banner("STAGE 6  compute_homolog_correlations + detect_substitutions")
homolog_correlations <- compute_homolog_correlations(
  metacells_a = mc_lemur, metacells_b = mc_human,
  metacell_correspondence = metacell_correspondence, homology_graph = hg_focal,
  species_a = LEMUR, species_b = HUMAN,
  focal_genes = unname(FOCAL_LEMUR),      # lemur (species-A) RAMP ids
  cor_method = COR_METHOD, min_pairs = MIN_PAIRS, verbose = TRUE
)
s3_write_rds(homolog_correlations,
             file.path(S3_RESULTS, "lemur_human_homolog_correlations.rds"))

substitutions <- detect_substitutions(
  homolog_correlations, delta_threshold = DELTA_THRESHOLD, require_ortholog = TRUE
)
s3_write_rds(substitutions, file.path(S3_RESULTS, "lemur_human_substitutions.rds"))

# =====================================================================
# REPRODUCTION GATE — published RAMP result (Ezran et al. 2025)
#   substitutions cols (real): gene, ortholog, best_paralog,
#                              r_ortholog, r_paralog, delta_r, flagged
#   `gene` holds the lemur (species-A) Ensembl id.
# =====================================================================
banner("REPRODUCTION GATE  (lemur RAMP1 flagged -> human RAMP2; RAMP2 conserved)")
subs <- as.data.frame(substitutions)
r1 <- subs[subs$gene == FOCAL_LEMUR[["RAMP1"]], , drop = FALSE]
r2 <- subs[subs$gene == FOCAL_LEMUR[["RAMP2"]], , drop = FALSE]

if (nrow(r1) == 1)
  message(sprintf("  RAMP1: r_ort=%.3f  r_par=%.3f  dr=%.3f  best_paralog=%s  flagged=%s",
                  r1$r_ortholog, r1$r_paralog, r1$delta_r, r1$best_paralog, r1$flagged))
if (nrow(r2) == 1)
  message(sprintf("  RAMP2: r_ort=%.3f  r_par=%.3f  dr=%.3f  flagged=%s",
                  r2$r_ortholog, r2$r_paralog, r2$delta_r, r2$flagged))

stopifnot(
  "RAMP1 present in substitutions"                 = nrow(r1) == 1,
  "RAMP1 must be flagged"                           = isTRUE(r1$flagged[1]),
  "RAMP1 delta_r reproduces ~0.60"                  = abs(r1$delta_r[1] - 0.60) < 0.12,
  "RAMP1 r_paralog exceeds r_ortholog by >0.30"     = (r1$r_paralog[1] - r1$r_ortholog[1]) > 0.30,
  "RAMP1 best paralog is human RAMP2"               = identical(as.character(r1$best_paralog[1]),
                                                                FOCAL_HUMAN[["RAMP2"]]),
  "RAMP2 must NOT be flagged (conserved)"           = (nrow(r2) == 0 || !isTRUE(r2$flagged[1]))
)
banner("GATE PASSED  \u2014 published RAMP substitution reproduced.")

# =====================================================================
# OPTIONAL — real plots. Write to FIG_DIR (git-ignored).
#   Reminder: the already-committed tests/manual/plots/*.png predate the
#   ignore rule and still need `git rm --cached` to untrack.
# =====================================================================
if (isTRUE(MAKE_PLOTS)) {
  banner("PLOTS -> ", FIG_DIR, "/")
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  save_png <- function(p, f, w = 7, h = 5) {
    ggplot2::ggsave(file.path(FIG_DIR, f), p, width = w, height = h, dpi = 150)
    message("  wrote ", file.path(FIG_DIR, f))
  }
  tryCatch(save_png(plot_substitution_elbow(substitutions),         "elbow.png"),
           error = function(e) message("  elbow skipped: ", conditionMessage(e)))
  tryCatch(save_png(plot_ortholog_vs_paralog(homolog_correlations), "ortholog_vs_paralog.png"),
           error = function(e) message("  ovp skipped: ", conditionMessage(e)))
  # plot_graining_curve(clusters, ...) is single-species; run on the clustered lemur obj
  tryCatch(save_png(plot_graining_curve(sc_lemur),                  "graining_curve.png"),
           error = function(e) message("  graining skipped: ", conditionMessage(e)))
}

banner("DONE  \u2014 all intermediates persisted under ", S3_RESULTS)
