#!/usr/bin/env Rscript
# =====================================================================
# data-raw/run_lemur_mouse.R
# ---------------------------------------------------------------------
# Lemur–mouse arm of the paralogswap validation. Regenerates Stages 0–6,
# persists every intermediate to S3, and asserts the RAMP checkpoint.
# Twin of run_lemur_human.R; differences are all mouse-specific:
#   * species_b = "mmusculus"; mouse symbols are title-case ("Ramp2")
#   * REQUIRES the RAMP2/VPS25 correction: Ensembl Compara mis-assigns mouse
#     Ramp2 to the VPS25 orthogroup (Ramp2/Vps25 adjacent on mouse chr11), so
#     lemur RAMP2 -> mouse returns NO ortholog (verified via Ensembl REST:
#     0 edges). eggNOG groups it correctly. Without the edge, lemur RAMP2 has
#     no ortholog and detect_substitutions drops it — losing the conserved
#     anchor and blocking triad closure (validation Section 3).
#
# Run from repo root inside the r-paralogswap env:
#     conda activate r-paralogswap
#     Rscript data-raw/run_lemur_mouse.R
# =====================================================================

set.seed(1)
stopifnot(dir.exists("R"))
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# =====================================================================
# CONFIG
# =====================================================================
LEMUR <- "mmurinus"    # species_a / anchor
MOUSE <- "mmusculus"   # species_b

ENSEMBL_RELEASE  <- 116
MIRRORS          <- c("useast", "asia", NA)  # cycle hosts; NA = default (www). Mirror
                                 # attempts drop the release pin (biomaRt forbids
                                 # version+mirror); current release == 116 anyway.
RESUME           <- TRUE
ID_TYPE          <- "symbol"     # atlases are symbol-keyed; lemur UPPER (RAMP1),
                                 # mouse title-case (Ramp1). Graph + focal follow suit.
MIN_CONFIDENCE   <- "low"
MIN_PERC_ID      <- NULL

NORMALIZATION <- "lognorm"
RESOLUTION    <- 0.8
N_PCS         <- 30
N_VAR_GENES   <- 2000
ASSAY         <- "RNA"
SLOT          <- "data"
COR_METHOD    <- "spearman"

GAMMA         <- 20
K_KNN         <- 5
MIN_CELLS     <- 10
MUTUAL_K      <- 5
MIN_PAIRS     <- 10
MIN_CORRESPONDENCE <- 0.3
DELTA_THRESHOLD    <- 0.3

MAKE_PLOTS      <- TRUE
FIG_DIR         <- "tests/manual/plots"        # git-ignored
SHOW_UNCORRECTED <- TRUE                        # demonstrate the RAMP2 gap first (Section 3)

# --- Focal genes. Anchored on lemur (species_a); mouse symbols title-case.
FOCAL_LEMUR <- c(RAMP1 = "RAMP1", RAMP2 = "RAMP2", RAMP3 = "RAMP3")
FOCAL_MOUSE <- c(RAMP1 = "Ramp1", RAMP2 = "Ramp2", RAMP3 = "Ramp3")

# --- The correction (VERIFIED needed via Ensembl REST r116: lemur RAMP2
#     ENSMICG00000028500 -> mouse returns 0 ortholog edges; mouse Ramp2 =
#     ENSMUSG00000001240, chr11). Symbol-keyed to match the symbol-space graph.
CORRECTIONS <- data.frame(
  gene_a        = "RAMP2",   # lemur
  gene_b        = "Ramp2",   # mouse
  species_a     = LEMUR,
  species_b     = MOUSE,
  relationship  = "ortholog",
  ortholog_type = "one2one",
  source        = "eggnog",
  stringsAsFactors = FALSE
)

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
s3_exists <- function(key) {
  st <- suppressWarnings(system2("aws", c("s3", "ls", key), stdout = TRUE, stderr = TRUE))
  isTRUE(length(st) > 0) && !any(grepl("error|not *found", st, ignore.case = TRUE))
}
resume_or <- function(key, build_fn) {
  if (isTRUE(RESUME) && s3_exists(key)) { message("  RESUME: loading ", key); return(s3_read_rds(key)) }
  obj <- build_fn(); s3_write_rds(obj, key); obj
}
banner <- function(...) message("\n", strrep("=", 62), "\n", ..., "\n", strrep("=", 62))
soft_check <- function(label, value, expected, tol = 0.15) {
  ok <- is.numeric(value) && !is.na(value) && abs(value - expected) <= tol * expected
  message(sprintf("  [%s] %s: got %s (expected ~%s)", if (ok) "ok" else "DRIFT", label, value, expected))
  invisible(ok)
}
# mirror-cycling homology build; drops the release pin on mirror attempts
build_graph_resilient <- function(..., mirrors = MIRRORS, tries_per = 2, wait = 8) {
  args <- list(...); pinned <- args$ensembl_release; last <- NULL; n <- 0L
  for (m in mirrors) for (t in seq_len(tries_per)) {
    n <- n + 1L; a <- args
    if (is.na(m)) { a$mirror <- NULL; a$ensembl_release <- pinned }
    else          { a$mirror <- m;    a$ensembl_release <- NULL }
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
                    conditionMessage(out))); Sys.sleep(wait)
  }
  stop(last)
}

# load lemur atlas early (bounds the Stage 1a pull; reused if lemur re-clustered)
lemur_raw <- readRDS("data-raw/lemur_lung.rds")
lemur_genes <- rownames(lemur_raw)
message("lemur atlas genes available to bridge: ", length(lemur_genes))

# =====================================================================
# STAGE 1 — homology graphs (broad orthologs + focal RAMP w/ correction)
# =====================================================================
banner("STAGE 1a  build_homology_graph  (atlas-restricted orthologs, ", LEMUR, " <-> ", MOUSE, ")")
hg_broad <- resume_or(file.path(S3_HOMOLOGY, "lemur_mouse_hg_broad.rds"), function()
  build_graph_resilient(
    species_a = LEMUR, species_b = MOUSE, genes = lemur_genes,
    ensembl_release = ENSEMBL_RELEASE, id_type = ID_TYPE, corrections = NULL,
    include_paralogs = FALSE, min_perc_id = MIN_PERC_ID,
    min_confidence = MIN_CONFIDENCE, verbose = TRUE))

# --- Section 3 narrative: show the RAMP2 gap BEFORE correcting it.
if (isTRUE(SHOW_UNCORRECTED)) {
  banner("STAGE 1b(i)  focal graph WITHOUT correction  (expect RAMP2 gap)")
  hg_focal_raw <- build_graph_resilient(
    species_a = LEMUR, species_b = MOUSE, genes = unname(FOCAL_LEMUR),
    ensembl_release = ENSEMBL_RELEASE, id_type = ID_TYPE, corrections = NULL,
    include_paralogs = TRUE, min_perc_id = MIN_PERC_ID,
    min_confidence = MIN_CONFIDENCE, verbose = TRUE)
  .r2_raw <- hg_focal_raw[hg_focal_raw$gene_a == "RAMP2" &
                          hg_focal_raw$relationship == "ortholog", ]
  message("  lemur RAMP2 ortholog edges WITHOUT correction: ", nrow(.r2_raw),
          if (nrow(.r2_raw) == 0) "  <-- reproduces the Compara mouse-Ramp2 gap" else "")
}

banner("STAGE 1b(ii)  focal graph WITH eggNOG correction")
hg_focal <- resume_or(file.path(S3_HOMOLOGY, "lemur_mouse_hg_focal.rds"), function()
  build_graph_resilient(
    species_a = LEMUR, species_b = MOUSE, genes = unname(FOCAL_LEMUR),
    ensembl_release = ENSEMBL_RELEASE, id_type = ID_TYPE, corrections = CORRECTIONS,
    include_paralogs = TRUE, min_perc_id = MIN_PERC_ID,
    min_confidence = MIN_CONFIDENCE, verbose = TRUE))

# correction landed? lemur RAMP2 <-> mouse Ramp2 ortholog edge, source=eggnog
.r2_fix <- hg_focal[hg_focal$gene_a == "RAMP2" & hg_focal$gene_b == "Ramp2" &
                    hg_focal$relationship == "ortholog", ]
message("  lemur RAMP2 -> mouse Ramp2 ortholog edge after correction: ", nrow(.r2_fix),
        if (nrow(.r2_fix) >= 1) paste0(" (source=", paste(unique(.r2_fix$source), collapse=","), ")")
        else "  <-- WARNING: correction did not merge")
.para_mouse <- hg_focal[hg_focal$relationship == "paralog" & hg_focal$species_a == MOUSE, ]
message("  mouse-side RAMP paralog edges in hg_focal: ", nrow(.para_mouse))
stopifnot("correction must place lemur RAMP2 -> mouse Ramp2" = nrow(.r2_fix) >= 1)

# =====================================================================
# STAGE 2 — clustering (per-species RESUME: lemur is already in S3)
# =====================================================================
banner("STAGE 2  cluster_species")
get_clustered <- function(sp_label, raw_getter) {
  key <- file.path(S3_RESULTS, paste0(sp_label, "_lung_clustered.rds"))
  if (isTRUE(RESUME) && s3_exists(key)) { message("  RESUME: ", sp_label, " clustered from S3"); return(s3_read_rds(key)) }
  raw <- raw_getter()
  scl <- cluster_species(setNames(list(raw), sp_label), species_col = "species",
                         assay = ASSAY, normalization = NORMALIZATION, resolution = RESOLUTION,
                         n_pcs = N_PCS, n_var_genes = N_VAR_GENES, verbose = TRUE)[[sp_label]]
  s3_write_rds(scl, key); scl
}
sc_lemur <- get_clustered("lemur", function() lemur_raw)  # already clustered in the human arm
sc_mouse <- get_clustered("mouse", function() s3_read_rds(file.path(S3_RESULTS, "mouse_lung.rds")))

# =====================================================================
# STAGE 3 — cluster correspondence (uses hg_broad)
# =====================================================================
banner("STAGE 3  match_clusters")
correspondence <- match_clusters(
  clusters_a = sc_lemur, clusters_b = sc_mouse, homology_graph = hg_broad,
  species_a = LEMUR, species_b = MOUSE, cor_method = COR_METHOD,
  min_correspondence = MIN_CORRESPONDENCE, assay = ASSAY, slot = SLOT, verbose = TRUE)
s3_write_rds(correspondence, file.path(S3_RESULTS, "lemur_mouse_cluster_correspondence.rds"))
cat("\n--- summary(cluster_correspondence) ---\n"); print(summary(correspondence))

# =====================================================================
# STAGE 4 — metacells per species
# =====================================================================
banner("STAGE 4  build_metacells  (gamma=", GAMMA, ")")
mc_lemur <- build_metacells(sc_lemur, gamma = GAMMA, assay = ASSAY, slot = SLOT,
                            n_var_genes = N_VAR_GENES, k_knn = K_KNN,
                            min_cells = MIN_CELLS, cluster_col = NULL, verbose = TRUE)
mc_mouse <- build_metacells(sc_mouse, gamma = GAMMA, assay = ASSAY, slot = SLOT,
                            n_var_genes = N_VAR_GENES, k_knn = K_KNN,
                            min_cells = MIN_CELLS, cluster_col = NULL, verbose = TRUE)
s3_write_rds(mc_lemur, file.path(S3_RESULTS, "lemur_lung_metacells.rds"))       # same as human arm
s3_write_rds(mc_mouse, file.path(S3_RESULTS, "mouse_lung_metacells.rds"))
soft_check("lemur metacells", ncol(mc_lemur$data), 1704)   # same lemur atlas -> should reproduce
message("  mouse metacells: ", ncol(mc_mouse$data), " (no prior checkpoint for this arm)")

# =====================================================================
# STAGE 5 — fine metacell correspondence (mutual-kNN; uses hg_broad)
# =====================================================================
banner("STAGE 5  match_metacells  (mutual_k=", MUTUAL_K, ")")
metacell_correspondence <- match_metacells(
  metacells_a = mc_lemur, metacells_b = mc_mouse,
  correspondence = correspondence, homology_graph = hg_broad,
  species_a = LEMUR, species_b = MOUSE, cor_method = COR_METHOD,
  mutual_k = MUTUAL_K, verbose = TRUE)
s3_write_rds(metacell_correspondence, file.path(S3_RESULTS, "lemur_mouse_metacell_correspondence.rds"))
message("  metacell pairs: ", nrow(metacell_correspondence), " (no prior checkpoint for this arm)")

# =====================================================================
# STAGE 6 — correlations (uses hg_focal, corrected) + detection
# =====================================================================
banner("STAGE 6  compute_homolog_correlations + detect_substitutions")
homolog_correlations <- compute_homolog_correlations(
  metacells_a = mc_lemur, metacells_b = mc_mouse,
  metacell_correspondence = metacell_correspondence, homology_graph = hg_focal,
  species_a = LEMUR, species_b = MOUSE, focal_genes = unname(FOCAL_LEMUR),
  cor_method = COR_METHOD, min_pairs = MIN_PAIRS, verbose = TRUE)
s3_write_rds(homolog_correlations, file.path(S3_RESULTS, "lemur_mouse_homolog_correlations.rds"))

substitutions <- detect_substitutions(homolog_correlations,
                                      delta_threshold = DELTA_THRESHOLD, require_ortholog = TRUE)
s3_write_rds(substitutions, file.path(S3_RESULTS, "lemur_mouse_substitutions.rds"))

# =====================================================================
# REPRODUCTION GATE — lemur RAMP1 substitutes to mouse Ramp2; RAMP2 conserved
#   No published r-values for this arm, so assert the PATTERN, not magnitudes.
# =====================================================================
banner("REPRODUCTION GATE  (lemur RAMP1 -> mouse Ramp2; RAMP2 conserved via correction)")
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
  "RAMP1 present in substitutions"               = nrow(r1) == 1,
  "RAMP1 must be flagged"                         = isTRUE(r1$flagged[1]),
  "RAMP1 best paralog is mouse Ramp2"             = identical(as.character(r1$best_paralog[1]),
                                                              FOCAL_MOUSE[["RAMP2"]]),
  "RAMP1 r_paralog exceeds r_ortholog by >0.30"   = (r1$r_paralog[1] - r1$r_ortholog[1]) > 0.30,
  "RAMP2 present (correction worked)"             = nrow(r2) == 1,
  "RAMP2 must NOT be flagged (conserved)"         = !isTRUE(r2$flagged[1])
)
banner("GATE PASSED  \u2014 lemur RAMP1 substitution + RAMP2 conservation reproduced (mouse arm).")

# =====================================================================
# OPTIONAL — plots (git-ignored). graining curve sweeps gammas; gamma=5 is the
#   multi-hour pass (thousands of metacells x per-metacell correlation), dropped.
# =====================================================================
if (isTRUE(MAKE_PLOTS)) {
  banner("PLOTS -> ", FIG_DIR, "/")
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  save_png <- function(p, f, w = 7, h = 5) {
    ggplot2::ggsave(file.path(FIG_DIR, f), p, width = w, height = h, dpi = 150)
    message("  wrote ", file.path(FIG_DIR, f))
  }
  tryCatch(save_png(plot_substitution_elbow(substitutions),         "lemur_mouse_elbow.png"),
           error = function(e) message("  elbow skipped: ", conditionMessage(e)))
  tryCatch(save_png(plot_ortholog_vs_paralog(homolog_correlations), "lemur_mouse_ortholog_vs_paralog.png"),
           error = function(e) message("  ovp skipped: ", conditionMessage(e)))
  tryCatch(save_png(plot_graining_curve(sc_mouse, gammas = c(10, 20, 50)), "mouse_graining_curve.png"),
           error = function(e) message("  graining skipped: ", conditionMessage(e)))
}

banner("DONE  \u2014 lemur-mouse intermediates persisted under ", S3_RESULTS)
