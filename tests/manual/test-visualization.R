# tests/manual/test-visualization.R
# Manual/interactive test of Stage 9 plots. Not a testthat file — these render
# to disk so you can eyeball them. Run from the package root after
# devtools::load_all() (so plot_* and .within_metacell_cor are available).

#library(paralogswap)   # or devtools::load_all(".")
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
outdir <- "tests/manual/plots"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(g, name, w = 6, h = 5)
  ggplot2::ggsave(file.path(outdir, name), g, width = w, height = h, dpi = 120)

# ---------------------------------------------------------------------------
# PATH A — real RAMP result (preferred). Set these to your objects in session.
#   hc_lh   : homolog_correlations from the lemur–human run
#   subs_lh : detect_substitutions(hc_lh, ...)
#   lemur   : the clustered lemur Seurat object (or cluster_species result)
# ---------------------------------------------------------------------------
have_real <- all(sapply(c("hc_lh", "subs_lh"), exists))

if (have_real) {
  message("Testing against real lemur–human objects.")

  g1 <- plot_substitution_elbow(subs_lh)                     # default threshold from attr
  save_plot(g1, "01_substitution_elbow.png")

  g2 <- plot_ortholog_vs_paralog(hc_lh, label_genes = c("RAMP1", "RAMP2", "RAMP3"))
  save_plot(g2, "02_ortholog_vs_paralog.png")

  # sanity checks: does the picture match the known result?
  df <- as.data.frame(subs_lh)
  ramp1 <- df[df$gene == "RAMP1", ]
  cat("\nRAMP1 row:\n"); print(ramp1)
  stopifnot(nrow(ramp1) == 1)
  cat("RAMP1 flagged:", ramp1$flagged,
      "| delta_r:", round(ramp1$delta_r, 3), "\n")
  cat("RAMP1 sits above diagonal:", ramp1$r_paralog > ramp1$r_ortholog, "\n")

  # threshold echo: plot line should follow whatever detect_substitutions used
  cat("stored threshold:", attr(subs_lh, "params")$delta_threshold, "\n")

} else {
  message("Real objects not found — using synthetic fixture for the 2 corr plots.")

  # planted substitution: gX's paralog tracks, ortholog doesn't
  set.seed(1)
  mk <- function(gene, r_o, r_p) data.frame(
    gene = gene, ortholog = paste0(gene, "_orth"),
    best_paralog = paste0(gene, "_para"),
    r_ortholog = r_o, r_paralog = r_p, delta_r = r_p - r_o,
    flagged = (r_p - r_o) >= 0.3, stringsAsFactors = FALSE)
  subs <- do.call(rbind, list(
    mk("RAMP1", -0.06, 0.54),   # the planted substitution
    mk("gA",     0.62, 0.30),   # conserved
    mk("gB",     0.40, 0.55),   # mild candidate
    mk("gC",     0.10, 0.12),   # noise
    mk("gD",     0.70, 0.20)))
  attr(subs, "params") <- list(delta_threshold = 0.3)
  class(subs) <- c("substitutions", "data.frame")

  # homolog_correlations shape: gene_a, gene_b, relationship, r, n_pairs
  hc <- do.call(rbind, lapply(seq_len(nrow(subs)), function(i) {
    s <- subs[i, ]
    rbind(
      data.frame(gene_a = s$gene, gene_b = s$ortholog,
                 relationship = "ortholog", r = s$r_ortholog, n_pairs = 40),
      data.frame(gene_a = s$gene, gene_b = s$best_paralog,
                 relationship = "paralog",  r = s$r_paralog,  n_pairs = 40))
  }))
  class(hc) <- c("homolog_correlations", "data.frame")

  g1 <- plot_substitution_elbow(subs)
  save_plot(g1, "01_substitution_elbow.png")

  g2 <- plot_ortholog_vs_paralog(hc, label_genes = c("RAMP1", "gB"))
  save_plot(g2, "02_ortholog_vs_paralog.png")

  # preview a different threshold without recomputing flags
  g2b <- plot_ortholog_vs_paralog(hc, delta_threshold = 0.5)
  save_plot(g2b, "02b_ortholog_vs_paralog_thresh0.5.png")

  cat("\nFixture: RAMP1 delta_r =", subs$delta_r[subs$gene == "RAMP1"],
      "(should flag at 0.3, not at 0.99)\n")
}

# ---------------------------------------------------------------------------
# PATH B — graining curve. Needs a clustered Seurat object + SuperCell.
#   Set `lemur` to your clustered lemur object (or cluster_species result).
#   This actually runs build_metacells at each gamma, so it's slow-ish.
# ---------------------------------------------------------------------------
if (exists("lemur")) {
  message("Testing plot_graining_curve on `lemur` — this grains at each gamma, be patient.")
  g3 <- plot_graining_curve(lemur, gammas = c(5, 10, 20, 50, 100), verbose = TRUE)
  save_plot(g3, "03_graining_curve.png", w = 6.5, h = 4.5)

  # quick internal check on the helper alone at one gamma
  mc <- build_metacells(lemur, gamma = 20, verbose = FALSE)
  dm <- Seurat::GetAssayData(.unwrap_clustered(lemur), assay = "RNA", layer = "data")
  coh <- paralogswap:::.within_metacell_cor(dm, mc$membership)
  cat("within-metacell coherence at gamma=20: mean",
      round(mean(coh, na.rm = TRUE), 3),
      "| NA (1-cell) metacells:", sum(is.na(coh)), "\n")
} else {
  message("`lemur` not in session — skipping graining curve. ",
          "Set it to your clustered lemur object to test plot_graining_curve.")
}

message("\nDone. Plots written to ", normalizePath(outdir))
