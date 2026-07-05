# Build a metacell grid + homology graph with a PLANTED substitution:
#   species A gene "gX" should track species B "gY_para" (paralog), NOT its
#   ortholog "gY_orth". So detect_substitutions must flag gX.

make_detection_fixture <- function(n_pairs = 40, seed = 1) {
  set.seed(seed)
  # shared latent pattern across the grid (e.g. an on/off compartment axis)
  pattern <- c(rep(3, n_pairs/2), rep(0, n_pairs/2)) + rnorm(n_pairs, 0, 0.3)

  # species A metacells
  genesA <- c("gX", "gA_orth_of_gY", "noise1", "noise2")
  Da <- rbind(
    gX               = pattern,                      # gX follows the pattern
    gA_orth_of_gY    = rev(pattern) + rnorm(n_pairs),# unrelated to pattern
    noise1           = rnorm(n_pairs),
    noise2           = rnorm(n_pairs)
  )
  colnames(Da) <- paste0("mcA", seq_len(n_pairs))

  # species B metacells
  Db <- rbind(
    gY_orth  = rev(pattern) + rnorm(n_pairs),   # gX's ORTHOLOG: does NOT track gX
    gY_para  = pattern + rnorm(n_pairs, 0, 0.3),# gX's PARALOG: DOES track gX
    noise1   = rnorm(n_pairs),
    noise2   = rnorm(n_pairs)
  )
  colnames(Db) <- paste0("mcB", seq_len(n_pairs))

  mca <- list(data = Da,
              meta = data.frame(metacell_id = colnames(Da),
                                cluster = "0", n_cells = 20,
                                stringsAsFactors = FALSE))
  mcb <- list(data = Db,
              meta = data.frame(metacell_id = colnames(Db),
                                cluster = "0", n_cells = 20,
                                stringsAsFactors = FALSE))
  class(mca) <- c("metacells","list"); class(mcb) <- c("metacells","list")

  # matched grid: pair mcA_i <-> mcB_i
  mmc <- data.frame(
    metacell_a = colnames(Da), metacell_b = colnames(Db),
    cluster_a = "0", cluster_b = "0",
    score = 0.7, stringsAsFactors = FALSE)
  class(mmc) <- c("metacell_correspondence","data.frame")

  # homology graph: gX(A) <-> gY_orth(B) is the ortholog;
  #   gY_orth <-> gY_para are within-B paralogs.
  hg <- as_homology_graph(rbind(
    data.frame(gene_a="gX", gene_b="gY_orth", species_a="spA", species_b="spB",
               relationship="ortholog", ortholog_type="one2one"),
    data.frame(gene_a="gY_orth", gene_b="gY_para", species_a="spB", species_b="spB",
               relationship="paralog", ortholog_type=NA)
  ), verbose = FALSE)

  list(mca = mca, mcb = mcb, mmc = mmc, hg = hg)
}

test_that("compute_homolog_correlations recovers ortholog and paralog r", {
  fx <- make_detection_fixture()
  hc <- compute_homolog_correlations(
    fx$mca, fx$mcb, fx$mmc, fx$hg,
    species_a = "spA", species_b = "spB",
    focal_genes = "gX", verbose = FALSE)
  expect_s3_class(hc, "homolog_correlations")
  # ortholog (gX vs gY_orth) should be LOW/negative; paralog (gX vs gY_para) HIGH
  r_orth <- hc$r[hc$relationship == "ortholog"]
  r_para <- hc$r[hc$relationship == "paralog"]
  expect_lt(r_orth, r_para)          # paralog tracks better
  expect_gt(r_para, 0.5)             # paralog genuinely high
})

test_that("detect_substitutions flags a planted substitution", {
  fx <- make_detection_fixture()
  hc <- compute_homolog_correlations(fx$mca, fx$mcb, fx$mmc, fx$hg,
                                     species_a="spA", species_b="spB",
                                     focal_genes="gX", verbose=FALSE)
  subs <- detect_substitutions(hc, delta_threshold = 0.3)
  expect_s3_class(subs, "substitutions")
  gx <- subs[subs$gene == "gX", ]
  expect_true(gx$flagged)
  expect_equal(gx$best_paralog, "gY_para")
  expect_gt(gx$delta_r, 0.3)
})

test_that("detect_substitutions does NOT flag a conserved gene", {
  # conserved: ortholog tracks best, no paralog beats it
  set.seed(2); n <- 40
  pattern <- c(rep(3, n/2), rep(0, n/2)) + rnorm(n, 0, 0.3)
  Da <- rbind(gC = pattern); colnames(Da) <- paste0("mcA", 1:n)
  Db <- rbind(gC_orth = pattern + rnorm(n,0,0.3),      # ortholog tracks well
              gC_para = rnorm(n))                       # paralog is noise
  colnames(Db) <- paste0("mcB", 1:n)
  mca <- structure(list(data=Da, meta=data.frame(metacell_id=colnames(Da),
              cluster="0", n_cells=20)), class=c("metacells","list"))
  mcb <- structure(list(data=Db, meta=data.frame(metacell_id=colnames(Db),
              cluster="0", n_cells=20)), class=c("metacells","list"))
  mmc <- structure(data.frame(metacell_a=colnames(Da), metacell_b=colnames(Db),
              cluster_a="0", cluster_b="0", score=0.7),
              class=c("metacell_correspondence","data.frame"))
  hg <- as_homology_graph(rbind(
    data.frame(gene_a="gC", gene_b="gC_orth", species_a="spA", species_b="spB",
               relationship="ortholog", ortholog_type="one2one"),
    data.frame(gene_a="gC_orth", gene_b="gC_para", species_a="spB", species_b="spB",
               relationship="paralog", ortholog_type=NA)), verbose=FALSE)
  hc <- compute_homolog_correlations(mca, mcb, mmc, hg, "spA","spB",
                                     focal_genes="gC", verbose=FALSE)
  subs <- detect_substitutions(hc, delta_threshold = 0.3)
  expect_false(subs$flagged[subs$gene == "gC"])
})

test_that("delta_threshold controls flagging", {
  fx <- make_detection_fixture()
  hc <- compute_homolog_correlations(fx$mca, fx$mcb, fx$mmc, fx$hg,
                                     "spA","spB", focal_genes="gX", verbose=FALSE)
  # gX has delta_r ~1.49 (r_para ~0.87 vs r_orth ~ -0.62)
  expect_false(detect_substitutions(hc, delta_threshold = 1.6)$flagged[1])
  expect_true(detect_substitutions(hc, delta_threshold = 1.0)$flagged[1])
})

test_that("min_pairs guard errors on too-small grid", {
  fx <- make_detection_fixture(n_pairs = 6)   # below default min_pairs=10
  expect_error(
    compute_homolog_correlations(fx$mca, fx$mcb, fx$mmc, fx$hg, "spA","spB",
                                 focal_genes="gX", verbose=FALSE),
    "Too few matched metacell pairs"
  )
})

test_that("NA-ortholog gene is dropped when require_ortholog=TRUE", {
  # gene whose ortholog has zero variance -> NA correlation -> dropped
  fx <- make_detection_fixture()
  fx$mcb$data["gY_orth", ] <- 5   # constant -> ortholog r = NA
  hc <- compute_homolog_correlations(fx$mca, fx$mcb, fx$mmc, fx$hg,
                                     "spA","spB", focal_genes="gX", verbose=FALSE)
  subs <- detect_substitutions(hc, require_ortholog = TRUE)
  expect_equal(nrow(subs[subs$gene=="gX", ]), 0)  # dropped: no valid ortholog
})

