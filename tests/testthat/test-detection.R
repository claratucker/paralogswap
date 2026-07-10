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
  # Derive the fixture's own statistics rather than hard-coding a threshold
  # tuned to whatever the fixture happened to produce. gX has a strongly
  # negative ortholog correlation and a positive paralog correlation, so
  # delta_z substantially exceeds delta_r.
  base <- detect_substitutions(hc, delta_threshold = -Inf, delta_scale = "r")
  d_r <- base$delta_r[1]; d_z <- base$delta_z[1]
  expect_gt(d_z, d_r)   # the transform widens the gap when r_ortholog < 0

  # r scale: flags iff delta_r clears the threshold
  expect_false(detect_substitutions(hc, delta_threshold = d_r + 0.1,
                                    delta_scale = "r")$flagged[1])
  expect_true(detect_substitutions(hc, delta_threshold = d_r - 0.1,
                                   delta_scale = "r")$flagged[1])

  # z scale: same, against delta_z
  expect_false(detect_substitutions(hc, delta_threshold = d_z + 0.1,
                                    delta_scale = "z")$flagged[1])
  expect_true(detect_substitutions(hc, delta_threshold = d_z - 0.1,
                                   delta_scale = "z")$flagged[1])

  # flagged must follow the scale it was thresholded on, never the other column
  t_mid <- (d_r + d_z) / 2   # lies between the two statistics
  expect_true(detect_substitutions(hc, delta_threshold = t_mid,
                                   delta_scale = "z")$flagged[1])
  expect_false(detect_substitutions(hc, delta_threshold = t_mid,
                                    delta_scale = "r")$flagged[1])
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


# ---- attr(x, "dropped"): why a focal gene could not be scored --------------
# Each fixture isolates one branch of the reason chain. Ordering matters: a gene
# whose every edge is NA must read focal_gene_invariant, not ortholog_invariant.

hc_fixture <- function(rows, min_pairs = 10, cor_method = "spearman") {
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "params") <- list(min_pairs = min_pairs, cor_method = cor_method)
  class(out) <- c("homolog_correlations", "data.frame")
  out
}
edge <- function(gene_a, gene_b, relationship, r, n_pairs = 100) {
  data.frame(gene_a = gene_a, gene_b = gene_b, relationship = relationship,
             r = r, n_pairs = n_pairs, stringsAsFactors = FALSE)
}
reason_for <- function(hc, gene) {
  d <- attr(detect_substitutions(hc), "dropped")
  d$reason[d$gene == gene]
}

test_that("a flat focal gene reads focal_gene_invariant, not ortholog_invariant", {
  # every edge undefined at full coverage: the species-A gene has no variance
  hc <- hc_fixture(list(
    edge("FLAT", "flat_orth", "ortholog", NA_real_),
    edge("FLAT", "para1",     "paralog",  NA_real_),
    edge("FLAT", "para2",     "paralog",  NA_real_)))
  expect_identical(reason_for(hc, "FLAT"), "focal_gene_invariant")

  d <- attr(detect_substitutions(hc), "dropped")
  expect_true(is.na(d$best_paralog[d$gene == "FLAT"]))   # the discriminator
})

test_that("a flat ortholog with live paralogs reads ortholog_invariant", {
  # full coverage, ortholog undefined, but a paralog correlates: the ortholog is
  # silenced in species B -- a substitution delta_r cannot score
  hc <- hc_fixture(list(
    edge("GENE", "orth",  "ortholog", NA_real_),
    edge("GENE", "para1", "paralog",  0.49),
    edge("GENE", "para2", "paralog",  0.41)))
  expect_identical(reason_for(hc, "GENE"), "ortholog_invariant")

  d <- attr(detect_substitutions(hc), "dropped")
  expect_identical(d$best_paralog[d$gene == "GENE"], "para1")
  expect_equal(d$r_paralog[d$gene == "GENE"], 0.49)
})

test_that("a sparse ortholog reads ortholog_too_sparse", {
  hc <- hc_fixture(list(
    edge("GENE", "orth",  "ortholog", NA_real_, n_pairs = 3),   # below min_pairs
    edge("GENE", "para1", "paralog",  0.5)))
  expect_identical(reason_for(hc, "GENE"), "ortholog_too_sparse")
})

test_that("all-NA edges below min_pairs prefer sparsity over invariance", {
  # focal_flat requires FULL coverage; sparse all-NA must not claim invariance
  hc <- hc_fixture(list(
    edge("GENE", "orth",  "ortholog", NA_real_, n_pairs = 2),
    edge("GENE", "para1", "paralog",  NA_real_, n_pairs = 2)))
  expect_identical(reason_for(hc, "GENE"), "ortholog_too_sparse")
})

test_that("a focal gene with no ortholog edge reads no_ortholog_edge", {
  hc <- hc_fixture(list(edge("GENE", "para1", "paralog", 0.6)))
  expect_identical(reason_for(hc, "GENE"), "no_ortholog_edge")
})

test_that("scoreable genes are absent from the dropped table", {
  hc <- hc_fixture(list(
    edge("OK",   "orth",  "ortholog", 0.1),
    edge("OK",   "para1", "paralog",  0.7),
    edge("FLAT", "orth2", "ortholog", NA_real_),
    edge("FLAT", "para2", "paralog",  NA_real_)))
  subs <- detect_substitutions(hc)
  expect_true("OK" %in% subs$gene)
  expect_false("FLAT" %in% subs$gene)
  expect_identical(attr(subs, "dropped")$gene, "FLAT")
  expect_true(subs$flagged[subs$gene == "OK"])           # delta_r = 0.6 > 0.3
})

test_that("params are absent on older objects without breaking the reasons", {
  hc <- hc_fixture(list(
    edge("FLAT", "orth",  "ortholog", NA_real_),
    edge("FLAT", "para1", "paralog",  NA_real_)))
  attr(hc, "params") <- NULL                              # pre-backfill object
  expect_identical(reason_for(hc, "FLAT"), "focal_gene_invariant")
})
