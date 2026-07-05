# Two toy species with a shared ortholog bridge and matching cluster structure.
# Species A clusters 0/1 should correspond to species B clusters 0/1 via orthologs.
make_matched_pair <- function(seed = 1) {
  set.seed(seed)
  ng <- 100; nc <- 80
  # shared "ortholog" genes g1..g100 (A) <-> h1..h100 (B), same expression program
  prog <- matrix(rpois(ng * nc, 2), ng, nc)
  prog[1:20, 1:40]  <- prog[1:20, 1:40]  + 12   # cluster-0 markers
  prog[21:40, 41:80] <- prog[21:40, 41:80] + 12 # cluster-1 markers

  mkA <- { m <- prog; rownames(m) <- paste0("g", 1:ng); colnames(m) <- paste0("a", 1:nc); m }
  mkB <- { m <- prog; rownames(m) <- paste0("h", 1:ng); colnames(m) <- paste0("b", 1:nc); m }

  oa <- Seurat::CreateSeuratObject(Matrix::Matrix(mkA, sparse = TRUE))
  ob <- Seurat::CreateSeuratObject(Matrix::Matrix(mkB, sparse = TRUE))
  list(a = oa, b = ob)
}

toy_graph <- function(ng = 100) {
  as_homology_graph(data.frame(
    gene_a = paste0("g", 1:ng), gene_b = paste0("h", 1:ng),
    species_a = "spA", species_b = "spB",
    relationship = "ortholog", ortholog_type = "one2one",
    stringsAsFactors = FALSE
  ), verbose = FALSE)
}

cluster_toys <- function(pair) {
  sc <- cluster_species(list(a = pair$a, b = pair$b),
                        n_pcs = 10, n_var_genes = 50, verbose = FALSE)
  sc
}

test_that("match_clusters returns cluster_correspondence with reciprocal matches", {
  pair <- make_matched_pair(); sc <- cluster_toys(pair); hg <- toy_graph()
  cc <- match_clusters(sc$a, sc$b, hg, species_a = "spA", species_b = "spB",
                       verbose = FALSE)
  expect_s3_class(cc, "cluster_correspondence")
  expect_true(all(c("cluster_a","cluster_b","score","reciprocal") %in% names(cc)))
  expect_true(any(cc$reciprocal))
})

test_that("correlation matrix is attached and sized correctly", {
  pair <- make_matched_pair(); sc <- cluster_toys(pair); hg <- toy_graph()
  cc <- match_clusters(sc$a, sc$b, hg, "spA", "spB", verbose = FALSE)
  cm <- attr(cc, "cor_matrix")
  expect_true(is.matrix(cm))
  expect_equal(nrow(cm), length(levels(Seurat::Idents(sc$a))))
})

test_that("min_correspondence gates reciprocal flag", {
  pair <- make_matched_pair(); sc <- cluster_toys(pair); hg <- toy_graph()
  cc_lo <- match_clusters(sc$a, sc$b, hg, "spA","spB",
                          min_correspondence = 0, verbose = FALSE)
  cc_hi <- match_clusters(sc$a, sc$b, hg, "spA","spB",
                          min_correspondence = 0.999, verbose = FALSE)
  expect_gte(sum(cc_lo$reciprocal), sum(cc_hi$reciprocal))
})

test_that("errors when no one-to-one orthologs bridge the objects", {
  pair <- make_matched_pair(); sc <- cluster_toys(pair)
  bad <- as_homology_graph(data.frame(
    gene_a = "x1", gene_b = "y1", species_a = "spA", species_b = "spB",
    relationship = "ortholog", ortholog_type = "one2one"), verbose = FALSE)
  expect_error(
    match_clusters(sc$a, sc$b, bad, "spA","spB", verbose = FALSE),
    "No one-to-one orthologs bridge"
  )
})

test_that("graph orientation flip works when species are reversed", {
  pair <- make_matched_pair(); sc <- cluster_toys(pair); hg <- toy_graph()
  # call with species reversed relative to the graph -> should flip internally
  cc <- match_clusters(sc$b, sc$a, hg, species_a = "spB", species_b = "spA",
                       verbose = FALSE)
  expect_s3_class(cc, "cluster_correspondence")
  expect_true(any(cc$reciprocal))
})
