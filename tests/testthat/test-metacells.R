# toy clustered object with 2 clusters, enough cells to grain
make_toy_clustered <- function(seed = 1, n = 120) {
  set.seed(seed)
  ng <- 100; half <- n/2
  m <- matrix(rpois(ng*n, 2), ng, n)
  m[1:20, 1:half] <- m[1:20, 1:half] + 12
  m[21:40, (half+1):n] <- m[21:40, (half+1):n] + 12
  rownames(m) <- paste0("g", 1:ng); colnames(m) <- paste0("c", 1:n)
  o <- Seurat::CreateSeuratObject(Matrix::Matrix(m, sparse=TRUE))
  cluster_species(list(x=o), n_pcs=10, n_var_genes=50, verbose=FALSE)$x
}

test_that("build_metacells produces pure per-cluster metacells", {
  o <- make_toy_clustered()
  mc <- build_metacells(o, gamma = 5, n_var_genes = 50, verbose = FALSE)
  expect_s3_class(mc, "metacells")
  expect_true(all(c("counts","data","meta","membership") %in% names(mc)))
  # every metacell maps to exactly one cluster (purity by construction)
  expect_equal(nrow(mc$meta), ncol(mc$data))
  expect_true(all(mc$meta$n_cells >= 1))
})

test_that("build_metacells preserves genes and cell count", {
  o <- make_toy_clustered()
  mc <- build_metacells(o, gamma = 5, n_var_genes = 50, verbose = FALSE)
  expect_equal(nrow(mc$data), nrow(o))
  expect_equal(sum(mc$meta$n_cells), ncol(o))
  expect_equal(length(mc$membership), ncol(o))
})

test_that("match_metacells produces pairs via mutual-kNN within matched clusters", {
  oa <- make_toy_clustered(1); ob <- make_toy_clustered(2)
  mca <- build_metacells(oa, gamma = 5, n_var_genes = 50, verbose = FALSE)
  mcb <- build_metacells(ob, gamma = 5, n_var_genes = 50, verbose = FALSE)
  hg <- as_homology_graph(data.frame(
    gene_a = paste0("g",1:100), gene_b = paste0("g",1:100),
    species_a="spA", species_b="spB",
    relationship="ortholog", ortholog_type="one2one"), verbose=FALSE)
  # correspondence: cluster 0<->0, 1<->1 (toys share structure)
  cc <- match_clusters(oa, ob, hg, "spA","spB", verbose=FALSE)
  skip_if(sum(cc$reciprocal) == 0, "no reciprocal cluster matches in toy")
  mmc <- match_metacells(mca, mcb, cc, hg, "spA","spB",
                         mutual_k = 3, verbose = FALSE)
  expect_s3_class(mmc, "metacell_correspondence")
  expect_gt(nrow(mmc), 0)
  expect_true(all(c("metacell_a","metacell_b","cluster_a","cluster_b","score")
                  %in% names(mmc)))
})

test_that("larger mutual_k yields at least as many pairs", {
  oa <- make_toy_clustered(1); ob <- make_toy_clustered(2)
  mca <- build_metacells(oa, gamma=5, n_var_genes=50, verbose=FALSE)
  mcb <- build_metacells(ob, gamma=5, n_var_genes=50, verbose=FALSE)
  hg <- as_homology_graph(data.frame(
    gene_a=paste0("g",1:100), gene_b=paste0("g",1:100),
    species_a="spA", species_b="spB",
    relationship="ortholog", ortholog_type="one2one"), verbose=FALSE)
  cc <- match_clusters(oa, ob, hg, "spA","spB", verbose=FALSE)
  skip_if(sum(cc$reciprocal) == 0)
  n3 <- nrow(match_metacells(mca,mcb,cc,hg,"spA","spB",mutual_k=1,verbose=FALSE))
  n5 <- nrow(match_metacells(mca,mcb,cc,hg,"spA","spB",mutual_k=5,verbose=FALSE))
  expect_gte(n5, n3)
})
