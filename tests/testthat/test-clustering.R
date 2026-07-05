# Small synthetic Seurat objects for fast, deterministic tests. Two "species"
# with structured expression so clustering finds >1 cluster.
make_toy <- function(n_cells = 120, n_genes = 200, seed = 1) {
  set.seed(seed)
  # two cell groups with different mean expression -> separable
  half <- n_cells / 2
  base <- matrix(rpois(n_genes * n_cells, 2), n_genes, n_cells)
  base[1:20, 1:half]            <- base[1:20, 1:half] + 15      # group A markers
  base[21:40, (half+1):n_cells] <- base[21:40, (half+1):n_cells] + 15  # group B
  rownames(base) <- paste0("g", seq_len(n_genes))
  colnames(base) <- paste0("c", seq_len(n_cells))
  Seurat::CreateSeuratObject(counts = Matrix::Matrix(base, sparse = TRUE))
}

test_that("cluster_species accepts a named list and returns species_clusters", {
  sc <- cluster_species(list(a = make_toy(), b = make_toy(seed = 2)),
                        n_pcs = 10, n_var_genes = 50, verbose = FALSE)
  expect_s3_class(sc, "species_clusters")
  expect_named(sc, c("a", "b"))
  expect_true(all(vapply(sc, inherits, logical(1), "Seurat")))
})

test_that("each returned object carries PCA and cluster assignments", {
  sc <- cluster_species(list(a = make_toy()),
                        n_pcs = 10, n_var_genes = 50, verbose = FALSE)
  o <- sc$a
  expect_true("pca" %in% Seurat::Reductions(o))
  expect_true("seurat_clusters" %in% colnames(o[[]]))
  expect_gte(length(unique(Seurat::Idents(o))), 1)
})

test_that("params are recorded as an attribute", {
  sc <- cluster_species(list(a = make_toy()), resolution = 0.5,
                        n_pcs = 10, n_var_genes = 50, verbose = FALSE)
  p <- attr(sc, "params")
  expect_equal(p$resolution, 0.5)
  expect_equal(p$normalization, "lognorm")
})

test_that("vector resolution produces per-resolution cluster columns", {
  sc <- cluster_species(list(a = make_toy()), resolution = c(0.5, 1.0),
                        n_pcs = 10, n_var_genes = 50, verbose = FALSE)
  o <- sc$a
  expect_true(all(c("clusters_res0.5", "clusters_res1") %in% colnames(o[[]])))
})

test_that("single merged object is split by species_col", {
  a <- make_toy(); b <- make_toy(seed = 2)
  # distinct cell names so merge doesn't rename
  a <- RenameCells(a, add.cell.id = "A")
  b <- RenameCells(b, add.cell.id = "B")
  a$species <- "sp_a"; b$species <- "sp_b"
  merged <- merge(a, b)
  # Seurat 5 splits counts into per-object layers after merge; rejoin if available
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    merged <- SeuratObject::JoinLayers(merged)
  }
  sc <- cluster_species(merged, species_col = "species",
                        n_pcs = 10, n_var_genes = 50, verbose = FALSE)
  expect_named(sc, c("sp_a", "sp_b"), ignore.order = TRUE)
})

test_that("unnamed list is rejected", {
  expect_error(
    cluster_species(list(make_toy(), make_toy()), verbose = FALSE),
    "must be NAMED"
  )
})

test_that("missing species_col on merged object errors informatively", {
  expect_error(
    cluster_species(make_toy(), species_col = "nonexistent", verbose = FALSE),
    "not found in object metadata"
  )
})
