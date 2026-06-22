test_that("as_homology_graph builds a valid object from minimal input", {
  edges <- data.frame(
    gene_a = c("Actb", "Gapdh", "Myl2"),
    gene_b = c("ACTB", "GAPDH", "MYL7"),
    relationship = c("ortholog", "ortholog", "paralog")
  )
  hg <- as_homology_graph(edges, verbose = FALSE)

  expect_s3_class(hg, "homology_graph")
  expect_s3_class(hg, "data.frame")
  expect_equal(nrow(hg), 3)
  expect_true(all(c("duplication_level", "seq_similarity") %in% colnames(hg)))
  expect_true(all(is.na(hg$seq_similarity)))
})

test_that("relationship is a three-level factor including unclassified", {
  edges <- data.frame(
    gene_a = c("a1", "a2", "a3"),
    gene_b = c("b1", "b2", "b3"),
    relationship = c("ortholog", "paralog", "unclassified")
  )
  hg <- as_homology_graph(edges, verbose = FALSE)

  expect_s3_class(hg$relationship, "factor")
  expect_identical(
    levels(hg$relationship),
    c("ortholog", "paralog", "unclassified")
  )
})

test_that("missing required columns raise an informative error", {
  bad <- data.frame(gene_a = "x", gene_b = "y")
  expect_error(
    as_homology_graph(bad, verbose = FALSE),
    "missing required column"
  )
})

test_that("invalid relationship values are rejected", {
  bad <- data.frame(
    gene_a = "x", gene_b = "y", relationship = "homolog"
  )
  expect_error(
    as_homology_graph(bad, verbose = FALSE),
    "ortholog.*paralog.*unclassified"
  )
})

test_that("non-data.frame input is rejected", {
  expect_error(
    as_homology_graph(list(gene_a = "x"), verbose = FALSE),
    "must be a data.frame"
  )
})

test_that("gene IDs supplied as factors are coerced to character", {
  edges <- data.frame(
    gene_a = factor(c("a1", "a2")),
    gene_b = factor(c("b1", "b2")),
    relationship = c("ortholog", "paralog"),
    stringsAsFactors = TRUE
  )
  hg <- as_homology_graph(edges, verbose = FALSE)
  expect_type(hg$gene_a, "character")
  expect_type(hg$gene_b, "character")
})

test_that("rows with missing gene IDs are dropped", {
  edges <- data.frame(
    gene_a = c("a1", NA, "a3"),
    gene_b = c("b1", "b2", NA),
    relationship = c("ortholog", "ortholog", "paralog")
  )
  hg <- as_homology_graph(edges, verbose = FALSE)
  expect_equal(nrow(hg), 1)
})

test_that("seq_similarity outside [0, 1] is rejected", {
  bad <- data.frame(
    gene_a = "x", gene_b = "y", relationship = "ortholog",
    seq_similarity = 1.5
  )
  expect_error(
    as_homology_graph(bad, verbose = FALSE),
    "seq_similarity"
  )
})

test_that("as_homology_graph is idempotent", {
  edges <- data.frame(
    gene_a = "a1", gene_b = "b1", relationship = "ortholog"
  )
  hg <- as_homology_graph(edges, verbose = FALSE)
  hg2 <- as_homology_graph(hg, verbose = FALSE)
  expect_identical(hg, hg2)
})

test_that("summary returns counts and percentages", {
  edges <- data.frame(
    gene_a = c("a1", "a2", "a3", "a4"),
    gene_b = c("b1", "b2", "b3", "b4"),
    relationship = c("ortholog", "ortholog", "paralog", "unclassified")
  )
  hg <- as_homology_graph(edges, verbose = FALSE)
  s <- summary(hg)

  expect_equal(s$n_edges, 4)
  expect_equal(s$n_ortholog, 2)
  expect_equal(s$n_paralog, 1)
  expect_equal(s$n_unclassified, 1)
  expect_equal(s$pct_ortholog, 50)
  expect_equal(s$pct_paralog, 25)
  expect_equal(s$pct_unclassified, 25)
})

test_that("summary reports multi-mapping genes", {
  edges <- data.frame(
    gene_a = c("a1", "a1", "a2"),
    gene_b = c("b1", "b2", "b3"),
    relationship = c("ortholog", "paralog", "ortholog")
  )
  hg <- as_homology_graph(edges, verbose = FALSE)
  s <- summary(hg)
  expect_equal(s$multi_a, 1)
  expect_equal(s$multi_b, 0)
})

test_that("empty graph summarizes without error and returns NA percentages", {
  empty <- data.frame(
    gene_a = character(0),
    gene_b = character(0),
    relationship = character(0)
  )
  hg <- as_homology_graph(empty, verbose = FALSE)
  expect_equal(nrow(hg), 0)
  s <- summary(hg)
  expect_equal(s$n_edges, 0)
  expect_true(is.na(s$pct_paralog))
})

test_that("unmapped_rate attribute is reported when present, NULL otherwise", {
  edges <- data.frame(
    gene_a = "a1", gene_b = "b1", relationship = "ortholog"
  )
  hg <- as_homology_graph(edges, verbose = FALSE)

  expect_null(summary(hg)$unmapped_rate)

  attr(hg, "unmapped_rate") <- 0.2
  expect_equal(summary(hg)$unmapped_rate, 0.2)
})
