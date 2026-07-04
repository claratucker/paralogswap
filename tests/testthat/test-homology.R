# Helper: a minimal valid edge table in the new schema.
valid_edges <- function() {
  data.frame(
    gene_a        = c("RAMP1", "RAMP2", "RAMP1"),
    gene_b        = c("Ramp1", "Ramp2", "RAMP3"),
    species_a     = c("hsapiens", "hsapiens", "hsapiens"),
    species_b     = c("mmusculus", "mmusculus", "hsapiens"),
    relationship  = c("ortholog", "ortholog", "paralog"),
    ortholog_type = c("one2one", "one2one", NA),
    stringsAsFactors = FALSE
  )
}

test_that("as_homology_graph builds a valid object from minimal input", {
  hg <- as_homology_graph(valid_edges(), verbose = FALSE)
  expect_s3_class(hg, "homology_graph")
  expect_s3_class(hg, "data.frame")
  expect_equal(nrow(hg), 3)
  expect_true(all(c("species_a", "species_b", "ortholog_type",
                    "ambiguous", "source") %in% colnames(hg)))
})

test_that("relationship is a two-level factor (ortholog, paralog)", {
  hg <- as_homology_graph(valid_edges(), verbose = FALSE)
  expect_s3_class(hg$relationship, "factor")
  expect_identical(levels(hg$relationship), c("ortholog", "paralog"))
})

test_that("missing required columns raise an informative error", {
  bad <- data.frame(gene_a = "x", gene_b = "y")  # no species/relationship/type
  expect_error(as_homology_graph(bad, verbose = FALSE),
               "missing required column")
})

test_that("invalid relationship values are rejected", {
  bad <- valid_edges(); bad$relationship[1] <- "homolog"
  expect_error(as_homology_graph(bad, verbose = FALSE),
               "ortholog.*paralog")
})

test_that("invalid ortholog_type values are rejected", {
  bad <- valid_edges(); bad$ortholog_type[1] <- "one2three"
  expect_error(as_homology_graph(bad, verbose = FALSE),
               "one2one.*one2many.*many2many")
})

test_that("ambiguous is derived from ortholog_type, not trusted from input", {
  e <- valid_edges()
  e$ortholog_type[1] <- "many2many"
  e$ambiguous <- FALSE          # deliberately wrong; must be overwritten
  hg <- as_homology_graph(e, verbose = FALSE)
  expect_true(hg$ambiguous[hg$ortholog_type == "many2many"][1])
  expect_false(any(hg$ambiguous[hg$relationship == "paralog"]))
})

test_that("source defaults to 'manual' when absent", {
  hg <- as_homology_graph(valid_edges(), verbose = FALSE)
  expect_true(all(hg$source == "manual"))
})

test_that("supplied source is preserved", {
  e <- valid_edges(); e$source <- "eggnog"
  hg <- as_homology_graph(e, verbose = FALSE)
  expect_true(all(hg$source == "eggnog"))
})

test_that("within/cross-species invariant is enforced", {
  # paralog edge marked cross-species -> error
  bad1 <- valid_edges(); bad1$species_b[3] <- "mmusculus"  # paralog now cross
  expect_error(as_homology_graph(bad1, verbose = FALSE),
               "paralog edges must have species_a == species_b")
  # ortholog edge marked within-species -> error
  bad2 <- valid_edges(); bad2$species_b[1] <- "hsapiens"   # ortholog now within
  expect_error(as_homology_graph(bad2, verbose = FALSE),
               "ortholog edges must have species_a != species_b")
})

test_that("gene IDs supplied as factors are coerced to character", {
  e <- valid_edges()
  e$gene_a <- factor(e$gene_a)
  hg <- as_homology_graph(e, verbose = FALSE)
  expect_type(hg$gene_a, "character")
})

test_that("perc_id outside [0,100] is rejected", {
  e <- valid_edges(); e$perc_id <- c(50, 150, NA)
  expect_error(as_homology_graph(e, verbose = FALSE), "perc_id")
})

test_that("as_homology_graph is idempotent", {
  hg  <- as_homology_graph(valid_edges(), verbose = FALSE)
  hg2 <- as_homology_graph(hg, verbose = FALSE)
  expect_identical(hg, hg2)
})

test_that("summary returns the finalized breakdown", {
  hg <- as_homology_graph(valid_edges(), verbose = FALSE)
  s  <- summary(hg)
  expect_equal(s$n_edges, 3)
  expect_equal(s$n_ortholog, 2)
  expect_equal(s$n_paralog, 1)
  expect_equal(s$n_one2one, 2)
})
