# Raw edge rows in the shape pull_orthologs / pull_paralogs produce, mirroring
# the verified RAMP result. These feed assemble_homology_edges directly (no net).

ortho_rows <- function() {
  # RAMP1: forward (human->mouse). RAMP3: reverse only (mouse->human), higher
  # perc_id, to force the direction-flip path. RAMP1 also appears reverse with
  # lower perc_id to test dedup preferring the better record.
  data.frame(
    gene_a = c("ENSG_RAMP1", "ENSMUS_RAMP1",  "ENSMUS_RAMP3"),
    gene_b = c("ENSMUS_RAMP1", "ENSG_RAMP1",  "ENSG_RAMP3"),
    species_a = c("hsapiens", "mmusculus", "mmusculus"),
    species_b = c("mmusculus", "hsapiens", "hsapiens"),
    relationship = "ortholog",
    ortholog_type = "one2one",
    duplication_level = NA_character_,
    perc_id = c(80, 40, 84),          # RAMP1 fwd(80)>rev(40); RAMP3 rev only(84)
    confidence = c(1, 1, 1),
    source = "ensembl",
    stringsAsFactors = FALSE
  )
}

para_rows <- function() {
  data.frame(
    gene_a = c("ENSG_RAMP1", "ENSMUS_RAMP1"),
    gene_b = c("ENSG_RAMP3", "ENSMUS_RAMP3"),
    species_a = c("hsapiens", "mmusculus"),
    species_b = c("hsapiens", "mmusculus"),
    relationship = "paralog",
    ortholog_type = NA_character_,
    duplication_level = NA_character_,
    perc_id = c(30, 25),
    confidence = NA,
    source = "ensembl",
    stringsAsFactors = FALSE
  )
}

test_that("dedup collapses both-direction ortholog pairs to one row", {
  hg <- assemble_homology_edges(ortho_rows(), para_rows(),
                                "hsapiens", "mmusculus", verbose = FALSE)
  orth <- hg[hg$relationship == "ortholog", ]
  # RAMP1 (2 input rows) + RAMP3 (1 row) -> 2 unique ortholog edges
  expect_equal(nrow(orth), 2)
})

test_that("ortholog direction is normalized to the focal species", {
  hg <- assemble_homology_edges(ortho_rows(), para_rows(),
                                "hsapiens", "mmusculus", verbose = FALSE)
  orth <- hg[hg$relationship == "ortholog", ]
  # This is the RAMP3-direction-bug regression test: EVERY ortholog edge must
  # have the focal species on side A, even though RAMP3 came in reversed.
  expect_true(all(orth$species_a == "hsapiens"))
  expect_true(all(orth$species_b == "mmusculus"))
})

test_that("dedup keeps the higher-confidence / higher-perc_id record", {
  hg <- assemble_homology_edges(ortho_rows(), para_rows(),
                                "hsapiens", "mmusculus", verbose = FALSE)
  ramp1 <- hg[hg$relationship == "ortholog" & hg$gene_a == "ENSG_RAMP1", ]
  expect_equal(nrow(ramp1), 1)
  expect_equal(ramp1$perc_id, 80)   # the forward record, not the 40 reverse
})

test_that("paralog edges are retained for both species", {
  hg <- assemble_homology_edges(ortho_rows(), para_rows(),
                                "hsapiens", "mmusculus", verbose = FALSE)
  para <- hg[hg$relationship == "paralog", ]
  expect_setequal(unique(para$species_a), c("hsapiens", "mmusculus"))
})

test_that("corrections override matching ensembl edges", {
  # Correction replaces the RAMP1 ortholog with an eggnog-sourced edge.
  corr <- data.frame(
    gene_a = "ENSG_RAMP1", gene_b = "ENSMUS_RAMP1",
    species_a = "hsapiens", species_b = "mmusculus",
    relationship = "ortholog", ortholog_type = "one2one",
    source = "eggnog", stringsAsFactors = FALSE
  )
  hg <- assemble_homology_edges(ortho_rows(), para_rows(),
                                "hsapiens", "mmusculus",
                                corrections = corr, verbose = FALSE)
  ramp1 <- hg[hg$relationship == "ortholog" & hg$gene_a == "ENSG_RAMP1", ]
  expect_equal(nrow(ramp1), 1)          # not duplicated
  expect_equal(as.character(ramp1$source), "eggnog")  # ensembl one removed
})

test_that("corrections add a genuinely new edge", {
  # RAMP2 absent from ensembl rows; correction supplies it (the real scenario).
  corr <- data.frame(
    gene_a = "ENSG_RAMP2", gene_b = "ENSMUS_RAMP2",
    species_a = "hsapiens", species_b = "mmusculus",
    relationship = "ortholog", ortholog_type = "one2one",
    source = "eggnog", stringsAsFactors = FALSE
  )
  hg <- assemble_homology_edges(ortho_rows(), para_rows(),
                                "hsapiens", "mmusculus",
                                corrections = corr, verbose = FALSE)
  expect_true(any(hg$gene_a == "ENSG_RAMP2" & hg$source == "eggnog"))
})

test_that("min_perc_id filters low-identity edges but keeps NA perc_id", {
  hg <- assemble_homology_edges(ortho_rows(), para_rows(),
                                "hsapiens", "mmusculus",
                                min_perc_id = 50, verbose = FALSE)
  # RAMP1 fwd (80) and RAMP3 (84) survive; paralogs (30,25) drop.
  expect_true(all(hg$perc_id >= 50 | is.na(hg$perc_id)))
})

test_that("min_confidence='high' drops low-confidence orthologs, keeps paralogs", {
  rows <- ortho_rows(); rows$confidence <- c(1, 1, 0)  # RAMP3 low-conf
  hg <- assemble_homology_edges(rows, para_rows(),
                                "hsapiens", "mmusculus",
                                min_confidence = "high", verbose = FALSE)
  orth <- hg[hg$relationship == "ortholog", ]
  expect_false(any(orth$gene_a == "ENSG_RAMP3"))  # dropped
  expect_true(any(hg$relationship == "paralog"))  # paralogs (NA conf) kept
})

# ---- corrections: add / remove -------------------------------------------

test_that("action='remove' retracts an edge instead of adding it", {
  o <- ortho_rows(); p <- para_rows()
  base <- assemble_homology_edges(o, p, species_a = "hsapiens",
                                  species_b = "mmusculus", verbose = FALSE)
  n0 <- nrow(base)

  rm_one <- data.frame(
    gene_a = "ENSG_RAMP1", gene_b = "ENSG_RAMP3",
    species_a = "hsapiens", species_b = "hsapiens",
    relationship = "paralog", ortholog_type = NA_character_,
    source = "eggnog", action = "remove", stringsAsFactors = FALSE)

  out <- assemble_homology_edges(o, p, species_a = "hsapiens",
                                 species_b = "mmusculus",
                                 corrections = rm_one, verbose = FALSE)

  expect_equal(nrow(out), n0 - 1L)
  d <- as.data.frame(out)
  expect_false(any(d$gene_a == "ENSG_RAMP1" & d$gene_b == "ENSG_RAMP3" &
                   d$relationship == "paralog"))
  ex <- attr(out, "excluded")
  expect_equal(nrow(ex), 1L)
  expect_identical(as.character(ex$source[1]), "ensembl")
})

test_that("removal matches regardless of edge orientation", {
  o <- ortho_rows(); p <- para_rows()
  rm_rev <- data.frame(
    gene_a = "ENSG_RAMP3", gene_b = "ENSG_RAMP1",
    species_a = "hsapiens", species_b = "hsapiens",
    relationship = "paralog", ortholog_type = NA_character_,
    source = "eggnog", action = "remove", stringsAsFactors = FALSE)

  out <- assemble_homology_edges(o, p, species_a = "hsapiens",
                                 species_b = "mmusculus",
                                 corrections = rm_rev, verbose = FALSE)
  expect_equal(nrow(attr(out, "excluded")), 1L)
})

test_that("removals are applied before additions for the same pair", {
  o <- ortho_rows(); p <- para_rows()
  corr <- rbind(
    data.frame(gene_a = "ENSG_RAMP1", gene_b = "ENSG_RAMP3",
               species_a = "hsapiens", species_b = "hsapiens",
               relationship = "paralog", ortholog_type = NA_character_,
               source = "eggnog", action = "remove", stringsAsFactors = FALSE),
    data.frame(gene_a = "ENSG_RAMP1", gene_b = "ENSG_RAMP3",
               species_a = "hsapiens", species_b = "hsapiens",
               relationship = "paralog", ortholog_type = NA_character_,
               source = "genetree", action = "add", stringsAsFactors = FALSE))

  out <- as.data.frame(assemble_homology_edges(
    o, p, species_a = "hsapiens", species_b = "mmusculus",
    corrections = corr, verbose = FALSE))

  row <- out[out$gene_a == "ENSG_RAMP1" & out$gene_b == "ENSG_RAMP3" &
             out$relationship == "paralog", ]
  expect_equal(nrow(row), 1L)
  expect_identical(as.character(row$source[1]), "genetree")
})

test_that("corrections without an action column still add (back-compat)", {
  o <- ortho_rows(); p <- para_rows()
  add_only <- data.frame(
    gene_a = "ENSG_RAMP2", gene_b = "ENSMUS_RAMP2",
    species_a = "hsapiens", species_b = "mmusculus",
    relationship = "ortholog", ortholog_type = "one2one",
    source = "eggnog", stringsAsFactors = FALSE)

  out <- assemble_homology_edges(o, p, species_a = "hsapiens",
                                 species_b = "mmusculus",
                                 corrections = add_only, verbose = FALSE)
  expect_true(any(as.data.frame(out)$source == "eggnog"))
  expect_null(attr(out, "excluded"))
})

test_that("an unrecognized action is rejected", {
  o <- ortho_rows(); p <- para_rows()
  bad <- data.frame(gene_a = "X", gene_b = "Y",
                    species_a = "hsapiens", species_b = "hsapiens",
                    relationship = "paralog", ortholog_type = NA_character_,
                    source = "manual", action = "delete", stringsAsFactors = FALSE)
  expect_error(
    assemble_homology_edges(o, p, species_a = "hsapiens", species_b = "mmusculus",
                            corrections = bad, verbose = FALSE),
    "add")
})

test_that("as_homology_graph refuses a corrections table containing removals", {
  rm_row <- data.frame(gene_a = "X", gene_b = "Y",
                       relationship = "paralog", ortholog_type = NA_character_,
                       action = "remove", stringsAsFactors = FALSE)
  expect_error(as_homology_graph(rm_row, verbose = FALSE), "removals")
})

# ---- b_ids sanitation (the empty-string filter bug) -----------------------

test_that(".clean_ids drops empty and NA ids that would unrestrict a pull", {
  expect_equal(paralogswap:::.clean_ids(c("Ramp1", "", "Ramp3", NA)),
               c("Ramp1", "Ramp3"))
  expect_null(paralogswap:::.clean_ids(c("", NA, "  ")))
  expect_null(paralogswap:::.clean_ids(NULL))
})
