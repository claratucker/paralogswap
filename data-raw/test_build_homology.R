# data-raw/test_build_homology.R
# Interactive check that build_homology_graph works on the RAMP genes.
# Confirms: (1) both-directions union recovers RAMP orthologs, (2) the Ensembl
# RAMP2 mis-assignment reproduces without correction, (3) the eggNOG correction
# repairs it. Run with devtools::load_all(".") first, from the repo root.
#
#   library(devtools); load_all("."); source("data-raw/test_build_homology.R")

# ---- Human RAMP Ensembl IDs (from the spike, CHECK 1) ------------------------
ramp_ids <- c(
  RAMP1 = "ENSG00000132329",
  RAMP2 = "ENSG00000131477",
  RAMP3 = "ENSG00000122679"
)

# ============================================================================
# 1. WITHOUT correction — expect the RAMP2 gap the spike found
# ============================================================================
cat("\n########## build WITHOUT correction (hsapiens -> mmusculus) ##########\n")
hg_raw <- build_homology_graph(
  species_a = "hsapiens",
  species_b = "mmusculus",
  genes     = ramp_ids
)

cat("\n--- summary(hg_raw) ---\n")
summary(hg_raw)

cat("\n--- all edges ---\n")
print(as.data.frame(hg_raw))

# Focus: did human RAMP2 (ENSG00000131477) get a mouse ortholog?
cat("\n--- RAMP2 ortholog edges (expect MISSING or wrong, per spike) ---\n")
ramp2_edges <- hg_raw[hg_raw$gene_a == "ENSG00000131477" &
                        hg_raw$relationship == "ortholog", ]
if (nrow(ramp2_edges) == 0) {
  cat("  RAMP2 has NO cross-species ortholog edge — reproduces the spike.\n")
} else {
  print(as.data.frame(ramp2_edges))
}

# ============================================================================
# 2. WITH the eggNOG correction — Ensembl IDs to match id_type = "ensembl"
# ============================================================================
# Mouse Ramp2 = ENSMUSG00000001240 (from the spike). The correction must be in
# the SAME id space (ensembl) as the rest of the graph, or it won't join.
ramp2_fix <- data.frame(
  gene_a        = "ENSG00000131477",   # human RAMP2
  gene_b        = "ENSMUSG00000001240", # mouse Ramp2
  species_a     = "hsapiens",
  species_b     = "mmusculus",
  relationship  = "ortholog",
  ortholog_type = "one2one",
  source        = "eggnog",
  stringsAsFactors = FALSE
)

cat("\n########## build WITH eggNOG correction ##########\n")
hg_fixed <- build_homology_graph(
  species_a   = "hsapiens",
  species_b   = "mmusculus",
  genes       = ramp_ids,
  corrections = ramp2_fix
)

cat("\n--- summary(hg_fixed) ---  (source breakdown should show eggnog=1)\n")
summary(hg_fixed)

cat("\n--- RAMP2 ortholog edge after correction (expect ONE, source=eggnog) ---\n")
ramp2_fixed <- hg_fixed[hg_fixed$gene_a == "ENSG00000131477" &
                          hg_fixed$relationship == "ortholog", ]
print(as.data.frame(ramp2_fixed))

# ============================================================================
# 3. Eyeball checks — read these, don't just run past them
# ============================================================================
cat("\n########## EYEBALL CHECKS ##########\n")

# (a) Both-directions union: RAMP1 and RAMP3 orthologs present, ONE row each
#     (not duplicated by the A->B and B->A pulls).
for (g in c("ENSG00000132329", "ENSG00000122679")) {  # RAMP1, RAMP3
  e <- hg_fixed[hg_fixed$gene_a == g & hg_fixed$relationship == "ortholog", ]
  cat(sprintf("  %s: %d ortholog edge(s) [expect 1]\n", g, nrow(e)))
}

# (b) Corrections override, not append: exactly ONE RAMP2 ortholog edge,
#     and its source is eggnog (Ensembl's version, if any, was removed).
cat(sprintf("  RAMP2: %d ortholog edge(s) [expect 1], source(s): %s\n",
            nrow(ramp2_fixed),
            paste(unique(ramp2_fixed$source), collapse = ", ")))

# (c) Paralog edges present for BOTH species (symmetric-X): RAMP family members
#     should appear as within-species paralogs on each side.
cat("  paralog edges by species:\n")
para <- hg_fixed[hg_fixed$relationship == "paralog", ]
print(table(para$species_a))

# (d) Invariant sanity: every paralog within-species, every ortholog cross.
ok_para  <- all(para$species_a == para$species_b)
orth     <- hg_fixed[hg_fixed$relationship == "ortholog", ]
ok_orth  <- all(orth$species_a != orth$species_b)
cat(sprintf("  invariant (paralog within / ortholog cross): %s / %s\n",
            ok_para, ok_orth))

cat("\n########## DONE ##########\n")
cat("Read: is RAMP2 missing pre-correction and present (eggnog) post? Are\n")
cat("RAMP1/RAMP3 single-rowed (union dedup working)? Paralogs on both sides?\n")
