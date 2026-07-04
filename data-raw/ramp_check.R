# data-raw/ramp_check.R
# SPIKE — go/no-go for the three-species RAMP validation.
# Question: do RAMP1/2/3 resolve as clean homologs across human/lemur/mouse in
# Ensembl, and what ortholog_type does each cross-species pair get?
# No atlas data, no package code — just biomaRt. Run interactively; read the
# CHECK lines as you go.

library(biomaRt)

# ---- CHECK 0: does Ensembl even have mouse lemur? -----------------------------
# Mouse lemur is Microcebus murinus -> dataset prefix "mmurinus".
# If it's absent from the release biomaRt hits, the 3-species plan is dead here
# and we fall back to mouse-human with lemur as documented future work.
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

lemur_hit <- datasets[grepl("murinus|lemur|microcebus", datasets$dataset,
                            ignore.case = TRUE) |
                      grepl("murinus|lemur|microcebus", datasets$description,
                            ignore.case = TRUE), ]
cat("\n=== CHECK 0: lemur in Ensembl? ===\n")
print(lemur_hit)
stopifnot("Mouse lemur dataset not found in this Ensembl release" =
            nrow(lemur_hit) > 0)

# Pin the exact dataset names we'll use. Adjust if CHECK 0 shows a different id.
ds_human <- "hsapiens_gene_ensembl"
ds_lemur <- lemur_hit$dataset[1]          # expected: "mmurinus_gene_ensembl"
ds_mouse <- "mmusculus_gene_ensembl"
cat("Using lemur dataset:", ds_lemur, "\n")

# Record the release so the whole spike is reproducible.
cat("Ensembl release:", listEnsembl()$version[listEnsembl()$biomart == "genes"], "\n")

# ---- Connect to the three marts ----------------------------------------------
mart_human <- useEnsembl("genes", dataset = ds_human)
mart_lemur <- useEnsembl("genes", dataset = ds_lemur)
mart_mouse <- useEnsembl("genes", dataset = ds_mouse)

# ---- CHECK 1: do the RAMP genes exist in each species? ------------------------
# Look them up by symbol. Lemur is the risk: symbols may be absent or "LOC…".
ramp_symbols <- c("RAMP1", "RAMP2", "RAMP3")

get_ramp <- function(mart, symbols, label) {
  # Try external_gene_name; report what comes back.
  res <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name",
                   "chromosome_name", "start_position"),
    filters = "external_gene_name",
    values = symbols,
    mart = mart
  )
  cat("\n--- RAMP genes in", label, "(by symbol) ---\n")
  print(res)
  res
}

cat("\n=== CHECK 1: RAMP genes present per species? ===\n")
ramp_human <- get_ramp(mart_human, ramp_symbols, "human")
ramp_lemur <- get_ramp(mart_lemur, c(ramp_symbols, tolower(ramp_symbols)),
                       "lemur")
ramp_mouse <- get_ramp(mart_mouse, c("Ramp1", "Ramp2", "Ramp3"), "mouse")

# ---- CHECK 2: human -> lemur and human -> mouse orthology + type --------------
# This is the core question. For each human RAMP, pull its lemur and mouse
# homologs WITH ortholog_type and % identity. ortholog_type is the field the
# whole many-to-many / substitution design depends on.
#
# Homolog attributes live on the HUMAN mart, named by the target species prefix.
# Target prefixes: lemur = mmurinus, mouse = mmusculus.

pull_homologs <- function(human_mart, human_ids, target_prefix, target_label) {
  attrs <- c(
    "ensembl_gene_id",                                   # human gene
    "external_gene_name",
    paste0(target_prefix, "_homolog_ensembl_gene"),      # target gene id
    paste0(target_prefix, "_homolog_associated_gene_name"),
    paste0(target_prefix, "_homolog_orthology_type"),    # <-- the key field
    paste0(target_prefix, "_homolog_perc_id"),           # % identity
    paste0(target_prefix, "_homolog_orthology_confidence")
  )
  res <- getBM(
    attributes = attrs,
    filters = "ensembl_gene_id",
    values = human_ids,
    mart = human_mart
  )
  cat("\n--- human -> ", target_label, " homologs (with ortholog_type) ---\n",
      sep = "")
  print(res)
  res
}

cat("\n=== CHECK 2: human->lemur and human->mouse orthology types ===\n")
human_ramp_ids <- ramp_human$ensembl_gene_id
stopifnot("No human RAMP Ensembl IDs found" = length(human_ramp_ids) > 0)

homo_lemur <- pull_homologs(mart_human, human_ramp_ids, "mmurinus", "lemur")
homo_mouse <- pull_homologs(mart_human, human_ramp_ids, "mmusculus", "mouse")

# ---- The verdict you actually read -------------------------------------------
cat("\n\n================= GO / NO-GO SUMMARY =================\n")
cat("For the 3-species RAMP validation to be reachable, you want to see:\n")
cat("  - each human RAMP1/2/3 with a lemur homolog and a mouse homolog\n")
cat("  - ortholog_type that is interpretable (one2one / one2many / many2many)\n")
cat("  - lemur gene names that are real symbols, not only LOC… placeholders\n")
cat("  - reasonable perc_id (not near-zero)\n\n")
cat("If lemur RAMP homologs are missing or all LOC/ambiguous: the finding may\n")
cat("not be cleanly reproducible from Ensembl homology alone — that is itself a\n")
cat("real result (it's what Ezran flag), but it changes the validation plan.\n")
cat("=====================================================\n")

# Save the raw pulls so we can inspect / reuse without re-querying.
saveRDS(
  list(ramp_human = ramp_human, ramp_lemur = ramp_lemur, ramp_mouse = ramp_mouse,
       homo_lemur = homo_lemur, homo_mouse = homo_mouse),
  "data-raw/ramp_check_result.rds"
)
cat("\nSaved raw pulls to data-raw/ramp_check_result.rds\n")
