# data-raw/load_mouse_lung.R
# Load the Tabula Muris Senis droplet lung h5ad into a Seurat object, using raw
# counts from raw.X. Older AnnData convention: obs/var are compound datasets,
# sparse matrices lack encoding/shape attributes, raw counts in raw.X.
#
#   setwd("~/staging"); source("~/paralogswap/data-raw/load_mouse_lung.R")

library(hdf5r)
library(Matrix)
library(Seurat)

h5_path <- "tms_lung_droplet.h5ad"
out_rds <- "mouse_lung.rds"

f <- H5File$new(h5_path, mode = "r")

# ---- dimensions from the tables (attrs are absent in this convention) -------
genes <- f[["var"]][]$index
cells <- f[["obs"]][]$index
n_genes <- length(genes)   # 20138
n_cells <- length(cells)   # 24540

# ---- read the sparse raw.X --------------------------------------------------
rx      <- f[["raw.X"]]
data    <- rx[["data"]][]
indices <- rx[["indices"]][]   # 0-based
indptr  <- rx[["indptr"]][]    # 0-based

# Infer orientation from indptr length. AnnData is cells x genes.
#   CSC over cells x genes: indptr length = n_genes + 1
#   CSR over cells x genes: indptr length = n_cells + 1
len_ip <- length(indptr)
if (len_ip == n_genes + 1L) {
  # CSC: columns are genes, indices are cells. Build cells x genes, transpose.
  mat_cxg <- sparseMatrix(i = indices + 1L, p = indptr, x = data,
                          dims = c(n_cells, n_genes), index1 = TRUE)
  counts <- t(mat_cxg)
} else if (len_ip == n_cells + 1L) {
  # CSR over cells x genes = CSC over genes x cells: columns are cells,
  # indices are genes. Build genes x cells directly.
  counts <- sparseMatrix(i = indices + 1L, p = indptr, x = data,
                         dims = c(n_genes, n_cells), index1 = TRUE)
} else {
  stop("indptr length ", len_ip, " matches neither n_genes+1 (", n_genes+1,
       ") nor n_cells+1 (", n_cells+1, "). Inspect raw.X manually.")
}

rownames(counts) <- genes
colnames(counts) <- cells
cat("counts dims (genes x cells):", nrow(counts), "x", ncol(counts), "\n")
stopifnot(nrow(counts) == n_genes, ncol(counts) == n_cells)

# ---- read compound obs table into metadata ---------------------------------
obs_df <- f[["obs"]][]   # data.frame with all fields
want <- c("tissue", "subtissue", "cell_ontology_class", "free_annotation",
          "method", "mouse.id", "age", "sex", "n_genes", "n_counts")
want <- want[want %in% colnames(obs_df)]
meta <- obs_df[, want, drop = FALSE]
rownames(meta) <- cells

# hdf5r reads string fields as character already; categorical codes are not used
# in this convention (fields are stored as plain strings/numbers). Spot-check:
f$close_all()

# ---- build the Seurat object ------------------------------------------------
mouse_lung <- CreateSeuratObject(
  counts    = counts,
  meta.data = meta,
  project   = "tabula_muris_senis_lung"
)
cat("\nSeurat object:\n"); print(mouse_lung)

# ---- confirm tissue + inspect annotations -----------------------------------
cat("\ntissue:\n"); print(table(mouse_lung$tissue))
cat("\ncell_ontology_class (top 10):\n")
print(head(sort(table(mouse_lung$cell_ontology_class), decreasing = TRUE), 10))

# ---- GO / NO-GO: RAMP counts (mouse symbols are Title-case) -----------------
ramp <- c("Ramp1", "Ramp2", "Ramp3")
ramp <- ramp[ramp %in% rownames(mouse_lung)]
cat("\nRAMP genes in object:", paste(ramp, collapse=", "), "\n")
cm <- GetAssayData(mouse_lung, layer = "counts")[ramp, , drop = FALSE]
cat("total raw counts per RAMP gene:\n"); print(Matrix::rowSums(cm))
cat("cells expressing each (nonzero):\n"); print(Matrix::rowSums(cm > 0))

# ---- save -------------------------------------------------------------------
saveRDS(mouse_lung, out_rds)
cat("\nSaved", out_rds, "\n")
cat("Stage:  aws s3 cp", out_rds, "s3://paralogswap-data/results/\n")
