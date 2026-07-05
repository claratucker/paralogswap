# data-raw/load_lemur_lung.R
# Load the Tabula Microcebus lung h5ad into a Seurat object, using raw_counts.
# Robust hdf5r path (no converter): read the CSC sparse matrix by hand, build
# genes x cells, wire raw counts into Seurat, decode categorical metadata.
# Verified: RAMP1/2/3 present with counts; endothelial compartment present.
#
# Run from the directory containing the h5ad (e.g. ~/staging):
#   setwd("~/staging"); source("~/paralogswap/data-raw/load_lemur_lung.R")

library(hdf5r)
library(Matrix)
library(Seurat)

h5_path  <- "Lung_FIRM_hvg.h5ad"
out_rds  <- "lemur_lung.rds"

f <- H5File$new(h5_path, mode = "r")

# ---- read the CSC sparse raw_counts ----------------------------------------
# AnnData stores cells x genes. This layer is csc_matrix, shape (n_cells, n_genes).
# The (data, indices, indptr) triplet is the column-compressed form of a
# cells x genes matrix; we build that, then transpose to genes x cells for Seurat.
rc      <- f[["layers/raw_counts"]]
data    <- rc[["data"]][]
indices <- rc[["indices"]][]     # 0-based cell (row) indices
indptr  <- rc[["indptr"]][]      # 0-based column (gene) pointers, length n_genes + 1
shape   <- h5attr(rc, "shape")   # c(n_cells, n_genes)
n_cells <- shape[1]
n_genes <- shape[2]

mat_cells_x_genes <- sparseMatrix(
  i     = indices + 1L,          # 0-based -> 1-based
  p     = indptr,                # column pointers (genes)
  x     = data,
  dims  = c(n_cells, n_genes),
  index1 = TRUE
)
counts <- t(mat_cells_x_genes)   # -> genes x cells

# ---- names ------------------------------------------------------------------
genes <- f[["var/_index"]][]
cells <- f[["obs/_index"]][]
stopifnot(length(genes) == n_genes, length(cells) == n_cells)
rownames(counts) <- genes
colnames(counts) <- cells

cat("counts dims (genes x cells):", nrow(counts), "x", ncol(counts), "\n")
stopifnot(nrow(counts) == n_genes, ncol(counts) == n_cells)  # transpose tripwire

# ---- read + DECODE obs metadata --------------------------------------------
# Handles both AnnData categorical conventions:
#   (a) newer: each column is a group with its own codes + categories
#   (b) older: integer codes in the column, labels in obs/__categories/<name>
obs_grp  <- f[["obs"]]
cats_grp <- if ("__categories" %in% names(obs_grp)) obs_grp[["__categories"]] else NULL

read_obs_col <- function(name) {
  node <- obs_grp[[name]]
  if (inherits(node, "H5Group")) {                 # convention (a)
    codes <- node[["codes"]][]
    cats  <- node[["categories"]][]
    return(cats[codes + 1L])
  }
  vals <- node[]
  if (!is.null(cats_grp) && name %in% names(cats_grp)) {  # convention (b)
    cats <- cats_grp[[name]][]
    return(cats[vals + 1L])                        # 0-based codes -> 1-based lookup
  }
  vals                                             # genuinely numeric (nCount_RNA, etc.)
}

want <- c("tissue", "subtissue", "cell_ontology_class_v1",
          "free_annotation_v1", "compartment_v1", "individual",
          "method", "sex", "age", "nCount_RNA", "nFeature_RNA")
want <- want[want %in% names(obs_grp)]

meta <- as.data.frame(setNames(lapply(want, read_obs_col), want),
                      stringsAsFactors = FALSE)
rownames(meta) <- cells

f$close_all()

# ---- build the Seurat object ------------------------------------------------
lemur_lung <- CreateSeuratObject(
  counts    = counts,
  meta.data = meta,
  project   = "tabula_microcebus_lung"
)
cat("\nSeurat object:\n"); print(lemur_lung)

# ---- verification -----------------------------------------------------------
cat("\ntissue:\n");      print(table(lemur_lung$tissue))
cat("\ncompartment:\n"); print(table(lemur_lung$compartment_v1))

ramp <- c("RAMP1", "RAMP2", "RAMP3")
ramp <- ramp[ramp %in% rownames(lemur_lung)]
cat("\nRAMP genes in object:", paste(ramp, collapse = ", "), "\n")
cm <- GetAssayData(lemur_lung, layer = "counts")[ramp, , drop = FALSE]
cat("total raw counts per RAMP gene:\n"); print(Matrix::rowSums(cm))
cat("cells expressing each (nonzero):\n"); print(Matrix::rowSums(cm > 0))

# ---- save (gitignored *.rds; also push to S3 manually) ----------------------
saveRDS(lemur_lung, out_rds)
cat("\nSaved", out_rds, "\n")
cat("Stage to S3:  aws s3 cp", out_rds,
    "s3://paralogswap-data/results/\n")
