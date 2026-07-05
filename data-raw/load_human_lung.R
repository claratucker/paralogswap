# data-raw/load_human_lung.R
# Load Tabula Sapiens v2 lung h5ad -> Seurat, using raw_counts.
# Subsets to 10X (drops Smart-seq2), subsamples to ~40k for cross-species
# balance with lemur (34k) / mouse (24k). Symbols as gene names, ensembl_id
# carried as feature-level alias. CSR sparse format.
#
#   setwd("~/staging"); source("~/paralogswap/data-raw/load_human_lung.R")

library(hdf5r)
library(Matrix)
library(Seurat)

h5_path       <- "ts_lung.h5ad"
out_rds       <- "human_lung.rds"
subsample_to  <- 40000
set.seed(1)

f <- H5File$new(h5_path, mode = "r")

# ---- decode a group-categorical or plain obs column -------------------------
read_obs <- function(name) {
  node <- f[[paste0("obs/", name)]]
  if (inherits(node, "H5Group")) {
    cats <- node[["categories"]][]; codes <- node[["codes"]][]
    cats[codes + 1L]
  } else node[]
}

# ---- cell selection: 10X only, then subsample -------------------------------
cells_all <- f[["obs/_index"]][]
method    <- read_obs("method")
is_10x    <- which(method == "10X")
cat("10X cells:", length(is_10x), "/ total", length(cells_all), "\n")

if (length(is_10x) > subsample_to) {
  keep_idx <- sort(sample(is_10x, subsample_to))
  cat("subsampled to", subsample_to, "for cross-species balance\n")
} else {
  keep_idx <- is_10x
  cat("kept all", length(is_10x), "10X cells (<= target)\n")
}
n_keep <- length(keep_idx)

# ---- read CSR raw_counts, build genes x cells, subset to kept cells ---------
# CSR over cells x genes: indptr length = n_cells+1, each row (cell) is a block.
# Reading whole then subsetting columns is simplest at this size (66k x 62k).
rc      <- f[["layers/raw_counts"]]
shape   <- h5attr(rc, "shape")          # c(n_cells, n_genes)
n_cells <- shape[1]; n_genes <- shape[2]
data    <- rc[["data"]][]
indices <- rc[["indices"]][]            # 0-based gene indices
indptr  <- rc[["indptr"]][]             # 0-based, length n_cells+1

# CSR(cells x genes): build as a cells x genes dgRMatrix-equivalent, i.e. use
# sparseMatrix with j = gene indices, p = cell pointers -> then it's cells x genes.
# sparseMatrix wants column-pointers, so build transposed directly: genes x cells
# is CSC with p over cells. CSR row-ptr over cells == CSC col-ptr over cells, so:
counts <- sparseMatrix(
  i     = indices + 1L,   # gene indices become ROW indices of genes x cells
  p     = indptr,         # cell pointers become COLUMN pointers
  x     = data,
  dims  = c(n_genes, n_cells),
  index1 = TRUE
)   # -> genes x cells, no transpose needed

genes <- f[["var/index"]][]
ens   <- f[["var/ensembl_id"]][]
rownames(counts) <- genes
colnames(counts) <- cells_all
stopifnot(nrow(counts) == n_genes, ncol(counts) == n_cells)

# subset to kept cells
counts <- counts[, keep_idx, drop = FALSE]
cat("counts after subset (genes x cells):", nrow(counts), "x", ncol(counts), "\n")

# ---- metadata for kept cells ------------------------------------------------
want <- c("tissue","compartment","cell_ontology_class","free_annotation",
          "broad_cell_class","method","assay","donor","age","sex","anatomical_position")
want <- want[want %in% names(f[["obs"]])]
meta <- as.data.frame(setNames(lapply(want, function(w) read_obs(w)[keep_idx]), want),
                      stringsAsFactors = FALSE)
rownames(meta) <- cells_all[keep_idx]

f$close_all()

# ---- build Seurat object ----------------------------------------------------
human_lung <- CreateSeuratObject(
  counts    = counts,
  meta.data = meta,
  project   = "tabula_sapiens_lung"
)
# carry ensembl_id as feature metadata alias
human_lung[["RNA"]] <- AddMetaData(human_lung[["RNA"]],
                                   metadata = setNames(ens, genes)[rownames(human_lung)],
                                   col.name = "ensembl_id")
cat("\nSeurat object:\n"); print(human_lung)

# ---- checks -----------------------------------------------------------------
cat("\ntissue:\n"); print(table(human_lung$tissue))
cat("\ncompartment:\n"); print(table(human_lung$compartment))

ramp <- c("RAMP1","RAMP2","RAMP3")
ramp <- ramp[ramp %in% rownames(human_lung)]
cat("\nRAMP genes:", paste(ramp, collapse=", "), "\n")
cm <- GetAssayData(human_lung, layer = "counts")[ramp, , drop = FALSE]
cat("total raw counts:\n"); print(Matrix::rowSums(cm))
cat("cells expressing (nonzero):\n"); print(Matrix::rowSums(cm > 0))

saveRDS(human_lung, out_rds)
cat("\nSaved", out_rds, "\n")
cat("Stage:  aws s3 cp", out_rds, "s3://paralogswap-data/results/\n")
