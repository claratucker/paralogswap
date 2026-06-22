#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="seurat-env"

echo "Rebuilding environment..."

# If lockfile exists, use it (recommended reproducibility path)
if [ -f conda-lock.yml ]; then
    echo "Using conda-lock.yml"
    conda-lock install -n "$ENV_NAME" conda-lock.yml
else
    echo "No lockfile found, falling back to environment.yml"
    conda env create -n "$ENV_NAME" -f environment.yml
fi

echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "Running Seurat test..."

Rscript -e '
suppressPackageStartupMessages({
  library(Seurat)
})

mat <- matrix(rpois(200, 5), nrow=20)
rownames(mat) <- paste0("gene", 1:20)
colnames(mat) <- paste0("cell", 1:10)

obj <- CreateSeuratObject(counts = mat)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 5)

print(obj)
cat("Seurat pipeline OK\n")
'
