# Environment

Reproducible single-cell RNA-seq environment for paralogswap (Seurat v5 workflows),
conda-managed, no renv.

## Reproducibility model

Two layers plus one GitHub-pinned package:

- **`environment.yml`** — human-readable spec: the direct dependencies, no build
  hashes. Read this to understand what the environment contains.
- **`conda-spec.txt`** — exact lock: every conda package pinned by URL and hash
  (`conda list --explicit`), linux-64. This is what reproduces the environment
  precisely.
- **SuperCell** — distributed only via GitHub (not on conda or CRAN), so it cannot
  live in the conda lock. Pinned by commit and installed separately (see below).

## Rebuild

```bash
# 1. Recreate the conda environment exactly (linux-64)
conda create --name seurat-env --file conda-spec.txt
conda activate seurat-env

# 2. Install SuperCell, pinned to the commit used here
Rscript -e 'remotes::install_github("GfellerLab/SuperCell@f0a041e")'
```

`conda-spec.txt` is a single-platform (linux-64) explicit spec, installed with
`conda create --file` — not the separate `conda-lock` tool. It reproduces the
conda layer exactly on linux-64.

## Key package versions

| Package    | Version | Source                          |
|------------|---------|---------------------------------|
| r-base     | 4.3.3   | conda-forge                     |
| Seurat     | 5.3.0   | conda-forge                     |
| biomaRt    | 2.58.0  | bioconda                        |
| SuperCell  | 1.0     | GitHub GfellerLab @ `f0a041e`   |
| igraph     | 2.1.4   | conda-forge                     |
| aws.s3     | 0.3.22  | conda-forge                     |

Full pinned set: see `conda-spec.txt`. Full session record: `ENVIRONMENT_sessionInfo.txt`.

## Key design decisions

- Fully conda-managed R ecosystem; **no renv**.
- `R_LIBS_USER` disabled to prevent user-library contamination.
- Compiled-heavy dependencies (igraph, nloptr, lme4, WeightedCluster, biomaRt)
  installed as conda-forge/bioconda **binaries**, not compiled from source — this
  avoids the source-build failures that occur when `install.packages`/`remotes`
  try to compile them against a fresh system.
- SuperCell installed from GitHub after the conda layer, so its dependencies are
  already present as binaries and it does not trigger a source-build cascade.

## Known issues

- **PCA error on tiny synthetic data** — too few cells/genes for the IRLBA SVD
  approximation. Fix: reduce `npcs` or use a larger dataset. Not an environment
  issue; relevant when writing small testthat fixtures.
- **tidyverse function masking (`filter`/`lag`)** — expected dplyr vs. stats
  namespace behavior, not an error.

## Validation

```r
library(Seurat)
library(tidyverse)
mat <- matrix(rpois(200, 5), nrow = 20)
obj <- CreateSeuratObject(mat)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 5)
DimPlot(obj, reduction = "pca")
```

Seurat v5 loads and the preprocessing pipeline runs cleanly in conda-managed
R 4.3.3, with no compilation fallback or missing-library errors.
# Paralogswap Environment

## System
- Conda environment: seurat-env
- R version: 4.3.3
- Seurat: 5.x (conda-forge)

## Reproducibility
All dependencies are installed via Conda.

## Activation
conda activate seurat-env
R

## Validation
Seurat pipeline test (CreateSeuratObject → RunPCA) passes successfully.
