paralogswap
===========

Detection of paralog substitution events in cross-species single-cell RNA-seq
data without cell-type annotation, native to the Seurat ecosystem.


Overview
--------

Cross-species single-cell comparison in Seurat requires reducing datasets to
one-to-one orthologs (Stuart et al. 2019; Hao et al. 2021). This step discards
paralogs and species-specific genes, losing track of paralog substitutions.
Paralog substitutions are cases where a gene's paralog tracks cross-species
expression more closely than its annotated ortholog, indicating that a cell
function has shifted to a different gene copy during evolution (Tarashansky et
al. 2021).

paralogswap detects these events by clustering each species independently,
deriving cross-species correspondence from expression, correlating each gene
against its ortholog and paralogs across matched cross-species metacells, and
flagging candidates.

The method is a focused detector for Seurat users investigating cross-species
single-cell data and flagging directions for evolutionary analysis, rather than
an atlas-alignment engine. For de novo cross-species cell and gene mapping at
large evolutionary distances, see SAMap (Tarashansky et al. 2021).

Note that `delta_r` subtracts a measured ortholog correlation from the best
paralog correlation. When the ortholog is unexpressed in the other species, that
correlation is undefined rather than zero and the gene cannot be scored. Such
genes are recorded in `attr(substitutions, "dropped")` with their paralog
correlation intact. The method is therefore least sensitive at the most complete
substitutions, where the ortholog has been lost or silenced rather than merely
diverged. Read the dropped table, not only the flagged one. See GETTING_STARTED
for scope and worked examples.


Relationship to existing tools
------------------------------

paralogswap adapts the paralog-substitution concept from SAMap (Tarashansky et
al. 2021), which learns cell and gene correspondences jointly through an
iterative graph and detects substitutions by correlating expression across that
aligned manifold. Both methods address the problem that single-cell counts are
too sparse for stable gene-gene correlation. SAMap dissolves the sparsity into a
joint, imputed cross-species manifold and measures genes over mutually-mapped
neighborhoods of individual cells. Because the cells themselves are what its
alignment matches, SAMap cannot pre-aggregate them without breaking the
correspondence it depends on. paralogswap instead resolves the sparsity by
coarse-graining each species into metacells (Bilous et al. 2022), matching
metacells across species within corresponding clusters, and correlating in gene
space over that grid. Metacells are a lightweight alternative that decouples
correspondence from measurement, keeping the pipeline transparent and inside R.
This trades some sensitivity for interpretability.

paralogswap and GeneSpectra (Song et al. 2024) both use metacells to correct
single-cell sparsity while preserving granularity, with different goals.
GeneSpectra classifies each gene by expression specificity and distribution,
then asks whether orthologs preserve their class across species, testing the
ortholog conjecture. paralogswap asks whether functional identity has swapped
from a gene to its paralog. GeneSpectra characterizes how an ortholog's
expression class is conserved or diverges; paralogswap flags the cases where a
paralog, not the ortholog, carries the conserved pattern.

Homology sourcing. paralogswap reads homology from curated Ensembl Compara calls
(Durinck et al. 2009) rather than re-deriving it from sequence alignment as
SAMap does, trading coverage on poorly-annotated genomes for transparency,
reproducibility, and speed. It keeps the full homology graph, including all
one-to-many and many-to-many edges, following the benchmark categorization of
Song et al. 2023 and SAMap's full-graph approach, rather than reducing to
one-to-one orthologs as in standard correlation analysis.

Curated orthology is fallible. Ensembl Compara release 116 mis-assigns mouse
Ramp2 to the VPS25 orthogroup, an adjacent-gene reconciliation error verifiable
in both query directions. Because orthology and paralogy are read from the same
gene tree, the error severs the gene's entire homology neighbourhood rather than
a single edge: Compara reports mouse Ramp2's only paralog as Vps25, and does not
consider it a paralog of Ramp1 or Ramp3. Repairing the gene therefore requires
restoring its paralogy alongside its ortholog.

paralogswap accepts user-supplied corrections through the `corrections`
argument, merged with provenance in the `source` column. The RAMP2 edges are
recovered from eggNOG orthogroups (Huerta-Cepas et al. 2019), which group the
gene correctly. An `action` column takes "add" or "remove". Removals are applied
before additions, match
irrespective of edge orientation, and are retained in `attr(x, "excluded")` with
their original source, so an override is recorded rather than hidden. Removal is
necessary rather than convenient: `detect_substitutions` takes `max(r_paralog)`,
so supplying correct edges cannot out-vote an incorrect one, since a spurious
paralog may win the maximum by chance.

Removal is declared rather than inferred. SAMap prunes its homology graph to
positively correlated gene pairs, which is available to a method that aligns
manifolds before measuring genes. paralogswap decouples correspondence from
measurement, so pruning edges by the correlations detection depends on would
tune the apparatus on the result.

The package fetches Ensembl automatically but does not attempt to detect or
repair such errors. `min_perc_id` is an optional absolute identity floor, rather
than SAMap's per-gene relative 0.25 cutoff, which prunes a BLAST bit-score graph
not constructed in this pipeline.


Installation
------------

    # SuperCell is distributed on GitHub only and must be installed first;
    # it is not on CRAN or Bioconductor and will not resolve as a dependency.
    remotes::install_github("GfellerLab/SuperCell")

    # biomaRt is on Bioconductor
    BiocManager::install("biomaRt")

    remotes::install_github("claratucker/paralogswap")

Imports: Seurat (Hao et al. 2021), biomaRt (Durinck et al. 2009), SuperCell
(Bilous et al. 2022), ggplot2, Matrix, rlang. Suggests: ggrepel.

A pinned conda environment is provided. See ENVIRONMENT.md for the exact
rebuild, including the SuperCell commit pin and the compiled dependencies
(WeightedCluster, weights) that must be installed as binaries rather than built
from source.


Quick start
-----------

The pipeline runs stage by stage. Two homology graphs are built per species
pair: a broad ortholog-only graph that bridges the matching stages, and a focal
graph carrying the paralog edges detection compares against. `genes` bounds the
paralog pull, which is otherwise genome-wide and will exceed biomaRt's request
ceiling.

    library(paralogswap)

    focal <- c("RAMP1", "RAMP2", "RAMP3")

    # Stage 0 -- homology. Species are biomaRt dataset prefixes.
    hg_broad <- build_homology_graph(
      species_a = "mmurinus", species_b = "hsapiens",
      genes = rownames(lemur_lung), id_type = "symbol",
      ensembl_release = 116, include_paralogs = FALSE
    )
    hg_focal <- build_homology_graph(
      species_a = "mmurinus", species_b = "hsapiens",
      genes = focal, id_type = "symbol",
      ensembl_release = 116, include_paralogs = TRUE
    )

    # Stage 1 -- independent clustering, no cross-species information
    sc <- cluster_species(list(lemur = lemur_lung, human = human_lung))

    # Stage 2 -- coarse cross-species cluster correspondence, from expression
    cc <- match_clusters(sc[["lemur"]], sc[["human"]], hg_broad,
                         species_a = "mmurinus", species_b = "hsapiens")

    # Stage 3 -- metacells, one species at a time
    mc_lemur <- build_metacells(sc[["lemur"]], gamma = 20)
    mc_human <- build_metacells(sc[["human"]], gamma = 20)

    # Stage 4 -- fine metacell correspondence within matched clusters
    mmc <- match_metacells(mc_lemur, mc_human, cc, hg_broad,
                           species_a = "mmurinus", species_b = "hsapiens")

    # Stages 5-6 -- correlate homologs, then flag substitutions
    hc <- compute_homolog_correlations(mc_lemur, mc_human, mmc, hg_focal,
                                       species_a = "mmurinus",
                                       species_b = "hsapiens",
                                       focal_genes = focal)
    subs <- detect_substitutions(hc)

    plot_ortholog_vs_paralog(hc)   # r_ortholog vs max r_paralog, with diagonal
    plot_substitution_elbow(subs)  # sorted delta_r; read the set off the bend
    attr(subs, "dropped")          # focal genes that could not be scored

`lemur_lung` and `human_lung` are user-supplied Seurat objects; see
data-raw/load_lemur_lung.R and data-raw/load_human_lung.R. Each stage returns an
inspectable object with a `summary()` method. data-raw/run_lemur_human.R is the
full worked pipeline, including persistence of every intermediate.

Note that the method is directional. `compute_homolog_correlations` scores focal
genes of species A against homologs in species B, so a run flags genes of its
anchor species only. A gene residing on the other species requires a run
anchored there, which may reuse the same metacells and grid.


Pipeline
--------

`build_homology_graph` performs ortholog and paralog assignment from Ensembl
Compara via biomaRt (Durinck et al. 2009), pulled in both directions. The full
graph is retained, and user-supplied corrections are added or removed with
provenance.

`cluster_species` clusters each species independently. Resolution defaults to
the Seurat convention of 0.8.

`match_clusters` derives coarse, interpretable cross-species cluster
correspondence from expression.

`build_metacells` performs SuperCell coarse-graining per species (Bilous et al.
2022).

`match_metacells` derives fine cross-species metacell correspondence within
matched clusters.

`compute_homolog_correlations` computes gene-space correlation per homolog pair
over the metacell grid.

`detect_substitutions` flags candidates where a paralog correlates better than
the ortholog, and records genes that cannot be scored.

`assess_resolution_stability` re-runs the pipeline across clustering resolutions
and graining levels, and reports `delta_r` and `flagged` at each setting (in
progress).

Diagnostic plots: `plot_substitution_elbow`, `plot_ortholog_vs_paralog`,
`plot_graining_curve`, `plot_correspondence_heatmap` (in progress). Use
`Seurat::ElbowPlot()` for the PCA scree, and `Seurat::DotPlot()` on metacells to
show a gene and its homologs across matched clusters.


Scope and assumptions
---------------------

Correspondence. Cross-species correspondence is derived from expression, not
supplied as a label.

Directionality. A run anchored on species A can flag genes of species A only.
Multi-species results are composed from runs anchored appropriately, and an
anchor determines not only which genes may be flagged but which are measurable
at all.

Full homology graph. All ortholog and paralog edges are retained and carry their
Ensembl `ortholog_type`. Genes in a many-to-many set are flagged `ambiguous`.
Nothing is silently reduced to one-to-one orthologs. `detect_substitutions`
reports `n_orthologs` per focal gene, so a gene whose `delta_r` depends on which
of several orthologs was taken is visible rather than hidden. Note that the
matching stages bridge on one-to-one orthologs only, so many-to-many edges enter
detection but not correspondence.

Gene identifiers. The homology graph and the expression objects must share an
identifier space, since edges join to rownames by string match and no conversion
is performed. `id_type = "ensembl"` builds the graph on Ensembl gene IDs;
`id_type = "symbol"` builds it on gene symbols. Atlases distributed with symbol
rownames, which is common, therefore require `id_type = "symbol"`, and any
`corrections` table must be written in the same space as the graph it joins.
Symbol conventions differ across species (human RAMP2, mouse Ramp2) and must
match the source annotation. Unnamed placeholder loci are common in
less-annotated genomes and are expected to drop out; the package reports the
unmapped rate.

Pairwise comparison. The current implementation compares two species per run.
Multi-species analysis runs as pairwise calls composed afterward, mirroring
SAMap's pairwise-plus-transitivity design (Tarashansky et al. 2021). The
validation vignette composes a three-species triad this way. A future release
will add a joint multi-species test that accounts for phylogeny, following the
strategies benchmarked in Song et al. 2023 and Wang et al. 2025.

Measurement substrate. Single-cell counts are too sparse for stable gene-gene
correlation, but collapsing each matched cluster to a single pseudobulk profile
leaves only a handful of observations per gene. paralogswap coarse-grains each
species into metacells with SuperCell (Bilous et al. 2022), matches metacells
across species within each corresponding cluster, and correlates over that fine
grid in gene space, giving many paired observations per gene while preserving
the manifold (Squair et al. 2021 motivates aggregation; SuperCell supplies the
granularity). In the validation, `gamma = 20` grains 34,084 lemur cells into
1,704 metacells, 40,000 human cells into 2,000, and 24,540 mouse cells into
1,228.

The graining level `gamma` trades resolution against robustness. It is chosen by
the user, by reading `plot_graining_curve` and taking the largest gamma still on
the correlation-preservation plateau. The package applies no automatic
selection, and `build_metacells` defaults to `gamma = 20` if none is given. The
curve reports only single-species fidelity and touches neither the homology
graph, the second species, nor the ortholog and paralog distinction, so it is
computable before detection exists and cannot be tuned by the number of
substitutions found.

Metacell matching. Cross-species metacell pairs are formed within reciprocally
matched clusters by mutual k-nearest-neighbor correlation on one-to-one
orthologs, not by reciprocal-best matching, which lets a few high-expression
metacells absorb the matches and collapses the grid. The mutual condition is
symmetric, so a grid computed with one species as anchor may be relabelled for a
run anchored on the other without recomputation. This assumes the cluster
correspondence is approximately correct, which the annotation check tests.

Clustering sensitivity. Derived correspondence depends on clustering resolution
and graining level. `assess_resolution_stability` re-runs the pipeline across a
small set of resolutions and gammas and reports `delta_r` and `flagged` at each
setting, with `delta_threshold` held fixed so the settings remain comparable (in
progress). Values are reported rather than only overlap counts, since a call
that clears the threshold narrowly at one setting is the case this check exists
to reveal.

Normalization. Log-normalization is the default (`normalization = "lognorm"`,
Seurat's `NormalizeData`). SCTransform is available (Hafemeister & Satija 2019).

Threshold. The default substitution threshold (`delta_r > 0.3`) matches
Tarashansky et al. 2021, so outputs are comparable to published results. The
threshold is a default, not a fixed rule. It is recorded on the returned object,
`plot_substitution_elbow` draws the line detection used against the ranked
`delta_r` distribution, and `plot_ortholog_vs_paralog` renders each gene's
paralog-versus-ortholog correlation directly, so the decision stays with the
user, following the `ElbowPlot` idiom in Seurat v5 (Hao et al. 2021).


Validation
----------

paralogswap reproduces a published finding from independent public data. The
primary validation is a three-species lung comparison: mouse (Tabula Muris
Senis, Almanzar et al. 2020), mouse lemur (Tabula Microcebus, Ezran et al.
2025), and human (Tabula Sapiens, Tabula Sapiens Consortium 2022).

The target is the RAMP expression-homologue triads reported by Ezran et al. 2025
in lung: human RAMP2, lemur RAMP1, mouse Ramp2; and human RAMP3, lemur RAMP2,
mouse Ramp2. Lemur RAMP1 and human RAMP3 are the evolutionary outliers, both
resembling the conserved RAMP2 pattern in lung endothelium.

The validation runs two lemur-anchored pairwise comparisons, plus a third focal
run anchored on human. The third is required because the second triad places its
outlier on human, and a lemur-anchored run cannot flag it.

Status of each validation target:

First triad, lemur to human: reproduced. Bridging on 13,518 one-to-one orthologs
and 192 mutual-kNN metacell pairs across 7 reciprocally matched cluster pairs,
lemur RAMP1 is flagged with r_ortholog = -0.058, r_paralog = 0.544, delta_r =
0.602, best paralog human RAMP2. Lemur RAMP2 is not flagged, with r_ortholog =
0.625 exceeding r_paralog = 0.486.

First triad, lemur to mouse: reproduced. Bridging on 11,873 one-to-one orthologs
and 203 mutual-kNN metacell pairs across 8 reciprocally matched cluster pairs,
lemur RAMP1 is flagged with r_ortholog = -0.595, r_paralog = 0.297, delta_r =
0.891, best paralog mouse Ramp2. This arm requires the Compara correction
described above, supplied as three added edges and two removals. The
substitution call does not depend on the correction, but the partner identity
does: without it, mouse Ramp2 is unreachable as a paralog and `best_paralog`
resolves to Ramp3 at 0.282, a margin of 0.015. Lemur RAMP2 is not flagged, but
at delta_r = 0.055 it clears the threshold narrowly, with mouse Ramp3 at 0.317
exceeding the ortholog at 0.262. The mouse arm supports the conserved anchor at
this threshold, not on the ordering of correlations.

Note that the two arms share the lemur clustering but not their reciprocal
cluster pairs, since reciprocity is a property of the species pair rather than
of lemur alone: 7 pairs against human, 8 against mouse. Restricting matching to
reciprocally matched clusters also leaves metacells unmatched. Of 1,704 lemur
metacells, 706 find no mouse partner, so the mouse grid rests on roughly 59
percent of the lemur metacells.

First triad closure: reproduced. Lemur RAMP1 is flagged in both arms and the
paralog carrying the conserved pattern is RAMP2 in each. The two runs share the
lemur atlas and its clustering, but the other species' data, homology graphs,
cluster correspondences, and metacell grids are independent, so closure is not
an artefact of a shared alignment. This is a single-instance transitivity check
on one gene family, following the criterion of Tarashansky et al. 2021, not a
distributional statistic.

Second triad: recovered, not flagged. Anchored on human, RAMP3 tracks lemur
RAMP2 at r = 0.486, above lemur RAMP1 at 0.407, and does not track its own
ortholog. This is the reported relationship. It is not a flagged substitution:
lemur RAMP3 shows no variance across the grid, expressed in 8 of 1,704 lemur
metacells atlas-wide, so delta_r has no baseline to subtract.
`detect_substitutions` records the gene in `attr(x, "dropped")` with reason
"ortholog_invariant". The same gene is recorded as "focal_gene_invariant" in the
lemur-anchored arms, where it is the focal gene and no correlation is defined at
all. The second triad therefore does not close on flags, and cannot.

One-to-one-only control: structurally satisfied. Every RAMP ortholog edge is
one-to-one (n_orthologs = 1 in all three runs), and the matching stages bridge
on one-to-one orthologs by construction, so spurious many-to-many linkage cannot
have produced this result. The control becomes an empirical test only for a
genome-wide run containing focal genes with n_orthologs > 1.

Cluster correspondence against published Cell Ontology annotations, used only as
a check on the correspondence step: pending.

Persistence of flagged substitutions across clustering resolutions: not yet run.
This matters most for the mouse arm, where lemur RAMP2 clears the conserved
anchor threshold by 0.055.

See vignettes/validation.Rmd.

Future work. Evolutionary stratification by paralog age requires a substitution
set spanning thousands of genes, and therefore a chunked genome-wide paralog
pull. It also requires care: because delta_r is undefined when the ortholog has
been silenced, and loss of one duplicate is a common fate of recent paralogs, a
rate estimated from scoreable genes alone would be biased against the hypothesis
it tests. `duplication_level` and `dn_ds` will be carried on paralog edges so the
duplication node of a flagged gene can be reported without the rate (in
progress). Additional
tissues, compartment-level enrichment, functional stratification, and the joint
multi-species test are also future work.


Troubleshooting
---------------

Ensembl's BioMart service is intermittently unavailable, and biomaRt has no
retry. A genome-wide pull may exceed its fixed request ceiling, and requests may
be redirected to status.ensembl.org. Pass `mirror = "useast"` or "asia" to route
around the main host. Note that biomaRt does not permit a mirror and a pinned
`ensembl_release` in the same call; when a mirror is used, the current release is
served. Bound large pulls with `genes` rather than relying on
`include_paralogs = FALSE` alone. See ENVIRONMENT.md.


Data
----

All datasets are public. Tabula Muris Senis: Figshare 8273102 (Almanzar et al.
2020). Tabula Sapiens: Figshare project 100973 (Tabula Sapiens Consortium 2022).
Tabula Microcebus: Figshare project 112227 (Ezran et al. 2025), released under a
data-use policy reproduced in full in DATA_POLICY.md.


License
-------

MIT. Data-use policy for Tabula Microcebus reproduced in DATA_POLICY.md.


References
----------

Almanzar, N. et al. 2020. A single-cell transcriptomic atlas characterizes
ageing tissues in the mouse. Nature 583:590.

Bilous, M. et al. 2022. Metacells untangle large and complex single-cell
transcriptome networks. BMC Bioinformatics 23:336.

Durinck, S. et al. 2009. Mapping identifiers for the integration of genomic
datasets with the R/Bioconductor package biomaRt. Nat Protoc 4:1184.

Ezran, C. et al. 2025. Mouse lemur cell atlas informs primate genes, physiology
and disease. Nature 644:185.

Hafemeister, C. and Satija, R. 2019. Normalization and variance stabilization of
single-cell RNA-seq data using regularized negative binomial regression. Genome
Biol 20:296.

Hao, Y. et al. 2021. Integrated analysis of multimodal single-cell data. Cell
184(13):3573-3587.e29.

Huerta-Cepas, J. et al. 2019. eggNOG 5.0: a hierarchical, functionally and
phylogenetically annotated orthology resource. Nucleic Acids Res 47:D309.

Song, Y., Miao, Z., Brazma, A. et al. 2023. Benchmarking strategies for
cross-species integration of single-cell RNA sequencing data. Nat Commun 14:6495.

Song, Y. et al. 2024. Revising the ortholog conjecture in cross-species
comparison of scRNA-seq data. bioRxiv 2024.06.21.600109.

Squair, J. W. et al. 2021. Confronting false discoveries in single-cell
differential expression. Nat Commun 12:5692.

Stuart, T. et al. 2019. Comprehensive integration of single-cell data. Cell
177:1888.

Tabula Sapiens Consortium. 2022. The Tabula Sapiens: a multiple-organ,
single-cell transcriptomic atlas of humans. Science 376:eabl4896.

Tarashansky, A. J. et al. 2021. Mapping single-cell atlases throughout Metazoa
unravels cell type evolution. eLife 10:e66747.

Wang, S. et al. 2025. Benchmarking cross-species single-cell RNA-seq data
integration methods: towards a cell type tree of life. Nucleic Acids Res
53:gkae1316.
