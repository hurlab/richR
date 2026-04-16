# richR <img src="man/figures/logo.png" align="right" height="120" />

[![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://img.shields.io/badge/version-0.1.4-green.svg)](https://github.com/hurlab/richR)
[![R-CMD-check](https://github.com/hurlab/richR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hurlab/richR/actions/workflows/R-CMD-check.yaml)
![](https://img.shields.io/github/languages/code-size/hurlab/richR)
[![DOI](https://zenodo.org/badge/243827597.svg)](https://zenodo.org/badge/latestdoi/243827597)

## Overview

**richR** is an R package for functional enrichment analysis and visualization.
It supports Over-Representation Analysis (ORA) via hypergeometric testing,
Gene Set Enrichment Analysis (GSEA), and kappa-score-based term clustering.
Built-in annotation builders cover GO, KEGG, Reactome, KEGG Module, and MSigDB,
and you can supply custom gene sets via GMT files or named lists.

### Key features

- **Enrichment**: `richGO()`, `richKEGG()`, `enrich()`, `richGSEA()`, `richDAVID()`
- **Clustering**: `richCluster()` groups terms by kappa similarity
- **20+ visualization functions** including bar, dot, lollipop, heatmap, network,
  volcano, circle bar, gene-dot, gene-heat, gene-bar, UpSet, and more
- **dplyr support**: `filter()`, `select()`, `mutate()`, `arrange()`, `group_by()`,
  `slice()`, `summarise()` work directly on enrichment result objects
- **Batch analysis**: `batchEnrich()` runs enrichment across multiple gene lists
- **Comparison**: `compareResult()` and `richCompareDot()` for multi-group comparisons
- **Publication-quality UpSet plots** with pixel-accurate panel alignment and
  intuitive color-blended multi-set bars (`richUpset()`)

## Installation

```r
# Install from GitHub
library(devtools)
install_github("hurlab/richR")
```

Bioconductor annotation packages are needed for building annotations:
```r
BiocManager::install(c("org.Hs.eg.db", "GO.db"))       # human GO
BiocManager::install("reactome.db")                      # Reactome (optional)
```

## Quick Start

### 1. Build annotation data

```r
library(richR)

# Check available species
showData()

# Build GO and KEGG annotations
hsago <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "GO")
hsako <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

# KEGG Module
hsakom <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGGM")

# MSigDB
hsamgi <- buildMSIGDB(species = "human", keytype = "SYMBOL", anntype = "GO")
```

### 2. Run enrichment analysis

```r
# GO enrichment (ORA)
gene <- sample(unique(hsago$GeneID), 1000)
resgo <- richGO(gene, godata = hsago, ontology = "BP")
head(resgo)

# KEGG enrichment
resko <- richKEGG(gene, kodata = hsako, pvalue = 0.05)

# GSEA (requires a named numeric vector)
genelist <- rnorm(1000)
names(genelist) <- sample(unique(hsako$GeneID), 1000)
res_gsea <- richGSEA(genelist, object = hsako)

# DAVID (online)
res_david <- richDAVID(gene, keytype = "ENTREZID", species = "human")
```

### 3. Custom gene sets

```r
# Import from GMT file (MSigDB, Enrichr, etc.)
annot <- readGMT("h.all.v2023.2.Hs.symbols.gmt", species = "human")
res <- enrich(my_genes, annot)

# From a named list
my_sets <- list(
  "Apoptosis"  = c("TP53", "BAX", "BCL2", "CASP3"),
  "Cell Cycle" = c("CDK1", "CDK2", "CCND1", "RB1")
)
annot <- buildAnnotFromList(my_sets)
res <- enrich(my_genes, annot)
```

### 4. Custom annotation with bioAnno

```r
# library(bioAnno)
# fromKEGG(species = "ath")
# athgo <- buildOwn(dbname = "org.ath.eg.db", anntype = "GO")
# athko <- buildOwn(dbname = "org.ath.eg.db", anntype = "KEGG")
# See https://github.com/guokai8/bioAnno for details
```

## Visualization

richR provides **20+ visualization functions**. All accept `richResult`,
`GSEAResult`, or `data.frame` objects and support saving to file via
`filename`, `width`, and `height` arguments. Each function has a `rich*`
primary name and backward-compatible `gg*` aliases.

### Term-level summary plots

Classic plots that rank and display enriched terms:

```r
# Bar plot — horizontal bars of gene count per term, colored by significance
richBar(resgo, top = 20, usePadj = FALSE)

# Dot plot — bubble plot: x = RichFactor, y = Term, size = gene count, color = -log10(p)
richDot(resko, top = 10, usePadj = FALSE)

# Lollipop — clean alternative to bar plot: x = RichFactor, color = significance
richLollipop(resko, top = 20)

# Circular bar plot — polar-coordinate bar chart of enrichment terms
richCircle(resko, top = 15)
```

<p align="center">
  <img src="man/figures/example_barplot.png" width="48%" alt="Bar plot example" />
  <img src="man/figures/example_dotplot.png" width="48%" alt="Dot plot example" />
</p>

### Scatter and volcano plots

Two-axis continuous plots for exploring enrichment structure:

```r
# Scatter — enrichment funnel plot (MA-plot analogy for enrichment)
#   x = log10(pathway size), y = log2(fold enrichment), color = -log10(p)
#   Small pathways (left) can have extreme FE by chance; large pathways
#   (right) with high FE are the most robust findings
richScatter(resgo, top = 50, label.top = 5)
richScatter(resko, usePadj = FALSE)

# Volcano — enrichment volcano plot
#   x = effect size (log2 FE for ORA, NES for GSEA)
#   y = -log10(p), color = effect (diverging gradient), size = -log10(p)
#   Labels auto-adjusted via ggrepel; works for both ORA and GSEA
richVolcano(resgo)                               # ORA: x = log2(Fold Enrichment)
richVolcano(res_gsea)                            # GSEA: x = NES
richVolcano(resgo, label.top = 10, short = TRUE) # customize labels

# Term similarity scatter — MDS projection of gene-set overlap (Jaccard)
#   Functionally related terms cluster together in 2D space
richTermSim(resko, top = 20)
```

### Gene-level plots

Show which genes drive the enrichment:

```r
# Gene-term dot plot — color by fold change (if provided) or -log10(p)
richGeneDot(resgo, fc = my_foldchanges)
richGeneDot(resgo)                              # no fc: color by -log10(Pvalue)

# Gene-term heatmap — pheatmap-style tile plot with borders
#   Missing gene-term combinations shown as grey; supports fold-change coloring
richGeneHeat(resgo, fc = my_foldchanges)
richGeneHeat(resgo, na.fill = "grey90", border.color = "grey40")

# Gene bar plot — bars per term with gene labels aligned inside
richGeneBar(resgo, top = 10)
```

### GSEA-specific plots

```r
# NES bar plot — up/down color-coded by enrichment direction
richNES(res_gsea, top = 20)

# Running enrichment score curve (classic GSEA plot)
richGSEAcurve(object = hsako, gseaRes = res_gsea)
richGSEAcurve(object = hsako, gseaRes = res_gsea, pathways = "hsa04110")

# ECDF step plot — cumulative rank distribution of gene sets
richECDF(res_gsea)
```

![GSEA enrichment curve](man/figures/example_gsea.png)

### Network and map plots

```r
# Gene-concept network — bipartite graph of terms and their genes
richNetplot(resko, top = 20)

# Term similarity network — edges weighted by gene-set overlap (kappa)
richNetwork(resgo, top = 20, weightcut = 0.01)

# Combined network map — overlay multiple enrichment results
richNetmap(list(resgo, resko), top = 50)
```

### Cluster and comparison plots

```r
# Kappa-based clustering of enrichment terms
resc <- richCluster(resgo)
richClusterDot(resc)

# Multi-group comparison
resko1 <- richKEGG(gene1, kodata = hsako)
resko2 <- richKEGG(gene2, kodata = hsako)
res_cmp <- compareResult(list(S1 = resko1, S2 = resko2))
richCompareDot(res_cmp)

# Comparison heatmap across multiple analyses
richHeatmap(list(GO = resgo, KEGG = resko), top = 50)
```

![Multi-group comparison dot plot](man/figures/example_comparedot.png)

### UpSet plot

`richUpset()` renders an UpSet-style plot (intersection bar chart +
dot matrix + set-size bars) for visualizing overlaps between named
gene lists. Built entirely with ggplot2 — no dependency on the UpSetR
package. Key features:

- **Pixel-accurate bar-to-dot alignment** — the intersection bar chart
  and the dot matrix are pre-aligned with `cowplot::align_plots()` and
  positioned via `ggdraw()` + `draw_plot()`, so the bars always line
  up with the dots beneath them regardless of y-axis label widths.
- **Blended multi-set bar colors** — single-set bars keep the set
  color; multi-set intersection bars are rendered in an RGB blend of
  the participating set colors, making each combination visually
  distinct. Set `main.bar.color` to force a uniform color.
- **Flexible input** — named list of character vectors, list of
  `richResult` / `GSEAResult` objects, or data frames with a `GeneID`
  / `leadingEdge` column.

```r
set.seed(123)
all_genes <- paste0("GENE", seq_len(250))
gene_lists <- list(
  "Treatment A" = sample(all_genes, 150),
  "Treatment B" = sample(all_genes, 130),
  "Control"     = sample(all_genes, 110),
  "Baseline"    = sample(all_genes, 100)
)

richUpset(gene_lists)

# Custom colors and ordering
richUpset(gene_lists,
          mycol = c("dodgerblue", "goldenrod1", "seagreen3", "orchid3"),
          order.by = "degree", nintersects = 20)

# Uniform bar color override
richUpset(gene_lists, main.bar.color = "steelblue")
```

![UpSet plot example](man/figures/example_upset.png)

### Plot summary table

| Function | Plot type | x-axis | y-axis | Color | Size |
|---|---|---|---|---|---|
| `richBar` | Bar | Gene count | Term | -log10(p) | — |
| `richDot` | Bubble | RichFactor | Term | -log10(p) | Gene count |
| `richLollipop` | Lollipop | RichFactor | Term | -log10(p) | — |
| `richCircle` | Circular bar | Gene count | Term (polar) | -log10(p) | — |
| `richScatter` | Funnel/MA | log10(Pathway Size) | log2(Fold Enrichment) | -log10(p) | — |
| `richVolcano` | Volcano | Effect (log2 FE / NES) | -log10(p) | Effect | -log10(p) |
| `richTermSim` | MDS scatter | Dim 1 | Dim 2 | -log10(p) | Gene count |
| `richNES` | Bar (GSEA) | NES | Term | Up/Down | — |
| `richECDF` | Step | NES rank | Cumulative fraction | — | — |
| `richGeneDot` | Dot matrix | Gene | Term | FC / -log10(p) | — |
| `richGeneHeat` | Heatmap | Gene | Term | FC / presence | — |
| `richGeneBar` | Stacked bar | RichFactor | Term | Gene labels | — |
| `richNetplot` | Network | — | — | — | — |
| `richNetwork` | Network | — | — | -log10(p) | Gene count |
| `richUpset` | UpSet | Intersection | Set size | Group | — |

## Working with results

```r
# Extract result table and detail table
result(resgo)
detail(resgo)

# Get genes for specific terms
getGenes(resgo, term = "GO:0006955")

# dplyr verbs work directly
library(dplyr)
resgo %>% filter(Padj < 0.01) %>% head()
resgo %>% select(Term, Pvalue, Padj)
resgo %>% arrange(Pvalue) %>% head(10)

# Batch enrichment across multiple gene lists
gene_lists <- list(GroupA = gene1, GroupB = gene2)
batch_res <- batchEnrich(gene_lists, annot = hsago)
```

## Function name aliases

All plotting functions use the `rich*` prefix as the primary name.
Old `gg*` names are kept as backward-compatible aliases:

| Primary name | Aliases |
|---|---|
| `richBar` | `ggbar` |
| `richDot` | `ggdot` |
| `richLollipop` | `gglollipop`, `ggLollipop` |
| `richGeneDot` | `gggenedot`, `ggGeneDot` |
| `richGeneHeat` | `gggeneheat`, `ggGeneHeat` |
| `richGeneBar` | `gggenebar`, `ggGeneBar` |
| `richCircle` | `ggcircbar`, `ggCircBar` |
| `richVolcano` | `ggvolcano`, `ggVolcano` |
| `richNES` | `ggNES`, `ggnes` |
| `richScatter` | `ggscatter`, `ggScatter` |
| `richECDF` | `ggecdf`, `ggECDF` |
| `richTermSim` | `ggtermsim`, `ggTermSim` |
| `richNetplot` | `ggnetplot` |
| `richNetwork` | `ggnetwork` |
| `richNetmap` | `ggnetmap` |
| `richUpset` | `ggupset` |
| `richHeatmap` | `ggheatmap` |
| `richClusterDot` | `ggcluster` |
| `richCompareDot` | `comparedot` |
| `richGSEAplot` | `ggGSEA` |
| `richGSEAcurve` | `plotGSEA` |

## Citation

If you use richR in your research, please cite:

```
Guo K and Hur J (2026). richR: Functional enrichment analysis and visualization.
R package version 0.1.4. https://github.com/hurlab/richR
```

Upstream releases and their historical DOIs are maintained at
[guokai8/richR](https://github.com/guokai8/richR); the hurlab fork
tracks upstream and ships fix-ahead patches (currently the
`richUpset()` alignment fix and multi-set color blending, proposed
upstream as PR #9 — pending merge at time of writing).

## Contact

For questions or bug reports on the hurlab fork, please open an
[issue on GitHub](https://github.com/hurlab/richR/issues).
For issues affecting upstream, please file at
[guokai8/richR](https://github.com/guokai8/richR/issues).
