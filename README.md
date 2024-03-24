
# scutilsR

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://img.shields.io/badge/devel%20version-0.1.0-green.svg)](https://github.com/jarninggau/scutilsR)

The package `scutilsR` is my code hub for utilities in scRNA-seq
analysis.

## Installation

You can install the development version of `scutilsR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# install DoubletFinder
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder", ref = "554097b")
# install CellChat
devtools::install_github("sqjin/CellChat", ref = "418b660")
devtools::install_github("JarningGau/scutilsR")
```

## Example

### IO

Read output from STARsolo

    Elife_2019_mouse_GSE113293/ # sample hub: contains a list of samples
    ├── GSM3102982              # sample ID:  contains a list of assays
    │   ├── corrected
    │   │   ├── barcodes.tsv
    │   │   ├── decontX.info.tsv
    │   │   ├── features.tsv
    │   │   └── matrix.mtx.gz
    │   ├── filtered            # output of STARsolo after cell calling
    │   │   ├── barcodes.tsv
    │   │   ├── features.tsv
    │   │   └── matrix.mtx.gz
    │   └── raw                 # output of STARsolo before cell calling
    │       ├── barcodes.tsv
    │       ├── features.tsv
    │       └── matrix.mtx.gz

``` r
mm <- ReadSolo(path = "/path/to/samplehub", assay = "assay")
```

Read output from [DNBelab
C4](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software)
pipeline

``` r
data <- ReadC4("/path/to/dnbc4.txt.gz") # C4 output is a compressed dense matrix
```

### Doublets Removal

Mark doublets rather than remove them via
[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder).

`MarkDoublets()` function run `DoubletsFinder` separately. All
parameters for `DoubletsFinder` are default.

- pK: auto selected by `FindOptimalpK()`
- pN: 0.25
- estimated percentage of doublets: 0.075

``` r
seu <- MarkDoublets(seu = seu, PCs = 1:10, split.by = "orig.ident")
```

Users can remove the marked doublets or clusters enriched marked
doublets.

### Remove Ambient RNAs & Calculate Contimination Rates

If you have the clusters information.

``` r
seu <- RemoveAmbientRNAs(seu, split.by = "orig.ident", cluster.name = "seurat_clusters")
```

else

``` r
seu <- RemoveAmbientRNAs(seu, split.by = "orig.ident", cluster.name = NULL)
```

### Find All Markers Parallelly

``` r
all.markers <- mcFindAllMarkers(seu = seu, do.flatten = T, only.pos = T, n.cores = 10) # returns a data.frame
all.markers <- mcFindAllMarkers(seu = seu, do.flatten = F, only.pos = T, n.cores = 10) # returns a list
```

### Annotate Cell Types According to Markers via Enrichment Analysis

collected gene sets

- [Mouse Cell Atlas (MCA)](https://doi.org/10.1016/j.cell.2018.02.001)
- [Human Cell Lanscape (HCL)](https://doi.org/10.1038/s41586-020-2157-4)

``` r
all.markers <- mcFindAllMarkers(seu.ds, do.flatten = F, only.pos = T, n.cores = 20)
all.markers <- lapply(all.markers, function(xx) subset(xx, p_val_adj < 1e-6 & avg_log2FC > log2(1.5)))

## enrichment analysis (human)
data("mca_hsa")
data("hcl_hsa")
t2g <- rbind(mca_hsa, hcl_hsa)

e.res <- pbapply::pblapply(all.markers, function(xx) {
  yy <- xx$Gene.name.uniq
  tmp <- clusterProfiler::enricher(yy, TERM2GENE = t2g)
  res <- tmp@result
  res$cluster <- xx$cluster[1]
  res
})

e.res.df <- do.call(rbind, e.res)
```

### Cell-Cell Communication

A [CellChat](https://github.com/sqjin/CellChat) wrapper
`CellChatHelper()` was implemented in `scutilsR`

``` r
DB <- subsetDB(CellChat::CellChatDB.human, search = "Secreted Signaling")

# outputs:
# - {out.dir}/cellchat.{name}.{mode}.rds
# - {out.dir}/cellchat.{name}.{mode}.pdf

CellChatHelper(seu = seu, 
               label.field = "seurat_clusters", 
               name = "run_name", 
               mode = "default", 
               DB = DB, 
               out.dir = getwd(), 
               cores = 10, 
               fig.width = 10, 
               fig.height = 8)
```
