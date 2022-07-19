#  Mark Doublets via DoubletFinder
#' @import Seurat
#' @import DoubletFinder
NULL

PreprocessSeurat <- function(seu, PCs=1:10) {
  message("  Normalize ...")
  seu <- NormalizeData(seu, verbose = F)
  message("  Find variable genes ...")
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = F)
  message("  Scale data ...")
  seu <- ScaleData(seu, verbose = F)
  message("  Run PCA ...")
  seu <- RunPCA(seu, verbose = F)
  message("  Find neighbors ...")
  k <- round(ncol(seu)*0.02)
  k <- ifelse(k < 20, 20, k)
  seu <- FindNeighbors(seu, dims = 1:10, reduction = "pca", k.param = k, verbose = F)
  message("  Find clusters ...")
  seu <- FindClusters(seu, resolution = 0.4, verbose = F)
  seu
}

FindOptimalpK <- function(seu, PCs=1:10, num.cores = 8) {
  sweep.res.list <- paramSweep_v3(seu, PCs = PCs, sct = F, num.cores = num.cores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn[which.max(bcmvn$BCmetric), ]$pK
  as.numeric(as.character(pK))
}


DF <- function(seu, PCs=1:10, auto.pK=TRUE, auto.cluster=TRUE) {
  message("1. Preprocessing ...")
  # this step is for getting clusters
  if (auto.cluster) {
    seu <- PreprocessSeurat(seu, PCs = PCs)
  }

  message("2. Find optimal pK ...")
  if (auto.pK) {
    optimal.pK <- FindOptimalpK(seu, PCs = PCs, num.cores = 8)
    message(paste0("   Optimal pK = ", optimal.pK))
  } else {
    optimal.pK <- 0.09 # default
  }

  message("3. Mark doublets ...")
  homotypic.prop <- modelHomotypic(seu$seurat_clusters)
  nExp_poi <- round(0.075*ncol(seu))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  pN = 0.25 # default 0.25
  pANN_name <- paste("pANN", pN, optimal.pK, nExp_poi, sep = "_")
  DF_name <- paste("DF.classifications", pN, optimal.pK, nExp_poi, sep = "_")
  DF_name.adj <- paste("DF.classifications", pN, optimal.pK, nExp_poi.adj, sep = "_")

  seu <- doubletFinder_v3(seu, PCs = PCs, pN = pN, pK = optimal.pK, nExp = nExp_poi, reuse.pANN = F, sct = F)
  seu <- doubletFinder_v3(seu, PCs = PCs, pN = pN, pK = optimal.pK, nExp = nExp_poi.adj, reuse.pANN = pANN_name, sct = F)
  # summarize results
  seu$DF.classifications <- ifelse(seu@meta.data[[DF_name.adj]] == "Doublet", "Doublet.hc",
                                   ifelse(seu@meta.data[[DF_name]] == "Doublet", "Doublet.lc", "Singlet"))
  seu@meta.data[[pANN_name]] <- NULL
  seu@meta.data[[DF_name]] <- NULL
  seu@meta.data[[DF_name.adj]] <- NULL
  seu
}

#' Mark doublets
#' @param seu Seurat object
#' @param PCs Vectors indicating used priciple components. Default: 1:10
#' @param split.by Name of a metadata column to split plot by. Default: NULL
#' @return Seurat object
#' @export
MarkDoublets <- function(seu, PCs=1:10, split.by=NULL) {
  if (is.null(split.by)) {
    seu.list <- list(seu)
  } else {
    seu.list <- SplitObject(seu, split.by = split.by)
  }
  names(seu.list) <- NULL
  metadata.new <- pbapply::pblapply(seu.list, function(xx) {
    seu.tmp <- DF(xx, PCs=PCs, auto.pK=TRUE, auto.cluster=TRUE)
    seu.tmp@meta.data
  })
  metadata.new <- do.call(rbind, metadata.new)
  seu$DF.classifications <- metadata.new[rownames(seu@meta.data), ]$DF.classifications
  return(seu)
}
