#' A wrapper of celda::decontX
#' @param seu Seurat object
#' @param split.by The grouping variable used to split the input object. Default: NULL.
#' @param cluster.name Name of cluster field in Seurat object. Default: NULL.
#' @return seurat object including a corrected counts matrix and contamination rates for each cell.
#' @export
RemoveAmbientRNAs <- function(seu, split.by = NULL, cluster.name = NULL) {
  if (is.null((cluster.name))) {
    stop("The clusters must be pre-defined.")
  }
  sce <- sceasy::convertFormat(seu, from="seurat", to="sce")
  if (is.null(split.by)) {
    sce <- celda::decontX(sce, z=sce[[cluster.name]])
  } else {
    sce <- celda::decontX(sce, z=sce[[cluster.name]], batch=sce[[split.by]])
  }
  seu[["decontX"]] <- CreateAssayObject(counts = round(celda::decontXcounts(sce), 0))
  seu$decontX_contamination <- SingleCellExperiment::colData(sce)$decontX_contamination
  return(seu)
}
