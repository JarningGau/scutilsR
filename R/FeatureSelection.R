message2 <- function(msg, show.time=T) {
  if (show.time) {
    message(sprintf("%s %s", Sys.time(), msg))
  } else {
    message(msg)
  }
}

#' Find Variable Features from Integrated Dataset
#'
#' @param object A Seurat object.
#' @param split.by The grouping variable used to split the input object.
#' @param do.norm A logical indicating whether the data should be normalized. Default is TRUE.
#' @param selection.method The method used to select variable features. Default is "vst".
#' @param nfeatures The number of variable features to select. Default is 2000.
#' @param min.cells The number of minimal cells for calculation. Default is 500.
#' @param ... Additional arguments passed to FindVariableFeatures.
#'
#' @return A vector of selected feature names.
#'
#' @importFrom Seurat SplitObject NormalizeData FindVariableFeatures SelectIntegrationFeatures
#'
#' @examples
#' \dontrun{
#' # Find variable features
#' vfeatures <- FindIntVarFeatures(object = seu, split.by = "Batch")
#' }
#'
#' @export
FindIntVarFeatures <- function(object, split.by, do.norm = TRUE,
                               selection.method = "vst", nfeatures = 2000,
                               min.cells = 500, ...) {
  message2("Split object ...")
  seu.list <- SplitObject(object, split.by = split.by)
  for (i in 1:length(seu.list)) {
    if (ncol(seu.list[[i]]) >= min.cells) {
      message2(sprintf("Find HVGs on %s dataset ...", names(seu.list)[i]))
      if (do.norm) {
        seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
      }
      seu.list[[i]] <- FindVariableFeatures(object = seu.list[[i]],
                                            selection.method = selection.method,
                                            nfeatures = 2000,
                                            verbose = FALSE,
                                            ...)
    } else {
      message2(sprintf("Skip %s dataset due to too less cells (%s cells).", names(seu.list)[i], ncol(seu.list[[i]])))
    }
  }
  message2("Integrate features ...")
  seu.list <- seu.list[sapply(seu.list, ncol) >= min.cells]
  assay <- DefaultAssay(object)
  features <- SelectIntegrationFeatures(object.list = seu.list,
                                        nfeatures = nfeatures,
                                        assay = rep(assay, length(x = seu.list)))
  message2("Done.")
  return(features)
}
