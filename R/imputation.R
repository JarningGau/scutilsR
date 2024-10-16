#' Impute gene expression data using Non-negative Matrix Factorization (NMF)
#'
#' This function performs data imputation on a Seurat object using NMF.
#'
#' @param seu A Seurat object containing the RNA assay with raw counts.
#' @param min_cells Minimum number of cells a gene must be expressed in to be included. Default is 50.
#' @param k Number of factors for NMF. Default is 100.
#' @param threads Number of threads to use for parallel computation. Default is 1.
#' @param seed Random seed for reproducibility. Default is 1024.
#'
#' @return A Seurat object with an additional "imputed" assay containing the imputed data.
#'
#' @details
#' The function performs the following steps:
#' 1. Filters genes based on the minimum number of cells they are expressed in.
#' 2. Extracts the log-normalized expression matrix from the RNA assay.
#' 3. Applies NMF using the RcppML package.
#' 4. Reconstructs the imputed matrix from the NMF factors.
#' 5. Adds the imputed matrix as a new assay to the Seurat object.
#'
#' @examples
#' seu <- qs::qread("output_03/06.HS_BM_donor1.seurat.lineage.qs")
#' seu <- impute_nmf(seu)
#'
#' @import Seurat
#' @import RcppML
#' @import Matrix
#' @importFrom parallel detectCores
#' @export
impute_nmf <- function(seu, min_cells = 50, k = 100, threads = 1, seed = 1024) {
  #### Method-1: NMF ####
  ## [g,c] ~ [g,k] %*% [k,k] %*% [k,c]
  ## [g,c] ~ [g,k] %*% [k,c] # k << c
  expr.in.cells <- Matrix::rowSums(seu[["RNA"]]@counts > 0)
  genes <- names(expr.in.cells)[expr.in.cells >= min_cells]
  expr.mat <- seu[["RNA"]]@data[genes, ] # log normalized matrix
  
  # Check and adjust threads parameter
  max_threads <- parallel::detectCores()
  threads <- min(max_threads %/% 2, threads)
  
  RcppML::setRcppMLthreads(threads)
  model <- RcppML::nmf(expr.mat, k = k, verbose = TRUE, seed = seed)
  imputed.mat <- model$w %*% diag(model$d) %*% model$h
  colnames(imputed.mat) <- colnames(expr.mat)
  rownames(imputed.mat) <- rownames(expr.mat)
  
  seu[["imputed"]] <- CreateAssayObject(data = imputed.mat)
  
  return(seu)
}
