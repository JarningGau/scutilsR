#' multi-core FindAllMarkers()
#' @param seu Seurat object
#' @param do.flatten Flatten the marker list into data.frame. Default: TRUE
#' @param only.pos Only positive markers. Default: TRUE
#' @param n.cores The number of threads. Default: 10
#' @param ... Arguments passed to Seurat::FindMarkers()
#' @return data.frame or list contains putative markers.
#' @references https://gist.github.com/diazdc/1735102c243cd16acb1b1f3fd09a26e1
#' @export
mcFindAllMarkers <- function(seu, do.flatten=TRUE, only.pos=TRUE, n.cores=10, ...) {
  n_clust <- levels(Seurat::Idents(seu))
  find.markers <- function(i){
    ident1 <- i
    ident2 <- n_clust[n_clust != i]
    table <- Seurat::FindMarkers(seu, ident.1 = ident1, ident.2 = ident2, only.pos = only.pos, ...)
    table$Gene.name.uniq <- rownames(table)
    table$cluster <- rep(i, nrow(table))
    return(table)
  }

  marker_results <- parallel::mclapply(n_clust, find.markers, mc.cores = n.cores)
  names(marker_results) <- n_clust
  # nice way to flatten list into a single DF
  if (do.flatten) {
    marker_results <- dplyr::bind_rows(marker_results)
  }
  return(marker_results)
}
