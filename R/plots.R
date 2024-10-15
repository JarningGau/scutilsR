#' Perform Batch Enrichment Analysis
#'
#' This function performs enrichment analysis on a list of marker sets using the clusterProfiler package.
#'
#' @param all.markers A list of marker sets, where each element is a vector of gene identifiers.
#' @param t2g A data frame with at least two columns: term and gene, for the enrichment analysis.
#'
#' @return A list of enrichment results (clusterProfiler::enricher object), with each element corresponding to an input marker set.
#'          NULL results are removed from the output.
#'
#' @examples
#' # Example usage:
#' # result <- enrich.batch(my_markers, my_term2gene)
#'
#' @export
enrich.batch <- function(all.markers, t2g) {
  e.res <- pbapply::pblapply(all.markers, function(xx) {
    clusterProfiler::enricher(xx, TERM2GENE = t2g)
  })
  names(e.res) <- names(all.markers)
  ## remove NULL
  no.result <- names(e.res)[is.null(e.res)]
  if (length(no.result) > 0) {
    no.result.cont <- paste(no.result, collapse = ",")
    message(glue::glue("{length(no.result)} with no enrichment results: {no.result.cont}."))
    e.res <- e.res[!is.null(e.res)]
  }
  return(e.res)
}

#' Create a Dot Plot for Enrichment Results
#'
#' This function generates a dot plot visualization for enrichment analysis results.
#'
#' @param e.res A list of enrichment results (clusterProfiler::enricher objects).
#' @param topN Integer. The number of top terms to include for each group. Default is 4.
#'
#' @return A ggplot object representing the dot plot of enrichment results.
#'
#' @details The function creates a dot plot where each dot represents an enriched term.
#' The size of the dot indicates the count of genes, and the color represents the
#' adjusted p-value (-log10 transformed).
#'
#' @examples
#' # Assuming 'enrichment_results' is your list of enrichment results
#' # plot <- enrich.dotplot(enrichment_results)
#' # print(plot)
#'
#' @import ggplot2
#' @importFrom stringr str_wrap
#'
#' @export
enrich.dotplot <- function(e.res, topN = 4) {
  data.use <- lapply(names(e.res), function(xx) {
    res <- head(e.res[[xx]]@result, 4)
    res$group <- xx
    res
  }) %>% Reduce(rbind, .)
  term.levels <- unique(data.use$Description)
  data.use <- lapply(names(e.res), function(xx) {
    res <- subset(e.res[[xx]]@result, Description %in% term.levels)
    res$group <- xx
    res
  }) %>% Reduce(rbind, .)
  data.use$Description <- factor(data.use$Description, levels = rev(term.levels))
  data.use$group <- factor(data.use$group, levels = names(e.res))
  ggplot(data.use, aes(group, Description)) +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), shape=21, color = "black") +
    scale_fill_gradientn(colors = c("white","red")) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) +
    labs(x = "") +
    theme_bw(base_size = 10) +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(color = "black"))
}