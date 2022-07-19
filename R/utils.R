#' Make directory if it not exists.
#' @param path A character vector containing a single path name
#' @export
safe_mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
}

#' Negate operation of \%in\%
#' @usage letters \%notin\% letters[1:4]
#' @export
`%notin%` <- Negate(`%in%`)


#' calculating the position of cluster labels
#' @param data A data.frame contains embeddings and clusters.
#' @param emb A string indicating the prefix of embedding names. Default: "tSNE"
#' @param group.by A string indicating the field of clusters. Default: "ClusterID"
#' @return A data.frame contains median coordinates of each cluster.
#' @export
get_label_pos <- function(data, emb = "tSNE", group.by="ClusterID") {
  new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by)]
  colnames(new.data) <- c("x","y","cluster")
  clusters <- names(table(new.data$cluster))
  new.pos <- lapply(clusters, function(i) {
    tmp.data = new.data[new.data$cluster == i, ]
    data.frame(
      x = stats::median(tmp.data$x),
      y = stats::median(tmp.data$y),
      label = i)
  })
  do.call(rbind, new.pos)
}
