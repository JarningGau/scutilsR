#' Main function of CellChat analysis
#' @param seu seurat object
#' @param label.field filed name in seurat meta.data contains cell labels.
#' @param name run name, the prefix of final rds file and pdf plot.
#' @param mode one of default, sensitive, population.
#' - default: default CellChat settings.
#' - sensitive: sensitive to low expressed genes.
#' - population: whether consider the proportion of cells in each group across all sequenced cells.
#' @param DB cellchat database. Default: NULL
#' @param out.dir output path, including pdf plots and cellchat rds file. Default: ./
#' @param cores threads used
#' @param fig.width figure width for pdf file. Default: 10
#' @param fig.height figure height for pdf file. Default: 8
#' @import patchwork
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @export
CellChatHelper <- function(seu, label.field, name="name", mode = "default", DB=NULL, out.dir=getwd(), cores=10, fig.width=10, fig.height=8) {
  data.input <- GetAssayData(seu, slot = "data", assay = "RNA")
  meta <- seu[[label.field]]
  if (is.null(DB)) {
    DB <- CellChat::subsetDB(CellChat::CellChatDB.human, search = "Secreted Signaling")
  }
  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = label.field)
  cellchat@DB <- DB

  cellchat <- CellChat::subsetData(cellchat)
  if (cores > 1) {
    future::plan("multiprocess", workers = cores)
  }
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  if (mode == "default") {
    cellchat <- CellChat::computeCommunProb(cellchat)
  } else if (mode == "population") {
    cellchat <- CellChat::computeCommunProb(cellchat, population.size = T)
  } else if (mode == "sensitive") {
    cellchat <- CellChat::computeCommunProb(cellchat, population.size = T, type = "truncatedMean", trim = 0.1)
  } else {
    stop("Error! Please set mode = one of [default, sensitive, population].")
  }
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  saveRDS(cellchat, file = file.path(out.dir, sprintf("cellchat.%s.%s.rds", name, mode)))

  pdf(file.path(out.dir, sprintf("cellchat.%s.%s.pdf", name, mode)), width = 10, height = 8)
  groupSize <- as.numeric(table(cellchat@idents))
  CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  mat <- cellchat@net$weight
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    print(CellChat::netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
  }

  for (i in cellchat@netP$pathways) {
    print(CellChat::netVisual_aggregate(cellchat, signaling = i, layout = "chord"))
  }

  for (i in cellchat@netP$pathways) {
    p <- CellChat::plotGeneExpression(cellchat, signaling = i)
    print(p + plot_annotation(title = paste(i, "pathway")) & theme(plot.title = element_text(hjust = .5)))
  }
  dev.off()
}
