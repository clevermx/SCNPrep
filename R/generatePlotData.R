ALL_REDUCTIONS <- c("FItSNE", "tsne", "umap", "pca")

#' Generate data.frame with plot data
#'
#' @param object SeuratObject
#' @param userAnnotations list of data.frames with user-defined annotation for cells
#' @param maxReductionDims maximum number of dimensions for each dim-reduction technique to take.
#' Mostly used for reducing number of PCs reported in the object.
#'
#' @return
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize_at
#' @importFrom dplyr sym
#' @importFrom Seurat Idents
#' @import jsonlite
#'
#' @examples
generatePlotData <- function(object, userAnnotations, maxReductionDims) {
  presentAssays <- Assays(object)
  reductions <- ALL_REDUCTIONS
  reductions <- reductions[reductions %in% names(object@reductions)]

  if (length(reductions) == 0) {
    stop("No TSNE, UMAP or PCA found, please make sure to run some dimensionality reduction first")
  }

  embeddings <- lapply(reductions, function(red) {
    emb <- object@reductions[[red]]@cell.embeddings
    dimMax <- min(ncol(emb), maxReductionDims)
    reduced <- emb[, 1:dimMax]
    reduced
  })

  dataForPlot <- as.data.frame(do.call(cbind, embeddings))

  if (length(levels(object$orig.ident)) > 1 || length(unique(object$orig.ident)) > 1) {
    dataForPlot$Sample <- object$orig.ident
    dataForPlot$Sample <- as.factor(dataForPlot$Sample)
  } else {
    samples <- gsub("(\\w+)(_|:)[NATGCx]+", "\\1", colnames(object))
    if (length(unique(samples)) > 1 && length(unique(samples)) < 100) {
      dataForPlot$Sample <- as.factor(samples)
    }
  }


  dataForPlot$Cluster <-  Seurat::Idents(object = object)
  metaColumns <- colnames(object@meta.data)
  umiGeneColumns <- c()

  for (assayIter in presentAssays) {
    umiColumn <- sprintf("nCount_%s", assayIter)
    geneColumn <- sprintf("nFeature_%s", assayIter)
    if (umiColumn %in% metaColumns) umiGeneColumns <- c(umiGeneColumns, umiColumn)
    if (geneColumn %in% metaColumns) umiGeneColumns <- c(umiGeneColumns, geneColumn)
  }

  for (column in umiGeneColumns) {
    log2column <- sprintf("%s_log2", column)
    dataForPlot[, column] <- object@meta.data[, column]
    dataForPlot[, log2column] <- log2(object@meta.data[, column])
  }

  clusterColnames <- grep("^(C|c)luster", colnames(dataForPlot), value = T)

  for (clusterColumn in clusterColnames) {
    dataForPlot[, clusterColumn] <- as.factor(dataForPlot[, clusterColumn])
  }

  for (userAnnotation in userAnnotations) {
    dataForPlot <- cbind(dataForPlot, userAnnotation[rownames(dataForPlot), ])
  }

  fields <- list()

  for (column in colnames(dataForPlot)) {
    thisColumn <- dataForPlot[, column]
    if (is.numeric(thisColumn)) {
      fields[[column]] = list(
        type=unbox("numeric"),
        range=c(min(thisColumn), max(thisColumn))
      )
    }
    if (is.factor(thisColumn)) {
      fields[[column]] = list(
        type=unbox("factor"),
        levels=levels(thisColumn)
      )
    }
  }


  annotations <- list()

  for (reduction in reductions) {
    for (clusterCol in clusterColnames) {
      dimColnames <- colnames(object@reductions[[reduction]]@cell.embeddings)[1:2]

      centers <- dataForPlot %>%
        dplyr::group_by(!!dplyr::sym(clusterCol)) %>%
        dplyr::summarize_at(dimColnames, median)
      centers <- as.data.frame(centers)
      centers$Text <- centers[, clusterCol]
      masks <- maskEstimator(dataForPlot, clusterCol, dimColumns = dimColnames)


      centersName = sprintf("%s_%s_centers", reduction, clusterCol)
      bordersName = sprintf("%s_%s_borders", reduction, clusterCol)

      annotations[[centersName]] <- list(
        type=unbox("text"),
        value=unbox(clusterCol),
        coords=dimColnames,
        data=centers
      )

      annotations[[bordersName]] <- list(
        type=unbox("polygon"),
        value=unbox("group"),
        coords=dimColnames,
        data=masks
      )
    }
  }

  plotDataForJson = list(
    fields = fields,
    data = dataForPlot,
    annotations = annotations
  )
  return(plotDataForJson)
}
